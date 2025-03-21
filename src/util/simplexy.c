/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "os-features.h"
#include "simplexy.h"
#include "dimage.h"
#include "simplexy-common.h"
#include "log.h"
#include "errors.h"
#include "resample.h"
#include "an-bool.h"

/*
 * simplexy.c
 *
 * Find sources in a float image.
 *
 * Algorithm outline:
 * 1. Estimate image noise
 * 2. Median filter and subtract to eliminate low-frequence sky gradient
 * 3. Find statistically significant pixels
 *    - Mask those pixels and a box around them
 * 4. Do connected components analysis on the resulting mask to find each
 *    component (source)
 * 5. Find the peaks in each object
 *    For each object:
 *       - Find the objects boundary
 *       - Cut it out
 *       - Smooth it
 *       - Find peaks in resulting cutout
 *       - Chose the most representative peak
 * 6. Extract the flux of each object as the value of the image at the peak
 *
 * BUGS:
 *
 * Note: To make simplexy() reentrant, do the following:
 * #define SIMPLEXY_REENTRANT
 * Or compile all the simplexy files with -DSIMPLEXY_REENTRANT
 *
 * Mike Blanton
 * 1/2006
 *
 */

/** Some testing:

 wget "http://antwrp.gsfc.nasa.gov/apod/image/0910/pleiades_andreo.jpg"
 jpegtopnm pleiades_andreo.jpg | pnmcut 0 0 500 400 | ppmtopgm | pnmtofits > img.fits
 an-fitstopnm -N 0 -X 255 -i img.fits | pnmtopng > img.png
 image2xy -O -p 5 -o img-u8.xy -v -S bgsub.fits -B bg.fits -M mask.fits -U smooth-u8.fits img.fits 
 an-fitstopnm -N 0 -X 255 -i bgsub.fits | pnmtopng > bgsub-u8.png
 an-fitstopnm -N 0 -X 255 -i bg.fits | pnmtopng > bg-u8.png
 an-fitstopnm -N 0 -X 255 -i mask.fits | pnmtopng > mask-u8.png
 an-fitstopnm -N 0 -X 255 -i smooth-u8.fits | pnmtopng > smooth-u8.png
 tabsort -d FLUX img-u8.xy sorted-u8.xy
 pngtopnm img.png | plotxy -i sorted-u8.xy -I - -x 1 -y 1 -N 100 -C red -P | plotxy -i sorted-u8.xy -I - -x 1 -y 1 -n 100 -C red -r 2 -o objs-u8.png

 image2xy -8 -O -p 5 -o img.xy -v -S bgsub.fits -B bg.fits -M mask.fits -U smooth-f.fits img.fits 
 an-fitstopnm -N 0 -X 255 -i bgsub.fits | pnmtopng > bgsub-f.png
 an-fitstopnm -N 0 -X 255 -i bg.fits | pnmtopng > bg-f.png
 an-fitstopnm -N 0 -X 255 -i mask.fits | pnmtopng > mask-f.png
 an-fitstopnm -N 0 -X 255 -i smooth-f.fits | pnmtopng > smooth-f.png

 python <<EOF
 import pyfits
 import sys
 from pylab import *

 figure(figsize=(5,4))
 for infn,outfn in [('smooth-f.fits', 'deltas-f.png'),
 ('smooth-u8.fits', 'deltas-u8.png')]:
 p=pyfits.open(infn)
 img = p[0].data
 (h,w) = img.shape
 sp = 10
 gridsize = 20
 deltas = (img[sp:h:gridsize, sp:w:gridsize]
 - img[0:h-sp:gridsize, 0:w-sp:gridsize]).ravel()
 clf()
 subplot(2,1,1)
 hist(deltas, bins=arange(0, 50, 1))
 title('pixel differences sampled by dsigma()')
 xticks([],[])
 gridsize = 1
 deltas2 = (img[sp:h:gridsize, sp:w:gridsize]
 - img[0:h-sp:gridsize, 0:w-sp:gridsize]).ravel()
 subplot(2,1,2)
 hist(deltas2, bins=arange(0, 50, 1))
 title('all pixel differences')
 savefig(outfn)
 EOF

 */



void simplexy_fill_in_defaults(simplexy_t* s) {
    if (s->dpsf == 0)
        s->dpsf = SIMPLEXY_DEFAULT_DPSF;
    if (s->plim == 0)
        s->plim = SIMPLEXY_DEFAULT_PLIM;
    if (s->dlim == 0)
        s->dlim = SIMPLEXY_DEFAULT_DLIM;
    if (s->saddle == 0)
        s->saddle = SIMPLEXY_DEFAULT_SADDLE;
    if (s->maxper == 0)
        s->maxper = SIMPLEXY_DEFAULT_MAXPER;
    if (s->maxsize == 0)
        s->maxsize = SIMPLEXY_DEFAULT_MAXSIZE;
    if (s->halfbox == 0)
        s->halfbox = SIMPLEXY_DEFAULT_HALFBOX;
    if (s->maxnpeaks == 0)
        s->maxnpeaks = SIMPLEXY_DEFAULT_MAXNPEAKS;
}

void simplexy_free_contents(simplexy_t* s) {
    free(s->image);
    s->image = NULL;
    free(s->image_u8);
    s->image_u8 = NULL;
    free(s->x);
    s->x = NULL;
    free(s->y);
    s->y = NULL;
    free(s->flux);
    s->flux = NULL;
    free(s->background);
    s->background = NULL;
}

int simplexy_run(simplexy_t* s) {
    int i;
    int nx = s->nx;
    int ny = s->ny;
    float limit;
    uint8_t* mask;
    // background-subtracted image.
    float* bgsub = NULL;
    int16_t* bgsub_i16 = NULL;
    // malloc'd background image to free.
    void* bgfree = NULL;
    // PSF-smoothed image.
    float* smoothed = NULL;
    // malloc'd smoothed image to free.
    void* smoothfree = NULL;
    // Connected-components image.
    int* ccimg = NULL;
    int nblobs;
 
    /* Exactly one of s->image and s->image_u8 should be non-NULL.*/
    assert(s->image || s->image_u8);
    assert(!s->image || !s->image_u8);

    logverb("simplexy: nx=%d, ny=%d\n", nx, ny);
    logverb("simplexy: dpsf=%f, plim=%f, dlim=%f, saddle=%f\n",
            s->dpsf, s->plim, s->dlim, s->saddle);
    logverb("simplexy: maxper=%d, maxnpeaks=%d, maxsize=%d, halfbox=%d\n",
            s->maxper, s->maxnpeaks, s->maxsize, s->halfbox);

    if (s->invert) {
        if (s->image) {
            for (i=0; i<nx*ny; i++)
                s->image[i] = -s->image[i];
        } else {
            for (i=0; i<nx*ny; i++)
                s->image_u8[i] = 255 - s->image_u8[i];
        }
    }

    if (s->nobgsub) {
        if (s->image)
            bgsub = s->image;
        else {
            bgsub_i16 = malloc((size_t)nx * (size_t)ny * sizeof(int16_t));
            bgfree = bgsub_i16;
            for (i=0; i<nx*ny; i++)
                bgsub_i16[i] = s->image_u8[i];
        }

    } else {
        // background subtraction via median smoothing.
        logverb("simplexy: median smoothing...\n");

        float* medianfiltered;
        medianfiltered = malloc((size_t)nx * (size_t)ny * sizeof(float));
        bgfree = medianfiltered;
        dmedsmooth(s->image, NULL, nx, ny, s->halfbox, medianfiltered);

        // subtract background from image, placing result in background.
        for (i=0; i<nx*ny; i++)
            medianfiltered[i] = s->image[i] - medianfiltered[i];
        bgsub = medianfiltered;
        medianfiltered = NULL;
    }

    if (s->dpsf > 0.0) {
        smoothed = malloc((size_t)nx * (size_t)ny * sizeof(float));
        smoothfree = smoothed;
        /* smooth by the point spread function (the optimal detection
         filter, since we assume a symmetric Gaussian PSF) */
        if (bgsub)
            dsmooth2(bgsub, nx, ny, s->dpsf, smoothed);
        else
            dsmooth2_i16(bgsub_i16, nx, ny, s->dpsf, smoothed);
    } else {
        if (bgsub)
            smoothed = bgsub;
        else {
            smoothed = malloc((size_t)nx * (size_t)ny * sizeof(float));
            smoothfree = smoothed;
            for (i=0; i<(nx*ny); i++)
                smoothed[i] = bgsub_i16[i];
        }
    }

    // estimate the noise in the image (sigma)
    if (s->sigma == 0.0) {
        logverb("simplexy: measuring image noise (sigma)...\n");
        dsigma(s->image, nx, ny, 5, 0, &(s->sigma));
        logverb("simplexy: found sigma=%g.\n", s->sigma);
    } else {
        logverb("simplexy: assuming sigma=%g.\n", s->sigma);
    }

    /* The noise in the psf-smoothed image is (approximately) 
     *    sigma / (2 * sqrt(pi) * dpsf)
     * This ignores the pixelization, replacing the sum by integral.
     *    The difference is only significant for small sigma, which
     *    would mean your image is undersampled anyway.
     */
    logverb("simplexy: finding objects...\n");
    limit = (s->sigma / (2.0 * sqrt(M_PI) * s->dpsf)) * s->plim;

    if (s->globalbg != 0.0) {
        limit += s->globalbg;
        logverb("Increased detection limit by %g to %g to compensate for global background level\n", s->globalbg, limit);
    }

    /* find pixels above the noise level, and flag a box of pixels around each one. */
    mask = malloc((size_t)nx*(size_t)ny);
    if (!dmask(smoothed, nx, ny, limit, s->dpsf, mask)) {
        FREEVEC(smoothfree);
        FREEVEC(bgfree);
        FREEVEC(mask);
        return 0;
    }
    FREEVEC(smoothfree);

    /* find connected-components in the mask image. */
    ccimg = malloc((size_t)nx * (size_t)ny * sizeof(int));
    dfind2_u8(mask, nx, ny, ccimg, &nblobs);
    FREEVEC(mask);
    logverb("simplexy: found %i blobs\n", nblobs);

    s->x = malloc(s->maxnpeaks * sizeof(float));
    s->y = malloc(s->maxnpeaks * sizeof(float));
	
    /* find all peaks within each object */
    logverb("simplexy: finding peaks...\n");
    if (bgsub)
        dallpeaks(bgsub, nx, ny, ccimg, s->x, s->y, &(s->npeaks), s->dpsf,
                  s->sigma, s->dlim, s->saddle, s->maxper, s->maxnpeaks, s->sigma, s->maxsize);
    else
        dallpeaks_i16(bgsub_i16, nx, ny, ccimg, s->x, s->y, &(s->npeaks), s->dpsf,
                      s->sigma, s->dlim, s->saddle, s->maxper, s->maxnpeaks, s->sigma, s->maxsize);
    logmsg("simplexy: found %i sources.\n", s->npeaks);
    FREEVEC(ccimg);

    s->x   = realloc(s->x, s->npeaks * sizeof(float));
    s->y   = realloc(s->y, s->npeaks * sizeof(float));
    s->flux       = malloc(s->npeaks * sizeof(float));
    s->background = malloc(s->npeaks * sizeof(float));

    for (i = 0; i < s->npeaks; i++) {
        // round
        int ix = (int)(s->x[i] + 0.5);
        int iy = (int)(s->y[i] + 0.5);
        Unused anbool finite;
        finite = isfinite(s->x[i]);
        assert(finite);
        finite = isfinite(s->y[i]);
        assert(finite);
        // these coordinates are now 0,0 is center of first pixel.
        assert(ix >= 0);
        assert(iy >= 0);
        assert(ix < nx);
        assert(iy < ny);
        if (bgsub) {
            s->flux[i]       = bgsub[ix + iy * nx];
            s->background[i] = s->image[ix + iy * nx] - s->flux[i];
        } else {
            s->flux[i]       = bgsub_i16[ix + iy * nx];
            s->background[i] = (float)s->image_u8[ix + iy * nx] - s->flux[i];
        }

        s->flux[i] -= s->globalbg;
        s->background[i] += s->globalbg;
    }

    FREEVEC(bgfree);

    return 1;
}

void simplexy_clean_cache() {
    dselip_cleanup();
}

