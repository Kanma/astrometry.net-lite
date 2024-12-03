/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "quadfile.h"
#include "starutil.h"
#include "ioutils.h"
#include "errors.h"
#include "an-endian.h"

#define CHUNK_QUADS 0

static quadfile_t* new_quadfile(const char* fn, fitsfile* fits) {
    quadfile_t* qf;
    qf = calloc(1, sizeof(quadfile_t));
    if (!qf) {
        SYSERROR("Couldn't malloc a quadfile struct");
        return NULL;
    }
    qf->healpix = -1;
    qf->hpnside = 1;

    qf->io = fits_open_fits(fn, fits);
    if (!qf->io) {
        ERROR("Failed to create fitsbin");
        return NULL;
    }


    qf->dimquads = 4;
    qf->numquads = -1;
    qf->numstars = -1;
    qf->index_scale_upper = -1.0;
    qf->index_scale_lower = -1.0;
    qf->indexid = 0;
    qf->healpix = -1;
    qf->hpnside = 1;

    int status = 0;
    fits_movabs_hdu(qf->io->fits, 1, NULL, &status);

    if (status == 0)
    {
        fits_read_key(qf->io->fits, TINT, "DIMQUADS", &qf->dimquads, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TINT, "NQUADS", &qf->numquads, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TINT, "NSTARS", &qf->numstars, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TDOUBLE, "SCALE_U", &qf->index_scale_upper, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TDOUBLE, "SCALE_L", &qf->index_scale_lower, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TINT, "INDEXID", &qf->indexid, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TINT, "HEALPIX", &qf->healpix, NULL, &status);

        status = 0;
        fits_read_key(qf->io->fits, TINT, "HPNSIDE", &qf->hpnside, NULL, &status);
    }

    if ((qf->numquads == -1) || (qf->numstars == -1) ||
        (qf->index_scale_upper == -1.0) || (qf->index_scale_lower == -1.0)) {
        ERROR("Couldn't find NQUADS or NSTARS or SCALE_U or SCALE_L entries in FITS header");
        {
            quadfile_close(qf);
            return NULL;
        }
    }

    fits_hdu_t primheader;
    primheader.fits = qf->io->fits;
    primheader.extension = 1;

    if (fits_check_endian(&primheader)) {
        ERROR("Quad file was written with the wrong endianness");
        {
            quadfile_close(qf);
            return NULL;
        }
    }

    return qf;
}

int quadfile_check(const quadfile_t* qf) {
    int q, i;
    if (qf->dimquads < 3 || qf->dimquads > DQMAX) {
        ERROR("Dimquads has illegal value %i", qf->dimquads);
        return -1;
    }
    for (q=0; q<qf->numquads; q++) {
        unsigned int stars[DQMAX];
        if (quadfile_get_stars(qf, q, stars)) {
            ERROR("Failed to get quad %i of %i", q, qf->numquads);
            return -1;
        }
        for (i=0; i<qf->dimquads; i++) {
            if (stars[i] >= qf->numstars) {
                ERROR("Star ID %i is out of bounds: num stars %i", stars[i], qf->numstars);
                return -1;
            }
        }
    }
    return 0;
}

int quadfile_dimquads(const quadfile_t* qf) {
    return qf->dimquads;
}

int quadfile_nquads(const quadfile_t* qf) {
    return qf->numquads;
}

fits_hdu_t* quadfile_get_header(const quadfile_t* qf) {
    return fits_get_primary_header(qf->io);
}

static quadfile_t* my_open(const char* fn, fitsfile* fits) {
    quadfile_t* qf = NULL;
 
    qf = new_quadfile(fn, fits);
    if (!qf)
        goto bailout;

    fits_read_chunk(qf->io, "quads", qf->dimquads * sizeof(uint32_t), &qf->numquads, &qf->quadarray);

    return qf;

 bailout:
    if (qf)
        quadfile_close(qf);
    return NULL;
}

char* quadfile_get_filename(const quadfile_t* qf) {
    return qf->io->filename;
}

quadfile_t* quadfile_open_fits(const char* filename, fitsfile* fits) {
    return my_open(filename, fits);
}

int quadfile_close(quadfile_t* qf) {
    int rtn;
    if (!qf) return 0;
    rtn = fits_io_close(qf->io);
    free(qf->quadarray);
    free(qf);
    return rtn;
}

double quadfile_get_index_scale_upper_arcsec(const quadfile_t* qf) {
    return rad2arcsec(qf->index_scale_upper);
}

double quadfile_get_index_scale_lower_arcsec(const quadfile_t* qf) {
    return rad2arcsec(qf->index_scale_lower);
}

int quadfile_get_stars(const quadfile_t* qf, unsigned int quadid, unsigned int* stars) {
    int i;
    if (quadid >= qf->numquads) {
        ERROR("Requested quadid %i, but number of quads is %i",	quadid, qf->numquads);
        assert(quadid < qf->numquads);
        return -1;
    }

    for (i=0; i<qf->dimquads; i++) {
        stars[i] = qf->quadarray[quadid * qf->dimquads + i];
    }
    return 0;
}



