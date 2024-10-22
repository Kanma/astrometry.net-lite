/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "starkd.h"
#include "kdtree.h"
#include "starutil.h"
#include "errors.h"
#include "tic.h"
#include "log.h"
#include "ioutils.h"


static startree_t* startree_alloc() {
    startree_t* s = calloc(1, sizeof(startree_t));
    if (!s) {
        fprintf(stderr, "Failed to allocate a star kdtree struct.\n");
        return NULL;
    }
    return s;
}


void startree_search_for_radec(const startree_t* s, double ra, double dec, double radius,
                               double** xyzresults, double** radecresults,
                               int** starinds, int* nresults) {
    double xyz[3];
    double r2;
    radecdeg2xyzarr(ra, dec, xyz);
    r2 = deg2distsq(radius);
    startree_search_for(s, xyz, r2, xyzresults, radecresults, starinds, nresults);
}

void startree_search_for(const startree_t* s, const double* xyzcenter, double radius2,
                         double** xyzresults, double** radecresults,
                         int** starinds, int* nresults) {
    kdtree_qres_t* res = NULL;
    int opts;
    double* xyz;
    int i, N;

    opts = KD_OPTIONS_SMALL_RADIUS;
    if (xyzresults || radecresults)
        opts |= KD_OPTIONS_RETURN_POINTS;

    res = kdtree_rangesearch_options(s->tree, xyzcenter, radius2, opts);
	
    if (!res || !res->nres) {
        if (xyzresults)
            *xyzresults = NULL;
        if (radecresults)
            *radecresults = NULL;
        if (starinds)
            *starinds = NULL;
        *nresults = 0;
        if (res)
            kdtree_free_query(res);
        return;
    }

    xyz = res->results.d;
    N = res->nres;
    *nresults = N;

    if (radecresults) {
        *radecresults = malloc(N * 2 * sizeof(double));
        for (i=0; i<N; i++)
            xyzarr2radecdegarr(xyz + i*3, (*radecresults) + i*2);
    }
    if (xyzresults) {
        // Steal the results array.
        *xyzresults = xyz;
        res->results.d = NULL;
    }
    if (starinds) {
        *starinds = malloc(res->nres * sizeof(int));
        for (i=0; i<N; i++)
            (*starinds)[i] = res->inds[i];
    }
    kdtree_free_query(res);
}


void startree_search(const startree_t* s, const double* xyzcenter, double radius2,
                     double** xyzresults, double** radecresults, int* nresults) {
    startree_search_for(s, xyzcenter, radius2, xyzresults, radecresults, NULL, nresults);
}

int startree_N(const startree_t* s) {
    return s->tree->ndata;
}

int startree_nodes(const startree_t* s) {
    return s->tree->nnodes;
}

int startree_D(const startree_t* s) {
    return s->tree->ndim;
}

fits_hdu_t* startree_header(const startree_t* s) {
    return s->header;
}

startree_t* startree_open_fits(const char* filename, fitsfile* fits) {
    struct timeval tv1, tv2;
    startree_t* s;
    bl* chunks;
    int i;
    fits_io_t* io;
    char* treename = STARTREE_NAME;

    s = startree_alloc();
    if (!s)
        return NULL;

    gettimeofday(&tv1, NULL);

    io = fits_open_fits(filename, fits);

    gettimeofday(&tv2, NULL);
    debug("kdtree_fits_open() took %g ms\n", millis_between(&tv1, &tv2));
    if (!io) {
        ERROR("Failed to open FITS file \"%s\"", filename);
        goto bailout;
    }

    gettimeofday(&tv1, NULL);
    if (!kdtree_fits_contains_tree(io, treename))
        treename = NULL;
    gettimeofday(&tv2, NULL);
    debug("kdtree_fits_contains_tree() took %g ms\n", millis_between(&tv1, &tv2));

    gettimeofday(&tv1, NULL);
    s->tree = kdtree_fits_read_tree(io, treename, &s->header);
    gettimeofday(&tv2, NULL);
    debug("kdtree_fits_read_tree() took %g ms\n", millis_between(&tv1, &tv2));
    if (!s->tree) {
        ERROR("Failed to read kdtree from file \"%s\"", filename);
        goto bailout;
    }

    // Check the tree dimensionality.
    // (because code trees can be confused...)
    if (s->tree->ndim != 3) {
        logverb("File %s contains a kd-tree with dim %i (not 3), named %s\n",
                filename, s->tree->ndim, treename);
        s->tree->io = NULL;
        goto bailout;
    }

    gettimeofday(&tv1, NULL);
    fits_read_chunk(io, "sweep", sizeof(uint8_t), &s->tree->ndata, &s->sweep);
    gettimeofday(&tv2, NULL);
    debug("reading chunks took %g ms\n", millis_between(&tv1, &tv2));

    // kdtree_fits_t is a typedef of fitsbin_t
    fits_io_close(io);

    return s;

 bailout:
    fits_io_close(io);
    startree_close(s);
    return NULL;
}

/*
 uint64_t startree_get_starid(const startree_t* s, int ind) {
 if (!s->starids)
 return 0;
 return s->starids[ind];
 }
 */
int startree_close(startree_t* s) {
    if (!s) return 0;
    if (s->inverse_perm)
        free(s->inverse_perm);
    if (s->header)
        free(s->header);
    if (s->tree) {
        kdtree_fits_close(s->tree);
    }
    free(s);
    return 0;
}

static int Ndata(const startree_t* s) {
    return s->tree->ndata;
}

int startree_check_inverse_perm(startree_t* s) {
    // ensure that each value appears exactly once.
    int i, N;
    uint8_t* counts;
    N = Ndata(s);
    counts = calloc(Ndata(s), sizeof(uint8_t));
    for (i=0; i<N; i++) {
        assert(s->inverse_perm[i] >= 0);
        assert(s->inverse_perm[i] < N);
        counts[s->inverse_perm[i]]++;
    }
    for (i=0; i<N; i++) {
        assert(counts[i] == 1);
    }
    return 0;
}

void startree_compute_inverse_perm(startree_t* s) {
    if (s->inverse_perm)
        return;
    // compute inverse permutation vector.
    s->inverse_perm = malloc(Ndata(s) * sizeof(int));
    if (!s->inverse_perm) {
        fprintf(stderr, "Failed to allocate star kdtree inverse permutation vector.\n");
        return;
    }
#ifndef NDEBUG
    {
        int i;
        for (i=0; i<Ndata(s); i++)
            s->inverse_perm[i] = -1;
    }
#endif
    kdtree_inverse_permutation(s->tree, s->inverse_perm);
#ifndef NDEBUG
    {
        int i;
        for (i=0; i<Ndata(s); i++)
            assert(s->inverse_perm[i] != -1);
    }
#endif
}

int startree_get_cut_nside(const startree_t* s) {
    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return -1;

    int res = -1;
    fits_read_key(s->header->fits, TINT, "CUTNSIDE", &res, NULL, &status);

    return res;
}

int startree_get_cut_nsweeps(const startree_t* s) {
    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return -1;

    int res = -1;
    fits_read_key(s->header->fits, TINT, "CUTNSWEP", &res, NULL, &status);

    return res;
}

double startree_get_cut_dedup(const startree_t* s) {
    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return -1;

    double res = 0.0;
    fits_read_key(s->header->fits, TDOUBLE, "CUTDEDUP", &res, NULL, &status);

    return res;
}

char* startree_get_cut_band(const startree_t* s) {
    static char* bands[] = { "R", "B", "J" };
    int i;

    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return NULL;

    char str[FITS_LINESZ + 1];
    fits_read_key(s->header->fits, TSTRING, "CUTBAND", str, NULL, &status);
    if (status != 0)
        return NULL;

    for (i=0; i<sizeof(bands) / sizeof(char*); i++) {
        if (streq(str, bands[i])) {
            return bands[i];
        }
    }

    return NULL;
}

int startree_get_cut_margin(const startree_t* s) {
    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return -1;

    int res = -1;
    fits_read_key(s->header->fits, TINT, "CUTMARG", &res, NULL, &status);

    return res;
}

double startree_get_jitter(const startree_t* s) {
    int status = 0;

    fits_movabs_hdu(s->header->fits, s->header->extension, NULL, &status);
    if (status != 0)
        return -1;

    double res = 0.0;
    fits_read_key(s->header->fits, TDOUBLE, "JITTER", &res, NULL, &status);

    return res;
}

int startree_get_sweep(const startree_t* s, int ind) {
    if (ind < 0 || ind >= Ndata(s) || !s->sweep)
        return -1;
    return s->sweep[ind];
}

int startree_get(startree_t* s, int starid, double* posn) {
    if (s->tree->perm && !s->inverse_perm) {
        startree_compute_inverse_perm(s);
        if (!s->inverse_perm)
            return -1;
    }
    if (starid >= Ndata(s)) {
        fprintf(stderr, "Invalid star ID: %u >= %u.\n", starid, Ndata(s));
        assert(0);
        return -1;
    }
    if (s->inverse_perm) {
        kdtree_copy_data_double(s->tree, s->inverse_perm[starid], 1, posn);
    } else {
        kdtree_copy_data_double(s->tree, starid, 1, posn);
    }
    return 0;
}

int startree_get_radec(startree_t* s, int starid, double* ra, double* dec) {
    double xyz[3];
    int rtn;
    rtn = startree_get(s, starid, xyz);
    if (rtn)
        return rtn;
    xyzarr2radecdeg(xyz, ra, dec);
    return rtn;
}
