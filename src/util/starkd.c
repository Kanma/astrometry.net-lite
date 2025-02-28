/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

startree_t* startree_open_fits(fits_file_t* fits) {
    startree_t* s;
    bl* chunks;
    int i;
    char* treename = STARTREE_NAME;

    s = startree_alloc();
    if (!s)
        return NULL;

    if (!kdtree_fits_contains_tree(fits, treename, NULL))
        treename = NULL;

    s->tree = kdtree_fits_read_tree(fits, treename, &s->header);
    if (!s->tree) {
        ERROR("Failed to read kdtree from file \"%s\"", fits->filename);
        goto bailout;
    }

    // Check the tree dimensionality.
    // (because code trees can be confused...)
    if (s->tree->ndim != 3) {
        logverb("File %s contains a kd-tree with dim %i (not 3), named %s\n",
                fits->filename, s->tree->ndim, treename);
        goto bailout;
    }

    fits_read_chunk(fits, "sweep", sizeof(uint8_t), &s->tree->ndata, &s->sweep, 1);

    return s;

 bailout:
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
    if (s->tree)
        kdtree_fits_close(s->tree);
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

void parse_startree_params(fitsfile* fits, fits_hdu_t* header) {
    int status = 0;

    fits_movabs_hdu(fits, header->extension, NULL, &status);
    if (status != 0)
        return;

    fits_read_key(fits, TINT, "CUTNSIDE", &header->fits->stars.cut_nside, NULL, &status);

    status = 0;
    fits_read_key(fits, TINT, "CUTNSWEP", &header->fits->stars.cut_nsweeps, NULL, &status);

    status = 0;
    fits_read_key(fits, TDOUBLE, "CUTDEDUP", &header->fits->stars.cut_dedup, NULL, &status);

    status = 0;

    static char* bands[] = { "R", "B", "J" };
    int i;
    char str[FITS_LINESZ + 1];

    fits_read_key(fits, TSTRING, "CUTBAND", str, NULL, &status);

    if (status == 0) {
        for (i=0; i<sizeof(bands) / sizeof(char*); i++) {
            if (strstr(str, bands[i]) == str) {
                header->fits->stars.cut_band = bands[i];
            }
        }
    }

    status = 0;
    fits_read_key(fits, TINT, "CUTMARG", &header->fits->stars.cut_margin, NULL, &status);

    status = 0;
    fits_read_key(fits, TDOUBLE, "JITTER", &header->fits->stars.jitter, NULL, &status);
}

int startree_get_cut_nside(const startree_t* s) {
    return s->header->fits->stars.cut_nside;
}

int startree_get_cut_nsweeps(const startree_t* s) {
    return s->header->fits->stars.cut_nsweeps;
}

double startree_get_cut_dedup(const startree_t* s) {
    return s->header->fits->stars.cut_dedup;
}

char* startree_get_cut_band(const startree_t* s) {
    return s->header->fits->stars.cut_band;
}

int startree_get_cut_margin(const startree_t* s) {
    return s->header->fits->stars.cut_margin;
}

double startree_get_jitter(const startree_t* s) {
    return s->header->fits->stars.jitter;
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
