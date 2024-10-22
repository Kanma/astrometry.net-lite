/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <stdlib.h>

#include "codekd.h"
#include "kdtree_fits_io.h"
#include "starutil.h"
#include "errors.h"

static codetree_t* codetree_alloc() {
    codetree_t* s = calloc(1, sizeof(codetree_t));
    if (!s) {
        fprintf(stderr, "Failed to allocate a code kdtree struct.\n");
        return NULL;
    }
    return s;
}

codetree_t* codetree_open_fits(const char* filename, fitsfile* fits) {
    codetree_t* s;
    fits_io_t* io;
    char* treename = CODETREE_NAME;

    s = codetree_alloc();
    if (!s)
        return s;

    io = fits_open_fits(filename, fits);
    if (!io) {
        ERROR("Failed to open FITS file");
        goto bailout;
    }

    if (!kdtree_fits_contains_tree(io, treename))
        treename = NULL;

    s->tree = kdtree_fits_read_tree(io, treename, &s->header);
    if (!s->tree) {
        ERROR("Failed to read code kdtree from file %s\n", filename);
        goto bailout;
    }

    // kdtree_fits_t is a typedef of fitsbin_t
    fits_io_close(io);

    return s;
 bailout:
    free(s);
    return NULL;
}

int codetree_close(codetree_t* s) {
    if (!s) return 0;
    if (s->header)
        free(s->header);
    if (s->tree)
        kdtree_fits_close(s->tree);
    free(s);
    return 0;
}
