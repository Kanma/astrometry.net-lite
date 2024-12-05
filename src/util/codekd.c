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

codetree_t* codetree_open_fits(const char* filename, fits_file_t* fits) {
    codetree_t* s;
    char* treename = CODETREE_NAME;

    s = codetree_alloc();
    if (!s)
        return s;

    if (!kdtree_fits_contains_tree(fits, treename))
        treename = NULL;

    s->tree = kdtree_fits_read_tree(fits, treename, &s->header);
    if (!s->tree) {
        ERROR("Failed to read code kdtree from file %s\n", filename);
        goto bailout;
    }

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
