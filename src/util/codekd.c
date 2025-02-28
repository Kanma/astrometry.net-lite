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

codetree_t* codetree_open_fits(fits_file_t* fits) {
    codetree_t* s;
    char* treename = CODETREE_NAME;

    s = codetree_alloc();
    if (!s)
        return s;

    if (!kdtree_fits_contains_tree(fits, treename, NULL))
        treename = NULL;

    s->tree = kdtree_fits_read_tree(fits, treename, &s->header);
    if (!s->tree) {
        ERROR("Failed to read code kdtree from file %s\n", fits->filename);
        goto bailout;
    }

    return s;
 bailout:
    free(s);
    return NULL;
}

int codetree_close(codetree_t* s) {
    if (!s) return 0;
    if (s->tree)
        kdtree_fits_close(s->tree);
    free(s);
    return 0;
}

void parse_codetree_params(fitsfile* fits, fits_hdu_t* header) {
    int status = 0;

    header->fits->code.circle = FALSE;
    header->fits->code.cx_less_than_dx = FALSE;
    header->fits->code.meanx_less_than_half = FALSE;

    fits_movabs_hdu(fits, header->extension, NULL, &status);
    if (status != 0)
        return;

    // check for CIRCLE field in ckdt header...
    fits_read_key(fits, TBYTE, "CIRCLE", &header->fits->code.circle, NULL, &status);

    // New indexes are cooked such that cx < dx for all codes, but not
    // all of the old ones are like this.
    status = 0;
    fits_read_key(fits, TBYTE, "CXDX", &header->fits->code.cx_less_than_dx, NULL, &status);

    status = 0;
    fits_read_key(fits, TBYTE, "CXDXLT1", &header->fits->code.meanx_less_than_half, NULL, &status);
}
