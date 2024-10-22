/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef CODE_KD_H
#define CODE_KD_H

#include "astrometry/kdtree.h"
#include "astrometry/kdtree_fits_io.h"

#define AN_FILETYPE_CODETREE "CKDT"

#define CODETREE_NAME "codes"

typedef struct {
    kdtree_t* tree;
    fits_hdu_t* header;
} codetree_t;

codetree_t* codetree_open_fits(const char* filename, fitsfile* fits);

int codetree_close(codetree_t* s);

#endif
