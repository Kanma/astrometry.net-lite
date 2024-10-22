/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef KDTREE_FITS_IO_H
#define KDTREE_FITS_IO_H

#include <stdio.h>

#include "astrometry/kdtree.h"
#include "astrometry/fits_io.h"


int kdtree_fits_contains_tree(const fits_io_t* io, const char* treename);

kdtree_t* kdtree_fits_read_tree(fits_io_t* io, const char* treename, fits_hdu_t** p_hdr);

int kdtree_fits_close(kdtree_t* io);


// names (actually prefixes) of FITS tables.
#define KD_STR_HEADER    "kdtree_header"
#define KD_STR_LR        "kdtree_lr"
#define KD_STR_PERM      "kdtree_perm"
#define KD_STR_BB        "kdtree_bb"
#define KD_STR_SPLIT     "kdtree_split"
#define KD_STR_SPLITDIM  "kdtree_splitdim"
#define KD_STR_DATA      "kdtree_data"
#define KD_STR_RANGE     "kdtree_range"

// is the given column name one of the above strings?
int kdtree_fits_column_is_kdtree(char* columnname);

#endif
