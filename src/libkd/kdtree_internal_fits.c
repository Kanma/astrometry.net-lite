/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef KDTREE_NO_FITS

#include "kdtree_fits_io.h"
#include "kdtree.h"
#include "ioutils.h"
#include "errors.h"

#define KDTREE_FITS_VERSION 1

static char* get_table_name(const char* treename, const char* tabname) {
    char* rtn;
    if (!treename) {
        return strdup_safe(tabname);
    }
    asprintf_safe(&rtn, "%s_%s", tabname, treename);
    return rtn;
}

int MANGLE(kdtree_read_fits)(fits_file_t* io, kdtree_t* kd) {
    char* tablename;

    // kd->lr
    tablename = get_table_name(kd->name, KD_STR_LR);
    fits_read_chunk(io, tablename, sizeof(u32), &kd->nbottom, &kd->lr, 0);
    free(tablename);

    // kd->perm
    tablename = get_table_name(kd->name, KD_STR_PERM);
    fits_read_chunk(io, tablename, sizeof(u32), &kd->ndata, &kd->perm, 0);
    free(tablename);

    // kd->bb
    kd->n_bb = 0;
    tablename = get_table_name(kd->name, KD_STR_BB);
    fits_read_chunk(io, tablename, sizeof(ttype) * kd->ndim * 2,  &kd->n_bb, &kd->bb.any, 0);
    free(tablename);

    // kd->split
    tablename = get_table_name(kd->name, KD_STR_SPLIT);
    fits_read_chunk(io, tablename, sizeof(ttype), &kd->ninterior, &kd->split.any, 0);
    free(tablename);

    // kd->splitdim
    tablename = get_table_name(kd->name, KD_STR_SPLITDIM);
    fits_read_chunk(io, tablename, sizeof(u8), &kd->ninterior, &kd->splitdim, 0);
    free(tablename);

    // kd->data
    tablename = get_table_name(kd->name, KD_STR_DATA);
    fits_read_chunk(io, tablename, sizeof(dtype) * kd->ndim, &kd->ndata, &kd->data.any, 0);
    free(tablename);

    // kd->minval/kd->maxval/kd->scale
    double* r;
    int nb = (kd->ndim * 2 + 1);
    tablename = get_table_name(kd->name, KD_STR_RANGE);
    if (fits_read_chunk(io, tablename, sizeof(double), &nb, &r, 1) == 0) {
        kd->minval = r;
        kd->maxval = r + kd->ndim;
        kd->scale  = r[kd->ndim * 2];
        kd->invscale = 1.0 / kd->scale;
    }
    free(tablename);

    if (!(kd->bb.any ||
          (kd->split.any && (TTYPE_INTEGER || kd->splitdim)))) {
        ERROR("kdtree contains neither bounding boxes nor split+dim data");
        return -1;
    }

    if ((TTYPE_INTEGER && !ETYPE_INTEGER) &&
        !(kd->minval && kd->maxval)) {
        ERROR("treee does not contain required range information");
        return -1;
    }

    if (kd->split.any) {
        if (kd->splitdim)
            kd->splitmask = UINT32_MAX;
        else
            compute_splitbits(kd);
    }

    return 0;
}

#endif
