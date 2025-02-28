/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#include "kdtree_fits_io.h"
#include "kdtree_internal.h"
#include "kdtree_mem.h"
#include "ioutils.h"
#include "errors.h"
#include "tic.h"
#include "log.h"


static void parse_tree_header(fitsfile* fits, fits_hdu_t* header, int oldstyle) {
    unsigned int ext_type, int_type, data_type;
    char str[FITS_LINESZ+1];
    int status = 0;

    if (header->hdutype != BINARY_TBL)
        return;

    if (oldstyle) {
        fits_read_key(fits, TINT, "NDIM", &header->tree.ndim, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "NDATA", &header->tree.ndata, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "NNODES", &header->tree.nnodes, NULL, &status);
    } else {
        fits_read_key(fits, TINT, "KDT_NDIM", &header->tree.ndim, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "KDT_NDAT", &header->tree.ndata, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "KDT_NNOD", &header->tree.nnodes, NULL, &status);
    }

    if (status != 0)
        return;

    fits_read_key(fits, TSTRING, "KDT_EXT", str, NULL, &status);
    ext_type = kdtree_kdtype_parse_ext_string(str);

    status = 0;
    fits_read_key(fits, TSTRING, "KDT_INT", str, NULL, &status);
    int_type = kdtree_kdtype_parse_tree_string(str);

    status = 0;
    fits_read_key(fits, TSTRING, "KDT_DATA", str, NULL, &status);
    data_type = kdtree_kdtype_parse_data_string(str);

    // default: external world is doubles.
    if (ext_type == KDT_NULL)
        ext_type = KDT_EXT_DOUBLE;

    header->tree.treetype = kdtree_kdtypes_to_treetype(ext_type, int_type, data_type);

    header->tree.has_linear_lr = 0;
    fits_read_key(fits, TBYTE, "KDT_LINL", &header->tree.has_linear_lr, NULL, &status);
}

static int is_tree_header_ok(fits_hdu_t* header) {
    if ((header->tree.ndim > -1) && (header->tree.ndata > -1) && (header->tree.nnodes > -1) &&
        (header->tree.treetype != KDT_NULL) &&
        (fits_check_endian(header) == 0)) {
        return 1;
    }
    return 0;
}

// declarations
KD_DECLARE(kdtree_read_fits, int, (fits_file_t* io, kdtree_t* kd));

static fits_hdu_t* find_tree(const char* treename, const fits_file_t* io, char** realname) {
    fits_hdu_t* header;
    int i;
    int status = 0;

    if (!treename) {
        // Look in the primary header...
        header = fits_get_primary_header(io);
        if (is_tree_header_ok(header)) {
            return header;
        }
    }

    // treat empty treename as NULL...
    if (treename && !treename[0])
        treename = NULL;

    // scan the extension headers, looking for one that contains a matching KDT_NAME entry.
    header = NULL;
    for (i = 1; i < io->nbHDUs; ++i) {
        status = 0;
        header = &io->hdus[i];

        // if the desired treename was specified, it must match;
        // bail out if this is not the case.
        // (if treename is NULL then anything matches.)
        if (treename && !((header->tree.name[0] != 0) && (strcmp(header->tree.name, treename) == 0))) {
            continue;
        }

        if (is_tree_header_ok(header)) {
            if (realname)
            {
                *realname = malloc(strlen(header->tree.name) + 1);
                strcpy(*realname, header->tree.name);
            }

            return header;
        }
    }

    return NULL;
}

void kdtree_parse(fitsfile* fits, fits_file_t* io) {
    int status = 0;

    fits_movabs_hdu(fits, 1, NULL, &status);
    parse_tree_header(fits, &io->hdus[0], 1);

    for (int i = 1; i < io->nbHDUs; ++i) {
        fits_hdu_t* header = &io->hdus[i];

        status = 0;
        fits_movabs_hdu(fits, i + 1, NULL, &status);
        fits_read_key(fits, TSTRING, "KDT_NAME", header->tree.name, NULL, &status);
        if (status != 0)
            continue;

        parse_tree_header(fits, header, 0);
    }
}

int kdtree_fits_contains_tree(const fits_file_t* io, const char* treename, fits_hdu_t** p_hdr) {
    fits_hdu_t* hdr = find_tree(treename, io, NULL);
    if (p_hdr)
        *p_hdr = hdr;
    return (hdr != NULL);
}

kdtree_t* kdtree_fits_read_tree(fits_file_t* io, const char* treename,
                                fits_hdu_t** p_hdr) {
    kdtree_t* kd = NULL;
    fits_hdu_t* header;
    int rtn = 0;
    int status = 0;

    kd = CALLOC(1, sizeof(kdtree_t));
    if (!kd) {
        SYSERROR("Couldn't allocate kdtree");
        return NULL;
    }

    header = find_tree(treename, io, &kd->name);
    if (!header) {
        // Not found.
        if (treename)
            ERROR("Kdtree header for a tree named \"%s\" was not found", treename);
        else
            ERROR("Kdtree header was not found");

        FREE(kd);
        return NULL;
    }

    kd->has_linear_lr = header->tree.has_linear_lr;

    if (p_hdr)
        *p_hdr = header;

    kd->ndata  = header->tree.ndata;
    kd->ndim   = header->tree.ndim;
    kd->nnodes = header->tree.nnodes;
    kd->nbottom = (header->tree.nnodes+1)/2;
    kd->ninterior = header->tree.nnodes - kd->nbottom;
    kd->nlevels = kdtree_nnodes_to_nlevels(header->tree.nnodes);
    kd->treetype = header->tree.treetype;

    KD_DISPATCH(kdtree_read_fits, kd->treetype, rtn = , (io, kd));

    if (rtn) {
        FREE(kd->name);
        FREE(kd);
        return NULL;
    }

    kdtree_update_funcs(kd);

    return kd;
}

int kdtree_fits_close(kdtree_t* kd) {
    if (!kd) return 0;
    FREE(kd->name);
    FREE(kd);
    return 0;
}
