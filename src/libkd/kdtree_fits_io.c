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


static int is_tree_header_ok(fits_hdu_t* header, int* ndim, int* ndata,
                             int* nnodes, unsigned int* treetype, int oldstyle) {
    unsigned int ext_type, int_type, data_type;
    char str[FITS_LINESZ+1];
    int status = 0;

    int hdutype;
    fits_movabs_hdu(header->fits, header->extension, &hdutype, &status);
    if ((status != 0) || (hdutype != BINARY_TBL))
        return 0;

    if (oldstyle) {
        fits_read_key(header->fits, TINT, "NDIM", ndim, NULL, &status);

        status = 0;
        fits_read_key(header->fits, TINT, "NDATA", ndata, NULL, &status);

        status = 0;
        fits_read_key(header->fits, TINT, "NNODES", nnodes, NULL, &status);
    } else {
        fits_read_key(header->fits, TINT, "KDT_NDIM", ndim, NULL, &status);

        status = 0;
        fits_read_key(header->fits, TINT, "KDT_NDAT", ndata, NULL, &status);

        status = 0;
        fits_read_key(header->fits, TINT, "KDT_NNOD", nnodes, NULL, &status);
    }

    if (status != 0)
        return 0;

    fits_read_key(header->fits, TSTRING, "KDT_EXT", str, NULL, &status);
    ext_type = kdtree_kdtype_parse_ext_string(str);

    status = 0;
    fits_read_key(header->fits, TSTRING, "KDT_INT", str, NULL, &status);
    int_type = kdtree_kdtype_parse_tree_string(str);

    status = 0;
    fits_read_key(header->fits, TSTRING, "KDT_DATA", str, NULL, &status);
    data_type = kdtree_kdtype_parse_data_string(str);

    // default: external world is doubles.
    if (ext_type == KDT_NULL)
        ext_type = KDT_EXT_DOUBLE;

    *treetype = kdtree_kdtypes_to_treetype(ext_type, int_type, data_type);

    if ((*ndim > -1) && (*ndata > -1) && (*nnodes > -1) &&
        (int_type != KDT_NULL) && (data_type != KDT_NULL) &&
        (fits_check_endian(header) == 0)) {
        return 1;
    }
    return 0;
}

// declarations
KD_DECLARE(kdtree_read_fits, int, (fits_file_t* io, kdtree_t* kd));

static fits_hdu_t* find_tree(const char* treename, const fits_file_t* io,
                               int* ndim, int* ndata, int* nnodes,
                               unsigned int* tt, char** realname) {
    fits_hdu_t* header;
    int i;
    int status = 0;
    char str[FITS_LINESZ+1];

    if (!treename) {
        // Look in the primary header...
        header = fits_get_primary_header(io);
        if (is_tree_header_ok(header, ndim, ndata, nnodes, tt, 1)) {
            return header;
        }
    }

    // treat empty treename as NULL...
    if (treename && !treename[0])
        treename = NULL;

    // scan the extension headers, looking for one that contains a matching KDT_NAME entry.
    header = NULL;
    for (i = 2; i <= io->nbHDUs; ++i) {
        status = 0;

        fits_movabs_hdu(io->fits, i, NULL, &status);
        fits_read_key(io->fits, TSTRING, "KDT_NAME", str, NULL, &status);
        if (status != 0)
            continue;

	    // if the desired treename was specified, it must match;
	    // bail out if this is not the case.
	    // (if treename is NULL then anything matches.)
    	if (treename && !((str[0] != 0) && (strcmp(str, treename) == 0))) {
            continue;
        }

        header = malloc(sizeof(fits_hdu_t));
        header->fits = io->fits;
        header->extension = i;

        if (is_tree_header_ok(header, ndim, ndata, nnodes, tt, 0)) {
            if (realname)
            {
                *realname = malloc(strlen(str) + 1);
                strcpy(*realname, str);
            }

            return header;
        }

        free(header);
    }

    return NULL;
}

int kdtree_fits_contains_tree(const fits_file_t* io, const char* treename) {
    int ndim, ndata, nnodes;
    unsigned int tt;

    fits_hdu_t* hdr;
    int rtn;

    hdr = find_tree(treename, io, &ndim, &ndata, &nnodes, &tt, NULL);

    rtn = (hdr != NULL);
    if (hdr != NULL)
        free(hdr);

    return rtn;
}

kdtree_t* kdtree_fits_read_tree(fits_file_t* io, const char* treename,
                                fits_hdu_t** p_hdr) {
    int ndim, ndata, nnodes;
    unsigned int tt;
    kdtree_t* kd = NULL;
    fits_hdu_t* header;
    int rtn = 0;
    int status = 0;

    kd = CALLOC(1, sizeof(kdtree_t));
    if (!kd) {
        SYSERROR("Couldn't allocate kdtree");
        return NULL;
    }

    header = find_tree(treename, io, &ndim, &ndata, &nnodes, &tt, &kd->name);
    if (!header) {
        // Not found.
        if (treename)
            ERROR("Kdtree header for a tree named \"%s\" was not found in file %s", treename, io->fits->Fptr->filename);
        else
            ERROR("Kdtree header was not found in file %s", io->fits->Fptr->filename);

        FREE(kd);
        return NULL;
    }

    kd->has_linear_lr = 0;
    fits_read_key(header->fits, TBYTE, "KDT_LINL", &kd->has_linear_lr, NULL, &status);

    if (p_hdr)
        *p_hdr = header;
    else
        free(header);

    kd->ndata  = ndata;
    kd->ndim   = ndim;
    kd->nnodes = nnodes;
    kd->nbottom = (nnodes+1)/2;
    kd->ninterior = nnodes - kd->nbottom;
    kd->nlevels = kdtree_nnodes_to_nlevels(nnodes);
    kd->treetype = tt;

    KD_DISPATCH(kdtree_read_fits, tt, rtn = , (io, kd));

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
