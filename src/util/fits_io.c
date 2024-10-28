/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#include "fits_io.h"
#include "kdtree_mem.h"
#include "ioutils.h"
#include "errors.h"
#include "tic.h"
#include "log.h"


fits_io_t* fits_open_fits(const char* filename, fitsfile* fits) {
    int status = 0;

    fits_io_t* io = malloc(sizeof(fits_io_t));
    io->fits = fits;

    fits_get_num_hdus(io->fits, &io->nbHDUs, &status);

    io->filename = malloc(strlen(filename) + 1);
    strcpy(io->filename, filename);

    return io;
}

fits_hdu_t* fits_get_primary_header(fits_io_t* io) {
    int status = 0;

    fits_hdu_t* header = malloc(sizeof(fits_hdu_t));
    header->fits = io->fits;
    header->extension = 1;

    return header;
}

int fits_read_chunk(fits_io_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data) {
    int status = 0;
    int hdutype = 0;
    char str[FITS_LINESZ + 1];

    int found = 0;
    for (int i = 2; i <= io->nbHDUs; ++i)
    {
        status = 0;
        fits_movabs_hdu(io->fits, i, &hdutype, &status);

        if ((status == 0) && (hdutype == BINARY_TBL))
        {
            fits_read_key(io->fits, TSTRING, "TTYPE1", str, NULL, &status);
            if ((status == 0) && (strcmp(tablename, str) == 0))
            {
                found = 1;
                break;
            }
        }
    }

    if (found == 0)
        return -1;

    if (*nbRows == 0)
    {
        fits_get_num_rows(io->fits, nbRows, &status);
        if (status != 0)
            return -1;
    }

    *data = malloc(itemSize * *nbRows);

    fits_read_tblbytes(io->fits, 1, 1, itemSize * *nbRows, *data, &status);
    if (status != 0)
    {
        free(*data);
        return -1;
    }

    return 0;
}

int fits_check_endian(const fits_hdu_t* header) {
    char filestr[FITS_LINESZ+1];
    char localstr[FITS_LINESZ+1];

    int status = 0;
    fits_movabs_hdu(header->fits, header->extension, NULL, &status);

    fits_read_key(header->fits, TSTRING, "ENDIAN", filestr, NULL, &status);
    if (status != 0) {
        // No ENDIAN header found.
        return 1;
    }

    uint32_t endian = ENDIAN_DETECTOR;
    unsigned char* cptr = (unsigned char*) &endian;
    sprintf(localstr, "%02x:%02x:%02x:%02x", (uint32_t)cptr[0], (uint32_t)cptr[1], (uint32_t)cptr[2], (uint32_t)cptr[3]);

    if (strcmp(filestr, localstr)) {
        fprintf(stderr, "File was written with endianness %s, this machine has endianness %s.\n", filestr, localstr);
        return -1;
    }

    return 0;
}

int fits_io_close(fits_io_t* io) {
    free(io->filename);
    free(io);
    return 0;
}
