/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef FITS_IO_H
#define FITS_IO_H

#include <fitsio.h>

#define FITS_LINESZ 80


struct fits_hdu_t
{
    fitsfile* fits;
    int extension;
};

typedef struct fits_hdu_t fits_hdu_t;


struct fits_io_t
{
    char* filename;
    fitsfile* fits;
    int nbHDUs;
};

typedef struct fits_io_t fits_io_t;


fits_io_t* fits_open_fits(const char* filename, fitsfile* fits);

int fits_read_chunk(fits_io_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data);

fits_hdu_t* fits_get_primary_header(fits_io_t* io);

int fits_check_endian(const fits_hdu_t* header);

int fits_io_close(fits_io_t* io);

#endif
