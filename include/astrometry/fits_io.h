/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef FITS_IO_H
#define FITS_IO_H

#include <fitsio.h>

#define FITS_LINESZ 80


struct fits_file_t
{
    fitsfile* fits;
    int nbHDUs;
#ifndef _WIN32
    int fd;
#else
    void* fd;
    void* mapping;
#endif
};

typedef struct fits_file_t fits_file_t;


struct fits_hdu_t
{
    fits_file_t* fits;
    int extension;
};

typedef struct fits_hdu_t fits_hdu_t;


fits_file_t* fits_open(const char* filename);

int fits_read_chunk(fits_file_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data);

fits_hdu_t* fits_get_primary_header(fits_file_t* io);

int fits_check_endian(const fits_hdu_t* header);

int fits_close(fits_file_t* io);

#endif
