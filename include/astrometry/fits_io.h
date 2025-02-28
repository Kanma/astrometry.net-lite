/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef FITS_IO_H
#define FITS_IO_H

#include <fitsio.h>
#include "an-bool.h"

#define FITS_LINESZ 80


struct fits_hdu_t;
typedef struct fits_hdu_t fits_hdu_t;


struct fits_file_t
{
    char* filename;

    int nbHDUs;
    fits_hdu_t* hdus;

    struct {
        unsigned int numquads;
        unsigned int numstars;
        int dimquads;
        double index_scale_upper;
        double index_scale_lower;
        int indexid;
        int healpix;
        int hpnside;
    } quads;

    struct {
        int cut_nside;
        int cut_nsweeps;
        double cut_dedup;
        char* cut_band;
        int cut_margin;
        double jitter;
    } stars;

    struct {
        anbool circle;
        anbool cx_less_than_dx;
        anbool meanx_less_than_half;
    } code;

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
    int hdutype;
    int endian;

    OFF_T headerStart;
    OFF_T dataStart;
    OFF_T dataEnd;

    struct table_t {
        char name[FITS_LINESZ + 1];
        int nbRows;
    } table;

    struct tree_t {
        char name[FITS_LINESZ+1];
        int ndim;
        int ndata;
        int nnodes;
        unsigned int treetype;
        int has_linear_lr;
    } tree;
};

//typedef struct fits_hdu_t fits_hdu_t;


fits_file_t* fits_open(const char* filename);

int fits_read_chunk(fits_file_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data, int closeFile);

fits_hdu_t* fits_get_primary_header(fits_file_t* io);

int fits_check_endian(const fits_hdu_t* header);

int fits_close(fits_file_t* io);

#endif
