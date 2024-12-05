/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef QUADFILE_H
#define QUADFILE_H

#include <sys/types.h>
#include <stdint.h>
#include "fits_io.h"

typedef struct {
    unsigned int numquads;
    unsigned int numstars;
    int dimquads;
    // upper bound of AB distance of quads in this index
    double index_scale_upper;
    // lower bound
    double index_scale_lower;
    // unique ID of this index
    int indexid;
    // healpix covered by this index
    int healpix;
    // Nside of the healpixelization
    int hpnside;

    fits_file_t* io;
    // when reading:
    uint32_t* quadarray;
} quadfile_t;

quadfile_t* quadfile_open_fits(const char* filename, fitsfile* fits);

char* quadfile_get_filename(const quadfile_t* qf);

int quadfile_close(quadfile_t* qf);

// Look at each quad, and ensure that the star ids it contains are all
// less than the number of stars ("numstars").  Returns 0=ok, -1=problem
int quadfile_check(const quadfile_t* qf);

// Copies the star ids of the stars that comprise quad "quadid".
// There will be qf->dimquads such stars.
// (this will be less than starutil.h : DQMAX, for ease of static
// allocation of arrays that will hold quads of stars)
int quadfile_get_stars(const quadfile_t* qf, unsigned int quadid,
                       unsigned int* stars);

int quadfile_dimquads(const quadfile_t* qf);

int quadfile_nquads(const quadfile_t* qf);

double quadfile_get_index_scale_upper_arcsec(const quadfile_t* qf);

double quadfile_get_index_scale_lower_arcsec(const quadfile_t* qf);

fits_hdu_t* quadfile_get_header(const quadfile_t* qf);

#endif
