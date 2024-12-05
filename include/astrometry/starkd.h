/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef STAR_KD_H
#define STAR_KD_H

#include <stdint.h>
#include <stdio.h>

#include "astrometry/kdtree.h"
#include "astrometry/kdtree_fits_io.h"
#include "astrometry/keywords.h"

#ifdef SWIG
// this keyword (from keywords.h) confuses swig
#define Malloc
#endif

#define AN_FILETYPE_STARTREE "SKDT"

#define AN_FILETYPE_TAGALONG "TAGALONG"

#define STARTREE_NAME "stars"

typedef struct {
    kdtree_t* tree;
    fits_hdu_t* header;
    int* inverse_perm;
    uint8_t* sweep;
} startree_t;

startree_t* startree_open_fits(const char* filename, fits_file_t* fits);

/**
 Searches for stars within a radius of a point.

 xyzcenter: double[3]: unit-sphere coordinates of point; see
 starutil.h : radecdeg2xyzarr() to convert RA,Decs to this form.

 radius2: radius-squared on the unit sphere; see starutil.h :
 deg2distsq() or arcsec2distsq().

 xyzresults: if non-NULL, returns the xyz positions of the stars that
 are found, in a newly-allocated array.

 radecresults: if non-NULL, returns the RA,Dec positions (in degrees)
 of the stars within range.

 starinds: if non-NULL, returns the indices of stars within range.
 This can be used to retrieve extra information about the stars, using
 the 'startree_get_data_column()' function.
 
 */
void startree_search_for(const startree_t* s, const double* xyzcenter, double radius2,
                         double** xyzresults, double** radecresults,
                         int** starinds, int* nresults);

/**
 RA, Dec, and radius in degrees.  Otherwise same as startree_search_for().
 */
void startree_search_for_radec(const startree_t* s, double ra, double dec, double radius,
                               double** xyzresults, double** radecresults,
                               int** starinds, int* nresults);

void startree_search(const startree_t* s, const double* xyzcenter, double radius2,
                     double** xyzresults, double** radecresults, int* nresults);

/**
 Reads a column of data from the "tag-along" table.

 Get the "inds" and "N" from "startree_search" or "startree_search_for".

 To get all entries, set "inds" = NULL and N = startree_N().

 The return value is a newly-allocated array of size N.  It should be
 freed using "startree_free_data_column"
 */
Malloc
double* startree_get_data_column(startree_t* s, const char* colname, const int* indices, int N);

/**
 Same as startree_get_data_column but for int64_t.  Don't you love C templating?
 */
Malloc
int64_t* startree_get_data_column_int64(startree_t* s, const char* colname, const int* indices, int N);

/**
 Reads a column of data from the "tag-along" table.

 The column may be an array (that is, each row contains multiple
 entries); the array size is placed in "arraysize".

 The array entries 
 */
Malloc
double* startree_get_data_column_array(startree_t* s, const char* colname, const int* indices, int N, int* arraysize);

void startree_free_data_column(startree_t* s, double* d);


/*
 Retrieve parameters of the cut-an process, if they are available.
 Older index files may not have these header cards.
 */
// healpix nside, or -1
int startree_get_cut_nside(const startree_t* s);
int startree_get_cut_nsweeps(const startree_t* s);
// in arcsec; 0 if none.
double startree_get_cut_dedup(const startree_t* s);
// band (one of several static strings), or NULL
char* startree_get_cut_band(const startree_t* s);
// margin, in healpix, or -1
int startree_get_cut_margin(const startree_t* s);

double startree_get_jitter(const startree_t* s);

//uint64_t startree_get_starid(const startree_t* s, int ind);

// returns the sweep number of star 'ind', or -1 if the index is out of bounds
// or the tree has no sweep numbers.
int startree_get_sweep(const startree_t* s, int ind);

int startree_N(const startree_t* s);

int startree_nodes(const startree_t* s);

int startree_D(const startree_t* s);

fits_hdu_t* startree_header(const startree_t* s);

int startree_get(startree_t* s, int starid, double *p_xyz);

int startree_get_radec(startree_t* s, int starid, double *p_ra, double *p_dec);

int startree_close(startree_t* s);

void startree_compute_inverse_perm(startree_t* s);

int startree_check_inverse_perm(startree_t* s);

#endif
