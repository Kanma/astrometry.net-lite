/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "quadfile.h"
#include "starutil.h"
#include "ioutils.h"
#include "errors.h"
#include "an-endian.h"

#define CHUNK_QUADS 0


void quadfile_parse(fitsfile* fits, fits_file_t* io) {
    io->quads.dimquads = 4;
    io->quads.numquads = -1;
    io->quads.numstars = -1;
    io->quads.index_scale_upper = -1.0;
    io->quads.index_scale_lower = -1.0;
    io->quads.indexid = 0;
    io->quads.healpix = -1;
    io->quads.hpnside = 1;

    int status = 0;
    fits_movabs_hdu(fits, 1, NULL, &status);

    if (status == 0)
    {
        fits_read_key(fits, TINT, "DIMQUADS", &io->quads.dimquads, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "NQUADS", &io->quads.numquads, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "NSTARS", &io->quads.numstars, NULL, &status);

        status = 0;
        fits_read_key(fits, TDOUBLE, "SCALE_U", &io->quads.index_scale_upper, NULL, &status);

        status = 0;
        fits_read_key(fits, TDOUBLE, "SCALE_L", &io->quads.index_scale_lower, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "INDEXID", &io->quads.indexid, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "HEALPIX", &io->quads.healpix, NULL, &status);

        status = 0;
        fits_read_key(fits, TINT, "HPNSIDE", &io->quads.hpnside, NULL, &status);
    }
}


static quadfile_t* new_quadfile(fits_file_t* io) {
    quadfile_t* qf;
    qf = calloc(1, sizeof(quadfile_t));
    if (!qf) {
        SYSERROR("Couldn't malloc a quadfile struct");
        return NULL;
    }

    qf->io = io;

    qf->dimquads = io->quads.dimquads;
    qf->numquads = io->quads.numquads;
    qf->numstars = io->quads.numstars;
    qf->index_scale_upper = io->quads.index_scale_upper;
    qf->index_scale_lower = io->quads.index_scale_lower;
    qf->indexid = io->quads.indexid;
    qf->healpix = io->quads.healpix;
    qf->hpnside = io->quads.hpnside;

    if ((qf->numquads == -1) || (qf->numstars == -1) ||
        (qf->index_scale_upper == -1.0) || (qf->index_scale_lower == -1.0)) {
        ERROR("Couldn't find NQUADS or NSTARS or SCALE_U or SCALE_L entries in FITS header");
        {
            quadfile_close(qf);
            return NULL;
        }
    }

    if (fits_check_endian(&io->hdus[0])) {
        ERROR("Quad file was written with the wrong endianness");
        {
            quadfile_close(qf);
            return NULL;
        }
    }

    return qf;
}

int quadfile_check(const quadfile_t* qf) {
    int q, i;
    if (qf->dimquads < 3 || qf->dimquads > DQMAX) {
        ERROR("Dimquads has illegal value %i", qf->dimquads);
        return -1;
    }
    for (q=0; q<qf->numquads; q++) {
        unsigned int stars[DQMAX];
        if (quadfile_get_stars(qf, q, stars)) {
            ERROR("Failed to get quad %i of %i", q, qf->numquads);
            return -1;
        }
        for (i=0; i<qf->dimquads; i++) {
            if (stars[i] >= qf->numstars) {
                ERROR("Star ID %i is out of bounds: num stars %i", stars[i], qf->numstars);
                return -1;
            }
        }
    }
    return 0;
}

int quadfile_dimquads(const quadfile_t* qf) {
    return qf->dimquads;
}

int quadfile_nquads(const quadfile_t* qf) {
    return qf->numquads;
}

fits_hdu_t* quadfile_get_header(const quadfile_t* qf) {
    return fits_get_primary_header(qf->io);
}

char* quadfile_get_filename(const quadfile_t* qf) {
    return qf->io->filename;
}

quadfile_t* quadfile_open_fits(fits_file_t* io) {
    quadfile_t* qf = NULL;
 
    qf = new_quadfile(io);
    if (!qf)
        goto bailout;

    fits_read_chunk(io, "quads", qf->dimquads * sizeof(uint32_t), &qf->numquads, &qf->quadarray, 1);

    return qf;

 bailout:
    if (qf)
        quadfile_close(qf);
    return NULL;
}

int quadfile_close(quadfile_t* qf) {
    int rtn = 0;
    if (!qf) return 0;
    free(qf);
    return rtn;
}

double quadfile_get_index_scale_upper_arcsec(const quadfile_t* qf) {
    return rad2arcsec(qf->index_scale_upper);
}

double quadfile_get_index_scale_lower_arcsec(const quadfile_t* qf) {
    return rad2arcsec(qf->index_scale_lower);
}

int quadfile_get_stars(const quadfile_t* qf, unsigned int quadid, unsigned int* stars) {
    int i;
    if (quadid >= qf->numquads) {
        ERROR("Requested quadid %i, but number of quads is %i",	quadid, qf->numquads);
        assert(quadid < qf->numquads);
        return -1;
    }

    for (i=0; i<qf->dimquads; i++) {
        stars[i] = qf->quadarray[quadid * qf->dimquads + i];
    }
    return 0;
}



