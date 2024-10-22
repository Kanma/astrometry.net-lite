/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include "quad-utils.h"

void quad_flip_parity(const double* code, double* flipcode, int dimcode) {
    int i;
    // swap CX <-> CY, DX <-> DY.
    for (i=0; i<dimcode/2; i++) {
        // use tmp in code "code" == "flipcode"
        double tmp;
        tmp = code[2*i+1];
        flipcode[2*i+1] = code[2*i+0];
        flipcode[2*i+0] = tmp;
    }
}
