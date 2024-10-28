/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef IOUTILS_H
#define IOUTILS_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <time.h>

#include "astrometry/an-bool.h"
#include "astrometry/bl.h"
#include "astrometry/keywords.h"

extern uint32_t ENDIAN_DETECTOR;

void QSORT_R(void* base, size_t nmembers, size_t member_size,
             void* token, int (*compar)(void *, const void *, const void *));

/**
   You should define the "comparison" function like this:

   static int QSORT_COMPARISON_FUNCTION(my_comparison, void* token, const void* v1, const void* v2) {
 */
#define QSORT_COMPARISON_FUNCTION(func, thunk, v1, v2) func(thunk, v1, v2)

/*
 It's not really _safe_ as such, it just prints an error message if it fails...
 */
void
ATTRIB_FORMAT(printf,2,3)
    asprintf_safe(char** strp, const char* format, ...);

anbool file_readable(const char* fn);

char* strdup_safe(const char* str);

#endif
