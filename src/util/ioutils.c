/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <errno.h>

#ifndef _WIN32
#  include <unistd.h>
#endif

#include "os-features.h"
#include "ioutils.h"
#include "errors.h"
#include "log.h"

uint32_t ENDIAN_DETECTOR = 0x01020304;

#include "qsort_reentrant.c"


#ifdef _WIN32
    #include <limits.h>

    int asprintf(char** strp, const char* fmt, ...)
    {
        int r;
        va_list ap;
        va_start(ap, fmt);
        r = vasprintf(strp, fmt, ap);
        va_end(ap);
        return(r);
    }

    int vasprintf(char** strp, const char* fmt, va_list ap)
    {
        int r = -1, size = _vscprintf(fmt, ap);

        if ((size >= 0) && (size < INT_MAX))
        {
            *strp = (char*)malloc(size + 1); //+1 for null
            if (*strp)
            {
                r = vsnprintf(*strp, size + 1, fmt, ap);  //+1 for null
                if ((r < 0) || (r > size))
                {
                    free(*strp);
                    *strp = 0;
                    r = -1;
                }
            }
        }
        else { *strp = 0; }

        return(r);
    }
#endif

void asprintf_safe(char** strp, const char* format, ...) {
    va_list lst;
    int rtn;
    va_start(lst, format);
    rtn = vasprintf(strp, format, lst);
    if (rtn == -1) {
        fprintf(stderr, "Error, vasprintf() failed: %s\n", strerror(errno));
        fprintf(stderr, "  (format: \"%s\")\n", format);
        assert(0);
        *strp = NULL;
    }
    va_end(lst);
}

anbool file_readable(const char* fn) {
#ifndef _WIN32
    return fn && (access(fn, R_OK) == 0);
#else
    return fn && (_access(fn, 04) == 0);
#endif
}

char* strdup_safe(const char* str) {
    char* rtn;
    if (!str) return NULL;
    rtn = strdup(str);
    if (!rtn) {
        fprintf(stderr, "Failed to strdup: %s\n", strerror(errno));
        assert(0);
    }
    return rtn;
}
