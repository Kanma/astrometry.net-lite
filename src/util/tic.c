/*
 # This file is part of the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifndef _WIN32
#  include <sys/time.h>
#  include <sys/resource.h>
#  include <unistd.h>
#else
#  include <winsock.h>
#endif

#include "tic.h"
#include "errors.h"
#include "log.h"

static time_t starttime;
static double starttime2;
static double startutime, startstime;

double timenow() {
#ifndef _WIN32
    struct timeval tv;
    if (gettimeofday(&tv, NULL)) {
        ERROR("Failed to get time of day");
        return -1.0;
    }
    return (double)(tv.tv_sec - 3600*24*365*30) + tv.tv_usec * 1e-6;
#else
    FILETIME ft_now;
    GetSystemTimeAsFileTime(&ft_now);
    LONGLONG ll_now = (LONGLONG)ft_now.dwLowDateTime + ((LONGLONG)(ft_now.dwHighDateTime) << 32LL);
    return (double)(ll_now) / 10000000.0;
#endif
}

int get_resource_stats(double* p_usertime, double* p_systime, long* p_maxrss) {
#ifndef _WIN32
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage)) {
        SYSERROR("Failed to get resource stats (getrusage)");
        return 1;
    }
    if (p_usertime) {
        *p_usertime = usage.ru_utime.tv_sec + 1e-6 * usage.ru_utime.tv_usec;
    }
    if (p_systime) {
        *p_systime = usage.ru_stime.tv_sec + 1e-6 * usage.ru_stime.tv_usec;
    }
    if (p_maxrss) {
        *p_maxrss = usage.ru_maxrss;
    }
    return 0;
#else
    clock_t elapsed = clock();
    *p_usertime = (double)(elapsed) / CLOCKS_PER_SEC;
    *p_systime = 0.0;
    *p_maxrss = 0.0;
    return 0;
#endif
}
