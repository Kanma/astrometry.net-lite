/*
 # This file is part of libkd.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

#ifndef _WIN32
    #include <sys/mman.h>
    #include <fcntl.h>
    #include <sys/errno.h>
    #include <unistd.h>
#else _WIN32
    #include <windows.h>
#endif

#include "fits_io.h"
#include "kdtree_mem.h"
#include "ioutils.h"
#include "errors.h"
#include "tic.h"
#include "log.h"


fits_file_t* fits_open(const char* filename) {
    int status = 0;

    fits_file_t* io = malloc(sizeof(fits_file_t));

    fits_open_file(&io->fits, filename, READONLY, &status);
    if (status != 0) {
        ERROR("Failed to open FITS file %s", filename);
        goto bailout;
    }

    fits_get_num_hdus(io->fits, &io->nbHDUs, &status);

#ifndef _WIN32
    io->fd = -1;
#else
    io->fd = INVALID_HANDLE_VALUE;
    io->mapping = INVALID_HANDLE_VALUE;
#endif

    return io;

bailout:
    free(io);
    return NULL;
}

fits_hdu_t* fits_get_primary_header(fits_file_t* io) {
    int status = 0;

    fits_hdu_t* header = malloc(sizeof(fits_hdu_t));
    header->fits = io;
    header->extension = 1;

    return header;
}

int fits_read_chunk(fits_file_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data) {
    int status = 0;
    int hdutype = 0;
    char str[FITS_LINESZ + 1];

    int found = 0;
    for (int i = 2; i <= io->nbHDUs; ++i)
    {
        status = 0;
        fits_movabs_hdu(io->fits, i, &hdutype, &status);

        if ((status == 0) && (hdutype == BINARY_TBL))
        {
            fits_read_key(io->fits, TSTRING, "TTYPE1", str, NULL, &status);
            if ((status == 0) && (strcmp(tablename, str) == 0))
            {
                found = 1;
                break;
            }
        }
    }

    if (found == 0)
        return -1;

    if (*nbRows == 0)
    {
        fits_get_num_rows(io->fits, nbRows, &status);
        if (status != 0)
            return -1;
    }

    OFF_T headerStart, dataStart, dataEnd;
    fits_get_hduoff(io->fits, &headerStart, &dataStart, &dataEnd, &status);
    if (status != 0)
        return -1;

#ifndef _WIN32
    if (io->fd == -1)
        io->fd = open(io->fits->Fptr->filename, O_RDONLY, 0);

    long int pageSize = sysconf(_SC_PAGESIZE);
    long int offset = dataStart % pageSize;

    *data = mmap(NULL, dataEnd - dataStart + offset, PROT_READ, MAP_PRIVATE, io->fd, dataStart - offset) + offset;
#else
    if (io->fd == INVALID_HANDLE_VALUE)
    {
        io->fd = CreateFileA(
            io->fits->Fptr->filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL
        );
    }

    if (io->mapping == INVALID_HANDLE_VALUE)
        io->mapping = CreateFileMappingA(io->fd, NULL, PAGE_READONLY, 0, 0, NULL);

    SYSTEM_INFO infos;
    GetSystemInfo(&infos);

    long int offset = dataStart % infos.dwAllocationGranularity;

    *data = (char*) MapViewOfFile(io->mapping, FILE_MAP_READ, 0, dataStart - offset, dataEnd - dataStart + offset) + offset;
#endif

    return 0;
}

int fits_check_endian(const fits_hdu_t* header) {
    char filestr[FITS_LINESZ+1];
    char localstr[FITS_LINESZ+1];

    int status = 0;
    fits_movabs_hdu(header->fits, header->extension, NULL, &status);

    fits_read_key(header->fits, TSTRING, "ENDIAN", filestr, NULL, &status);
    if (status != 0) {
        // No ENDIAN header found.
        return 1;
    }

    uint32_t endian = ENDIAN_DETECTOR;
    unsigned char* cptr = (unsigned char*) &endian;
    sprintf(localstr, "%02x:%02x:%02x:%02x", (uint32_t)cptr[0], (uint32_t)cptr[1], (uint32_t)cptr[2], (uint32_t)cptr[3]);

    if (strcmp(filestr, localstr)) {
        fprintf(stderr, "File was written with endianness %s, this machine has endianness %s.\n", filestr, localstr);
        return -1;
    }

    return 0;
}

int fits_close(fits_file_t* io) {
#ifndef _WIN32
    if (io->fd >= 0)
        close(io->fd);
#else
    if (io->mapping != INVALID_HANDLE_VALUE)
        CloseHandle(io->mapping);

    if (io->fd != INVALID_HANDLE_VALUE)
        CloseHandle(io->fd);
#endif

    if (io->fits)
    {
        int status = 0;
        fits_close_file(io->fits, &status);
    }

    free(io);

    return 0;
}
