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
#include "kdtree_fits_io.h"
#include "kdtree_mem.h"
#include "quadfile.h"
#include "starkd.h"
#include "codekd.h"
#include "ioutils.h"
#include "errors.h"
#include "tic.h"
#include "log.h"


static void parse_header(fitsfile* fits, fits_hdu_t* header) {
    char str[FITS_LINESZ+1];
    char localstr[FITS_LINESZ+1];

    int status = 0;
    fits_movabs_hdu(fits, header->extension, &header->hdutype, &status);

    if ((status == 0) && (header->hdutype == BINARY_TBL))
    {
        fits_read_key(fits, TSTRING, "TTYPE1", header->table.name, NULL, &status);
        fits_get_num_rows(fits, &header->table.nbRows, &status);
    }

    status = 0;

    fits_get_hduoff(fits, &header->headerStart, &header->dataStart, &header->dataEnd, &status);

    fits_read_key(fits, TSTRING, "ENDIAN", str, NULL, &status);
    if (status != 0) {
        // No ENDIAN header found.
        header->endian = 1;
    } else {
        uint32_t endian = ENDIAN_DETECTOR;
        unsigned char* cptr = (unsigned char*) &endian;
        sprintf(localstr, "%02x:%02x:%02x:%02x", (uint32_t)cptr[0], (uint32_t)cptr[1], (uint32_t)cptr[2], (uint32_t)cptr[3]);

        if (strcmp(str, localstr)) {
            fprintf(stderr, "File was written with endianness %s, this machine has endianness %s.\n", str, localstr);
            header->endian = -1;
        } else {
            header->endian = 0;
        }
    }
}


fits_file_t* fits_open(const char* filename) {
    int status = 0;

    fits_file_t* io = malloc(sizeof(fits_file_t));
    memset(io, 0, sizeof(fits_file_t));

    fitsfile* fits = NULL;

    fits_open_file(&fits, filename, READONLY, &status);
    if (status != 0) {
        ERROR("Failed to open FITS file %s", filename);
        goto bailout;
    }

    io->filename = strdup(filename);

    fits_get_num_hdus(fits, &io->nbHDUs, &status);
    if (status != 0) {
        ERROR("Failed to retrieve the number of HDUs in the FITS file %s", filename);
        goto bailout;
    }

    io->hdus = malloc(io->nbHDUs * sizeof(fits_hdu_t));
    memset(io->hdus, 0, io->nbHDUs * sizeof(fits_hdu_t));

    for (int i = 0; i < io->nbHDUs; ++i) {
        io->hdus[i].fits = io;
        io->hdus[i].extension = i + 1;
        io->hdus[i].extension = i + 1;
        io->hdus[i].tree.ndim = -1;
        io->hdus[i].tree.ndata = -1;
        io->hdus[i].tree.nnodes = -1;

        parse_header(fits, &io->hdus[i]);
    }

    io->stars.cut_nside = -1;
    io->stars.cut_nsweeps = -1;
    io->stars.cut_dedup = 0.0;
    io->stars.cut_band = NULL;
    io->stars.cut_margin = -1;
    io->stars.jitter = 0.0;

#ifndef _WIN32
    io->fd = -1;
#else
    io->fd = INVALID_HANDLE_VALUE;
    io->mapping = INVALID_HANDLE_VALUE;
#endif

    quadfile_parse(fits, io);
    kdtree_parse(fits, io);

    fits_hdu_t* header = NULL;
    if (kdtree_fits_contains_tree(io, STARTREE_NAME, &header))
        parse_startree_params(fits, header);
    else
        parse_startree_params(fits, &io->hdus[0]);

    if (kdtree_fits_contains_tree(io, CODETREE_NAME, &header))
        parse_codetree_params(fits, header);
    else
        parse_codetree_params(fits, &io->hdus[0]);

    fits_close_file(fits, &status);

    return io;

bailout:
    if (io->hdus)
        free(io->hdus);

    if (io->filename)
        free(io->filename);

    free(io);

    if (fits)
        fits_close_file(fits, &status);

    return NULL;
}

fits_hdu_t* fits_get_primary_header(fits_file_t* io) {
    return &io->hdus[0];
}

int fits_read_chunk(fits_file_t* io, const char* tablename, size_t itemSize, int* nbRows, void** data, int closeFile) {
    int status = 0;
    int hdutype = 0;
    char str[FITS_LINESZ + 1];

    fits_hdu_t* header = NULL;

    int found = 0;
    for (int i = 1; i <= io->nbHDUs; ++i)
    {
        header = &io->hdus[i];

        if (header->hdutype == BINARY_TBL)
        {
            if ((status == 0) && (strcmp(tablename, header->table.name) == 0))
            {
                found = 1;
                break;
            }
        }
    }

    if (found == 0)
        return -1;

    if (*nbRows == 0)
        *nbRows = header->table.nbRows;

#ifndef _WIN32
    if (io->fd == -1) {
        io->fd = open(io->filename, O_RDONLY, 0);

        if (io->fd == -1) {
            ERROR("Failed to open the file '%s', error=%d", io->filename, errno);
            return -1;
        }
    }

    long int pageSize = sysconf(_SC_PAGESIZE);
    long int offset = header->dataStart % pageSize;

    *data = mmap(NULL, header->dataEnd - header->dataStart + offset, PROT_READ, MAP_PRIVATE, io->fd, header->dataStart - offset);
    if (*data != MAP_FAILED) {
        *data = (char*) *data + offset;
    } else {
        *data = NULL;
        ERROR("Failed to mmap, error=%d", errno);
    }

    if (closeFile) {
        close(io->fd);
        io->fd = -1;
    }
#else
    if (io->fd == INVALID_HANDLE_VALUE)
    {
        io->fd = CreateFileA(
            io->filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL
        );

        if (io->fd == INVALID_HANDLE_VALUE) {
            ERROR("Failed to open the file '%s', error=%d", io->filename, GetLastError());
            return -1;
        }
    }

    if (io->mapping == INVALID_HANDLE_VALUE) {
        io->mapping = CreateFileMappingA(io->fd, NULL, PAGE_READONLY, 0, 0, NULL);

        if (io->mapping == NULL) {
            ERROR("Failed to CreateFileMappingA, error=%d", GetLastError());
            CloseHandle(io->fd);
            io->mapping = INVALID_HANDLE_VALUE;
            io->fd = INVALID_HANDLE_VALUE;
            return -1;
        }
    }

    SYSTEM_INFO infos;
    GetSystemInfo(&infos);

    long int offset = header->dataStart % infos.dwAllocationGranularity;

    *data = (char*) MapViewOfFile(io->mapping, FILE_MAP_READ, 0, header->dataStart - offset, header->dataEnd - header->dataStart + offset);
    if (*data != NULL) {
         *data = (char*) *data + offset;
    } else {
        ERROR("Failed to MapViewOfFile, error=%d", GetLastError());
    }

     + offset;

    if (closeFile) {
        CloseHandle(io->mapping);
        CloseHandle(io->fd);
        io->mapping = INVALID_HANDLE_VALUE;
        io->fd = INVALID_HANDLE_VALUE;
    }
#endif

    if (*data != NULL)
        return 0;

    return -1;
}

int fits_check_endian(const fits_hdu_t* header) {
    return header->endian;
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

    if (io->hdus)
        free(io->hdus);

    if (io->filename)
        free(io->filename);

    free(io);

    return 0;
}
