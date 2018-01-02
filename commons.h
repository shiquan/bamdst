/*
   Copyright (c) 2013-2014 Beijing Genomics Institution (BGI)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Quan Shi (shiquan@genomics.cn) */

#ifndef COMMONS_H
#define COMMONS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <sys/types.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <zlib.h>
#include <getopt.h>
//#include <emmintrin.h>

#ifdef _WIN32
#include <mingw/math.h>
#define drand48() ((double)rand() / RAND_MAX)
#endif

#ifndef NAN
#define NAN (0.0 / 0.0)
#endif

/* Self-defined boolean type */
#define TRUE 1
#define FALSE 0
#define boolean int

#ifndef	__cplusplus
#ifndef bool
# define bool char
#endif
#endif

#if defined(__va_copy) && !defined(va_copy)
#   define va_copy __va_copy
#endif
#if !defined(va_copy)
#   define va_copy(to, from) ((to) = (from))
#endif

#define ArraySize(a) (sizeof(a)/sizeof((a)[0]))

#define isNull(a) ((a) == 0)
#define isZero(x) ((x) == 0)

#define STREQ(a, b) (strcmp((a), (b)) == 0)


void writeout(char *format, ...);
void warnings(char *format, ...);
void errabort(char *format, ...);
void LOG(char *format, ...);
void debug(char const *format, ...);
FILE *stream_open(char const *file, char const *mode);
FILE *open_wfile(char const *file);
FILE *open_rfile(char const *file);
gzFile safe_gzopen(char const *file);
void write_file(FILE *fp, char *format, ...);

void *needmem(size_t size);
void freemem(void *pt);
void mustfree(void *pt);

#define mustfree(pt) do {                                               \
        if (isNull(pt)) {                                               \
            errabort("%s : %d Trying to free an empty point!", __FILE__, __LINE__); \
        }                                                               \
        free(pt);                                                       \
    } while(0)

void *needlargemem(size_t size);
void *resizemem(void *vp, size_t size);
void *resizelargemem(void *vp, size_t size);
void *enlarge_empty_mem(void *vp, size_t old_size, size_t new_size);
void *enlarge_empty_largemem(void *vp, size_t old_size, size_t new_size);

#define errassert(assertion, ...)  do {                                 \
        if (assertion) {                                                \
            fprintf(stderr, __VA_ARGS__);				\
            fprintf(stderr, "\n");                                      \
            assert(0);                                                  \
        } while(0)

#endif
