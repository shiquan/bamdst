/* The MIT License

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

/* The library is used for count numbers. 
 * Support 32bit for normal number, and 64bit for bigger number
 */
#ifndef COUNT_H
#define COUNT_H
#include "commons.h"

#define COUNTINT 8
#define COUNTINT64 16

typedef struct
  {
  uint32_t m, n;
  uint32_t *a;
  }
count32_t;

typedef struct
  {
  uint32_t m, n;
  uint64_t *a;
  }
count64_t;

#define count32_init(c)	do {			\
 (c) = (count32_t *)malloc(sizeof(count32_t));	\
 (c)->m = 0;					\
 (c)->n = 2;					\
 (c)->a = malloc(sizeof(uint32_t)*2);		\
 memset((c)->a, 0, (c)->n *sizeof(uint32_t));	\
 } while(0)

#define count64_init(c)	do {			\
 (c) = (count64_t*)malloc(sizeof(count64_t));	\
 (c)->m = 0;					\
 (c)->n = 2;					\
 (c)->a = malloc(sizeof(uint64_t)*2);		\
 memset((c)->a, 0, (c)->n *sizeof(uint64_t));	\
 } while(0)

#define count_destroy(c) { free((c)->a); free(c); }

#define count_zero(c) do {			\
 int i;						\
 for (i = 0; i < (c)->n; ++i) (c)->a[i] = 0;	\
 } while(0)


// increase c once in postion d
#define count_increase(c, d, type) do {				\
 if ((c)->n == (c)->m)						\
   {								\
   if ((c)->m)							\
     {								\
     (c)->n = (c)->m << 1;					\
     (c)->a = realloc((c)->a, sizeof(type) * (c)->n);		\
     memset((c)->a + (c)->m, 0, ((c)->n-(c)->m)*sizeof(type));	\
     }								\
   }								\
 if((c)->n <= d+1)						\
   {								\
   (c)->n = (d+1)<<1;						\
   (c)->a = realloc((c)->a, sizeof(type) *(c)->n);		\
   memset((c)->a+(c)->m, 0, ((c)->n - (c)->m)*sizeof(type));	\
   (c)->m = d+1;						\
   }								\
 else if((c)->m < d+1) (c)->m = d+1;				\
 (c)->a[d]++;							\
 } while(0)

#define count_resize(c, l, type) do {		\
 if((c)->n < l)					\
   {						\
   (c)->a = realloc((c)->a, sizeof(type) * l);	\
   int i;					\
   for (i = (c)->n; i < l; ++i) (c)->a[i] = 0;	\
   (c)->n = l;					\
   }						\
 } while(0);

#define count_increaseN(c, d, cnt, type) do {			\
 if((c)->n <= d+1)						\
   {								\
   (c)->n = (d+1)<<1;						\
   (c)->a = realloc((c)->a, sizeof(type) *(c)->n);		\
   memset((c)->a+(c)->m, 0, ((c)->n - (c)->m)* sizeof(type));	\
   (c)->m = d+1;						\
   }								\
 else if((c)->m < d+1) (c)->m = d+1;				\
 (c)->a[d] += cnt;						\
 } while(0)

#endif
