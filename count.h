/* count.h - macros for */
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

#define count32_init(c)	do {								\
 (c) = (count32_t *)malloc(sizeof(count32_t));				\
 (c)->m = (c)->n = 0;										\
 } while(0)

#define count64_init(c)	do {						\
 (c) = (count64_t*)malloc(sizeof(count64_t));		\
 (c)->m = (c)->n = 0;								\
 } while(0)

#define count_destroy(c) { free((c)->a); free(c); }

#define count_zero(c) do {						\
 int i;											\
 for (i = 0; i < (c)->n; ++i) (c)->a[i] = 0;	\
 } while(0)

#define count_merge(ca, cb, type) do {					\
 int i;													\
 if ((ca)->m < (cb)->m)									\
	 {													\
	 (ca)->n = (cb)->m + 1;								\
	 (ca)->a = realloc((a)->a, (a)->n * sizeof(type));	\
	 for (i = (ca)->m; i < (cb)->n; ++i) (ca)->a[i] = 0;	\
	 (ca)->m = (cb)->m;									\
	 }													\
 for (i = 0; i < (cb)->m; ++i)	(ca)->[i] += (cb)->a[i];	\
 } while(0)

#define count_increase(c, d, type) do {								\
 if ((c)->n == (c)->m)												\
	 {																\
	 if ((c)->m)													\
		 {															\
		 (c)->n = (c)->m << 1;										\
		 (c)->a = realloc((c)->a, sizeof(type) * (c)->n);			\
		 memset((c)->a + (c)->m, 0, ((c)->n-(c)->m)*sizeof(type));	\
		 }													\
	 else													\
		 {													\
		 (c)->n = 2;										\
		 (c)->a = malloc(sizeof(type)*2);					\
		 memset((c)->a, 0, (c)->n *sizeof(type));			\
		 }														\
	 }															\
 if((c)->n <= d)												\
	 {															\
	 (c)->n = d<<1;												\
	 (c)->a = realloc((c)->a, sizeof(type) *(c)->n);			\
	 memset((c)->a+(c)->m, 0, ((c)->n - (c)->m)*sizeof(type));	\
	 (c)->m = d+1;											\
	 }															\
 else if((c)->m < d) (c)->m = d+1;								\
 (c)->a[d]++;													\
 } while(0)

#define count_resize(c, l, type) do {				\
 if((c)->n < l)										\
	 {												\
	 (c)->a = realloc((c)->a, sizeof(type) * l);	\
	 int i;											\
	 for (i = (c)->n; i < l; ++i) (c)->a[i] = 0;	\
	 (c)->n = l;									\
	 }												\
 } while(0);

#define count_increaseN(c, d, cnt, type) do {			\
 if((c)->n <= d)										\
	 {													\
	 (c)->n = d<<1;										\
	 (c)->a = realloc((c)->a, sizeof(type) *(c)->n);	\
	 int i;												\
	 for (i = (c)->m; i < (c)->n; ++i)	(c)->a[i] = 0;	\
	 (c)->m = d;										\
	 }													\
 else if((c)->m < d) (c)->m = d;						\
 (c)->a[d] += cnt;										\
 } while(0)

#endif
