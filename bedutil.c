//
//  Created by shiquan on 5/18/14.
//  Copyright (c) 2014 ___Beijing Genomics Institution (BGI)___.
//

#include "commons.h"
#include "bedutil.h"
#include "ksort.h"
#include "kseq.h"

KSORT_INIT_GENERIC(uint64_t)

KSTREAM_INIT(gzFile, gzread, 8192)

#define LIDX_SHIFT 13
#define swapvalue(a, b, t) {t c = a; a = b; b = c;}

static inf_t *inf_init()
  {
  inf_t *inf;
  inf = (inf_t*)needmem(sizeof(inf_t));
  inf->total = inf->length = inf->region = 0;
  return inf;
  }

static void sort_reg(bedreglist_t *reg)
  {
  ks_introsort(uint64_t, reg->m, reg->a);
  }

static inf_t *inf_stat(regHash_t *rghsh)
  {
  int i;
  khiter_t k;
  inf_t *inf = inf_init();
  for (k = kh_begin(rghsh); k < kh_end(rghsh); ++k)
    {
    if (kh_exist(rghsh, k))
      {
      bedreglist_t *bed = &kh_val(rghsh, k);
      for (i = 0; i < bed->m; ++i)
	{
	uint32_t beg = bed->a[i] >> 32;
	uint32_t end = (uint32_t)bed->a[i];
	assert(end >= beg);
	//uint32_t length = end - beg == 0 ? 1 : end - beg ; // 0based
	inf->length += end -beg;
	}

      inf->total += bed->m;
      }
    }
  return inf;
  }
/*
 * return the names of the chromosomes and the total number 
 * Inspired by Petr Danecek's regidx.c
 */
/* static char **seq_names(regHash_t *reghash, int *n) */
/*   { */
/*   int m = 2; */
/*   char **name = needmem(m, sizeof(char*)); */

/*   } */

/*
 * read bed file (0-based chromosome begin end) and tsv file (1-based chromosome posotion)
 */
static void
bed_read(const char *fn, regHash_t * reghash, int add1, int add2, int *ret) 
  {
  if (add1 < 0) errabort("the left flank size must be greater than 0! maybe you want try `trim` instead of `merge`...");
  if (add2 < 0) errabort("the right flank size must be greater than 0! maybe you want try `trim` instead of `merge`...");
  gzFile fp = safe_gzopen(fn);
  kstream_t * ks;
  int dret;
  kstring_t *str;
  str = (kstring_t*)needmem(sizeof(kstring_t));
  ks = ks_init(fp);
  int line = 0;
  while (ks_getuntil(ks, 0, str, &dret) >= 0)
    {
    int32_t beg = 0, end = 0;
    line++;
    if (isNull(str->s))
      {
      if (dret != '\n') while ((ks_getc(ks)) > 0 && dret != '\n');
      warnings("%s: line %d is empty! skip... ", fn, line);
      continue;
      }
    khiter_t k;
    k= kh_get(reg, reghash, str->s);
    if (k == kh_end(reghash))
      {
      int ret;
      k = kh_put(reg, reghash, strdup(str->s), &ret);
      bedreglist_t *b;
      b = (bedreglist_t*)needmem(sizeof(bedreglist_t));
      kh_val(reghash, k) = *b;
      }
    bedreglist_t *p = &kh_val(reghash, k);

    if (dret != '\n')
      {
      if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0]))
	{
	beg = atoi(str->s);
	if (dret != '\n')
	  {
	  if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0]))
	    {
	    end = atoi(str->s);	  
	    while (dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0); // skip all other parts
	    }
	  }
	}
      }
    if (beg < 0) beg = 0;
    if (end == 0) // treat as tsv file
      {
      // this is different from bed_read in bedidx.c
      end = beg;
      beg = beg < 1 ? 0 : beg-1;
      }
    
    if (end == 0)
      {
      warnings("%s: line %d is malformed! skip... ", fn, line);
      continue;
      }
    if (end < beg) swapvalue(beg, end, uint32_t);
    if (add1) beg -= add1;
    if (add2) end += add2; // end should be check the edge of the chromose
    if (beg <= 0) beg = 1;
    if (beg == end)
      {
      //warnings("begin == stop, only happened in insert variation. "
      //"make sure your bed file is 0-based! : %u\t%u", beg, end);
      *ret = 1;
      }
    if (p->m == p->n)
      {
      p->n = p->n ? p->n <<1 : 4;
      p->a = (uint64_t*)realloc(p->a, p->n*sizeof(uint64_t));
      }

    p->a[p->m++] = (uint64_t)beg<<32 | (uint32_t)end;
    //kh_val(reghash, k).a = p->a;
    }
  ks_destroy(ks);
  gzclose(fp);
  freemem(str->s); freemem(str); 
  }
	
static void bed_destroy(regHash_t *reghash, bedvoid_destroy func) 
  {
  khint_t k;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      bedreglist_t *bed = &kh_val(reghash, k);
      freemem(bed->a);
      freemem(bed->idx);
      if (bed->data) func(bed->data);
      freemem((char*)kh_key(reghash, k));
      kh_del(reg, reghash, k);
      //freemem(bed);
      }
    }
  kh_destroy(reg, reghash);
  }

void destroy_void(void *data)
{
}

static void reg_destroy(bedreglist_t *reg, bedvoid_destroy func) 
{
 freemem(reg->a);
 freemem(reg->idx);
 func(reg->data);
 freemem(reg);
}

static bedreglist_t * bed_directadd(bedreglist_t *bed1, bedreglist_t *bed2)
{
 bedreglist_t * bed;
 bed = calloc(1, sizeof(bedreglist_t));
 bed->n = bed1->m + bed2->m + 1;
 bed->m = bed1->m + bed2->m;
 bed->a = (uint64_t*)needmem( sizeof(uint64_t) * bed->n);
 memcpy(bed->a, bed1->a, bed1->m * sizeof(uint64_t));
 memcpy(bed->a + bed1->m, bed2->a, bed2->m * sizeof(uint64_t));
 return bed;
}

static void regcore_merge(bedreglist_t *bed) 
  {
  int i, m = 0;
  uint32_t lastbeg = 0, lastend = 0;
  if (bed->m == 1)
    {
    lastbeg = bed->a[0]>>32; lastend = (uint32_t)bed->a[0];
    bed->flag = 0;
    return;
    }
  sort_reg(bed);
  uint64_t *b;
  b = (uint64_t*) needmem((bed->m+1) * sizeof(uint64_t));
  for (i = 0; i < bed->m; ++i)
    {
    uint32_t beg, end;
    beg = bed->a[i]>>32; end = (uint32_t)bed->a[i];
    if (lastend < 1)
      {
      lastend = end;
      lastbeg = beg;
      continue;
      }
    if (lastend + 1> beg)
      {
      if (lastend < end) lastend = end;
      }
    else //if (lastend < beg)
      {
      b[m++] = (uint64_t) lastbeg<<32 |lastend;
      lastbeg = beg; lastend = end;
      }
    }
  if (lastend > 0)
    {
    b[m++] = (uint64_t) lastbeg<<32 | lastend;
    }
  memset(bed->a, 0, bed->m*sizeof(uint64_t));
  memcpy(bed->a, b, m*sizeof(uint64_t));
  bed->m = m;
  freemem(b);
  bed->flag = 0;
  }

static void bed_merge(regHash_t *reghash) 
  {
  khiter_t k;
  for (k = 0; k < kh_end(reghash); ++k)
    if (kh_exist(reghash, k)) regcore_merge(&kh_val(reghash, k));
  }

static bedreglist_t * regcore_uniq(bedreglist_t * bed1, bedreglist_t * bed2) 
  {
  int i= 0, j = 0;
  bedreglist_t * bed = bed_directadd(bed1, bed2);
  sort_reg(bed);
  uint32_t lastbeg = 0;
  uint32_t lastend = 0;
  uint64_t *b;
  b = (uint64_t*)needmem((bed->m+1) * sizeof(uint64_t));
  for (i = 0; i < bed->m; ++i)
    {
    int beg, end;
    beg = bed->a[i]>>32; end = (uint32_t)bed->a[i];

    /* condition 1:  init beg and end
     *
     *         lastbeg        lastend
     *         |              |
     *         ===============
     *                            ---------------------
     *                           |                    |
     *                           beg                  end
     *
     */
    if (lastend <= beg || lastend < 1)
      {
      lastbeg = beg; lastend = end;
      continue;
      }
    
    /* condition 2:  init beg and end
     *
     *         lastbeg                   lastend
     *         |                              |
     *         ================================
     *              ---------------------
     *             |                    |
     *             beg                  end
     *
     */
    
    if (lastend > end)
      {
      b[j++] = (uint64_t) beg <<32 | (uint32_t)end;
      lastbeg = end;
      continue;
      }
    else
      {
      /* condition 3:  init beg and end
       *
       *         lastbeg        lastend
       *         |              |
       *         ===============
       *                 ---------------------
       *                |                    |
       *                beg                  end
       *
       */
      b[j++] = (uint64_t) beg <<32 | (uint32_t)lastend;
      lastbeg = lastend;
      lastend = lastend == end ? 0 : end;
      }
    }
  bed->n = bed->m+1;
  bed->m = j;
  freemem(bed->a);
  bed->a = b;
  return bed;
  }

static bedreglist_t * regcore_diff(bedreglist_t * bed1, bedreglist_t * bed2)
  {
  int i, j = 0;
  bedreglist_t * uniq = regcore_uniq(bed1, bed2);
  regcore_merge(uniq);
  bedreglist_t * reg = bed_directadd(bed1, uniq);
  sort_reg(reg);
  uint64_t *b;
  b = (uint64_t*)needmem((reg->m+1) * sizeof(uint64_t));
  uint32_t lastbeg = 0; uint32_t lastend = 0;
  for (i = 0; i < reg->m; ++i)
    {
    uint32_t beg, end;
    beg = reg->a[i]>>32; end = (uint32_t)reg->a[i];

    /* condition 1:  init beg and end
     *
     *         beg             end
     *         |              |
     *         ===============
     *
     */
    if (lastend < 1)
      {
      lastbeg = beg; lastend = end;
      continue;
      }

    /* condition 2: 
     *
     *         lastbeg     lastend
     *         |          |
     *  region -----------
     *                 (no overlap)
     *                           ===============  new region
     *                          |              |
     *                          beg            end
     */
    if (lastend <= beg)
      {
      b[j++] = (uint64_t) lastbeg << 32| lastend;
      lastbeg = beg; lastend = end;
      continue;
      }
    if (lastbeg == beg)
      {
      
      /* condition 3: 
       *
       *         lastbeg     lastend
       *         |               |
       *  region ----------------
       *         |||||||||||||||| (equal)
       *         ===============  new region
       *        |              |
       *        beg            end
       */
      if (end == lastend)
	{
	lastbeg = lastend = 0;
	continue;
	}

      /* condition 4: 
       *
       *         lastbeg     lastend
       *         |          |
       *  region -----------
       *         | overlap |~~~~~~~
       *         ==================  new region
       *        |                 |
       *        beg            end
       */
      
      if (lastend < end)
	{
	//lastbeg = lastend + 1;
	lastbeg = lastend;
	lastend =  end;
	}
      else
	{
	// this condition will be not happened unless bed file is not sorted!
	/* condition 5: 
	 *
	 *         lastbeg                    lastend
	 *         |                          |
	 *  region ---------------------------
	 *         | overlap |||||||||
	 *         ==================  new region
	 *         |                |
	 *         beg            end
	 */
	errabort("condition5 error: please cantact the developer when you see this message!");
	lastbeg = end;
	}

      assert(lastbeg < lastend);
      continue;
      }
    
    // lastend > beg come here
    if (lastbeg < beg)
      {
      //b[j++] = (uint64_t) lastbeg << 32| (beg -1);
      b[j++] = (uint64_t) lastbeg << 32| beg ; // 0based beg
      }
    else
      {
      errabort("FIXME: not properly sorted!"
	       "Contact developer if you see this message!");
      }
    
    /* condition 6: 
     *
     *         lastbeg              lastend
     *         |                   |
     *  region --------------------
     *         ***********||||||||~~~~~~~~~~
     *                    ==================  new region
     *                   |                |
     *                   beg            end
     */
    if (lastend < end)
      {
      //lastbeg = lastend + 1;
      lastbeg = lastend;
      lastend = end;
      continue;
      }

    /* condition 7: 
     *
     *         lastbeg                 lastend
     *         |                            |
     *  region -----------------------------
     *         ***********||||||||||||||||||
     *                    ==================  new region
     *                   |                |
     *                   beg            end
     */
    else if (lastend == end)
      {
      lastbeg = lastend = 0;
      continue;
      }
    else
      {
      /* condition 8: 
       *
       *         lastbeg                                lastend
       *         |                                            |
       *  region ---------------------------------------------
       *         ***********||||||||||||||||||~~~~~~~~~~~~~~~~
       *                    ==================  new region
       *                   |                |
       *                   beg            end
       */
      lastbeg = end;
      }
    }
  reg->m = j;
  freemem(reg->a);
  reg->a = b;
  return reg;
  }

static bedreglist_t * regcore_trim(bedreglist_t *bed1, int trim1, int trim2)
  {
  regcore_merge(bed1);
  bedreglist_t *bed;
  bed = (bedreglist_t *)needmem(sizeof(bedreglist_t));
  bed->n = bed1->m+1;
  bed->a = (uint64_t*)needmem(bed->n * sizeof(uint64_t));
  bed->m = bed1->m;
  int i;
  for (i = 0; i < bed1->m; ++i)
    {
    uint32_t beg, end;
    beg = bed1->a[i] >> 32;
    end = (uint32_t)bed1->a[i];
	 
    if (end -beg <= trim1 + trim2)
      {
      warnings("[%u\t%u] The region is too short to trim! Skip... use base0to1 to trans 0-based file to 1-based file!", beg, end);
      }
    else
      {
      beg += trim1;
      end -= trim2;
      assert(beg < end);
      }
    bed->a[i] = (uint64_t)beg<<32 | end;
    }
  regcore_merge(bed);
  return bed;
  }

static void bed_uniq(regHash_t * reghash1, regHash_t * reghash2) 
  {
  khiter_t k, l;
  for (k = 0; k < kh_end(reghash1); ++k)
    {
    if (kh_exist(reghash1, k))
      {
      l = kh_get(reg, reghash2, kh_key(reghash1, k));
      if (l != kh_end(reghash2))
	{
	bedreglist_t * bed;
	bed = regcore_uniq( &kh_val(reghash1, k), &kh_val(reghash2, l));
	freemem(kh_val(reghash1, k).a);
	freemem(kh_val(reghash1, k).idx);
	memcpy(&kh_val(reghash1, k), bed, sizeof(bedreglist_t));
	freemem(bed);
	}
      else
	{
	freemem(kh_val(reghash1, k).a);
	freemem(kh_val(reghash1, k).idx);
	freemem((char*)kh_key(reghash1, k));
	kh_del(reg, reghash1, k);
	}
      }
    }
  }

static void bed_diff(regHash_t *reghash1, regHash_t *reghash2)
  {
  khiter_t k, l;
  for (k = 0; k < kh_end(reghash1); ++k)
    {
    if (kh_exist(reghash1, k))
      {
      l = kh_get(reg, reghash2, kh_key(reghash1, k));
      if (l != kh_end(reghash2))
	{
	bedreglist_t * bed;
	bed = regcore_diff(&kh_val(reghash1, k), &kh_val(reghash2, l));
	freemem(kh_val(reghash1, k).a);
	freemem(kh_val(reghash1, k).idx);
	memcpy(&kh_val(reghash1, k), bed, sizeof(bedreglist_t));
	freemem(bed);
	} 
      }
    }
  }


static void bed_trim(regHash_t *rghsh, int trim1, int trim2)
  {
  khiter_t k;
  for (k = 0; k < kh_end(rghsh); ++k)
    {
    if (kh_exist(rghsh, k))
      {
      bedreglist_t *bed = regcore_trim(&kh_val(rghsh, k), trim1, trim2);
      freemem(kh_val(rghsh, k).a);
      freemem(kh_val(rghsh, k).idx);
      memcpy(&kh_val(rghsh, k), bed, sizeof(bedreglist_t));
      freemem(bed);
      }
    }
  }

static void bed_view(regHash_t *reghash, char *reg)
  {
  bed_merge(reghash);
  
  }

static void bed_1to0(regHash_t *reghash)
  {
  khiter_t k;
  int i;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      for (i = 0; i < kh_val(reghash, k).m; ++i)
	{
	uint32_t beg = (kh_val(reghash, k).a[i] >>32) -1;
	uint32_t end = (uint32_t)kh_val(reghash, k).a[i];
	//if (isZero(beg)) errabort("the beg is 0: %u\t%u", beg, end);
	kh_val(reghash, k).a[i] = (uint64_t)beg <<32|end;
	}
      }
    }
  }

static void bed_0to1(regHash_t *reghash)
  {
  khiter_t k;
  int i;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      for (i = 0; i < kh_val(reghash, k).m; ++i)
	{
	uint32_t beg = (kh_val(reghash, k).a[i] >>32) +1;
	uint32_t end = (uint32_t)kh_val(reghash, k).a[i];
	if (beg > end) errabort("the beg is greater than end: %u\t%u ! it's not a 0 based file.", beg, end);
	kh_val(reghash, k).a[i] = (uint64_t)beg <<32|end;
	}
      }
    }
  }

static void bed_save(const char *fn, regHash_t *reghash) 
  {
  FILE *fp;
  fp = open_wfile(fn);
  khiter_t k;
  int i;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      for (i = 0; i < kh_val(reghash, k).m; ++i)
	{
	uint32_t beg = kh_val(reghash, k).a[i] >>32;
	uint32_t end = (uint32_t)kh_val(reghash, k).a[i];
	fprintf(fp, "%s\t%u\t%u\n", kh_key(reghash, k), beg, end);
	}
      }
    }
  fclose(fp);
  }

static void bed_pipeout(regHash_t *reghash)
  {
  khiter_t k;
  int i;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      for (i = 0; i < kh_val(reghash, k).m; ++i)
	{
	uint32_t beg = kh_val(reghash, k).a[i] >>32;
	uint32_t end = (uint32_t)kh_val(reghash, k).a[i];
	printf("%s\t%u\t%u\n", kh_key(reghash, k), beg, end);
	}
      }
    }
  fflush(stdout);
  }

static int bed_check_chromosome_length (regHash_t *reghash, chrHash_t *chrhash)
  {
  khiter_t k, l;
  int i;
  for (k = 0; k < kh_end(reghash); ++k)
    {
    if (kh_exist(reghash, k))
      {
      char *key = (char*)kh_key(reghash, k);
      l = kh_get(chr, chrhash, key);
      if (l == kh_end(chrhash))
	{
	warnings("There is no chromosome %s in this bam file!", key);
	continue;
	}
      uint32_t length = kh_val(chrhash, l).length;
      for (i = 0; i < kh_val(reghash, k).m; ++i)
	{
	uint32_t beg = kh_val(reghash, k).a[i] >>32;
	if (beg >= length)
	  {
	  kh_val(reghash, k).m = i;
	  return (int)l;
	  }
	uint32_t end = (uint32_t)kh_val(reghash, k).a[i];
	if (end > length)
	  {
	  kh_val(reghash, k).a[i] = (uint64_t)beg << 32 | length;
	  kh_val(reghash, k).m = i > 1 ? i -1 : 1;
	  return (int)l;
	  }
	}
      }
    }
  return -1;
  }

static bedHandle_t defaultBedHandler =
  {
  inf_init,
  bed_read,
  bed_merge,
  bed_uniq,
  bed_diff,
  bed_trim,
  bed_directadd,
  regcore_merge,
  regcore_uniq,
  regcore_diff,
  regcore_trim,
  sort_reg,
  bed_save,
  bed_destroy,
  inf_stat,
  bed_1to0,
  bed_0to1,
  bed_pipeout,
  bed_check_chromosome_length,
  };

bedHandle_t const *bedHand = &defaultBedHandler;
