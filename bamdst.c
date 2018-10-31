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

/* Author: Quan Shi (shiquan@genomics.cn) */

#include "commons.h"
#include "count.h"
#include "bedutil.h"

// bam.h and sam_header.h are standard header from samtools
#include "bam.h"
#include "sam_header.h"

// khash, kstring and knetfile are standard utils of klib
#include "khash.h" 
#include "kstring.h"
#include "knetfile.h"
#include "bgzf.h" // write tabix-able depth.gz file

static char const *program_name = "bamdst";
static char const *Version = "1.0.9";

/* flank region will be stat in the coverage report file,
 * this value can be set by -f / --flank */
static int flank_reg = 200;

/* extern bedHand from bedutil.c, it is a collection of functions*/
extern bedHandle_t *bedHand;

/* only accepted one stdin pipeline */
static bool stdin_lock = FALSE;

/* the bed file is zero based */
static bool zero_based = TRUE;

/* The number of threads after which there are 
   diminishing performance gains. */
//enum { DEFAULT_MAX_THREADS = 8 };

/* duplicate will be removed if rmdup_mark is TRUE 
 * this option is removed since version 1.0.0
 */
// static bool rmdup_mark = FALSE;

/* export target reads to a specitified bam file */
static char* export_target_bam = NULL;
bamFile bamoutfp;

int check_filename_isbam(char *name)
{
    int length = strlen(name);
    if (strcmp(name + length - 4, ".bam")) return 1;
    return 0;
}

static const int WINDOW_SIZE = 102400;
// force replace the existed files
//static bool is_forced = TRUE;

static char *outdir = NULL;

// init hash struct to store uncover regions
static int uncover_cutoff = 5;
static regHash_t *h_uncov;

void h_uncov_init()
{
    h_uncov = kh_init(reg);
}

// init hash struct to store chromosome length
static chrHash_t *h_chrlen;
void h_chrlength_init()
{
    h_chrlen = kh_init(chr);
}

// warp bamheader to retrieve chromosome legth hash
void header2chrhash(bam_header_t *h)
{
    int i, n, ret;
    khiter_t k;
    h->dict = sam_header_parse2(h->text);
    const char *tags[] = {"SN", "LN", "UR", "M5", NULL};
    char **tbl = sam_header2tbl_n(h->dict, "SQ", tags, &n);
    for (i = 0; i < n; i++)
    {
	k = kh_put(chr, h_chrlen, strdup(tbl[4*i]), &ret);
	kh_val(h_chrlen, k).length = atoi(tbl[4*i+1]);
    }
    if (tbl) free(tbl);
}

void chrhash_destroy()
{
    khiter_t k;
    if (h_chrlen)
    {
	for (k = 0; k < kh_end(h_chrlen); ++k)
	{
	    if (kh_exist(h_chrlen, k)) freemem((char*)kh_key(h_chrlen, k));
	    // skip destroy void data...
	}
	kh_destroy(chr, h_chrlen);
    }
}

/* @FLANK_REGION coverage of flank region list in report.
   @INSERTSIZE_LIMIT the insert size bigger than this will not be calculated.
   @MAXDEPTH_LIMIT the depth bigger than this will not be calculated, 
   0 means no limit.
   @CUTOFF  if you want know coverage of specity depth, set this value.
   @MAPQ_LIMIT  >=mapQ_limit list in report.
*/
struct opt_aux
{
    int nfiles;
    char **inputs;
    int isize_lim;
    int cutoff;
    int mapQ_lim;
    int maxdepth;
};


struct depnode
{
    unsigned len; // len will be set to 0, if not allocated memory
    unsigned start;
    unsigned stop;
    unsigned *vals; // raw depth
    unsigned *cnts; // reads count in this position
    unsigned *rmdupdep; // clean depth, rmdup, mapQ > 20, primary hit
    unsigned *covdep; // coverage depth
    struct depnode *next;
};

// debug the depnode list init
static const char * init_debugmsg[] =
    {
	"Success",
	"Trying to allocated an unempty node",
	"Trying to allocated a zero memory",
	"END of list"
    };

// debug macro, I think it is a good way to find memeory problem
#define INIT_DEBUG(x) do {						\
	int _a = x;							\
	if(_a > 0 && _a < 3)						\
	{								\
	    warnings("%s : %d %s",__FILE__, __LINE__, init_debugmsg[_a]); \
	}								\
    } while(0)

/* node->vals is the depth value of each loc 
 * node->cnts is the count of covered reads  
 * init the memory before use it */
static int depnode_init(struct depnode *node)
{
    if ( isNull(node) ) return 3;
    if ( node->len ) return 1;
    node->len = node->stop - node->start + 1;
    if ( isZero(node->len)) return 2;
    node->vals = (unsigned*)needmem((node->len) * sizeof(unsigned));
    node->cnts = (unsigned*)needmem((node->len) * sizeof(unsigned));
    node->rmdupdep = (unsigned*)needmem((node->len) * sizeof(unsigned));
    node->covdep = (unsigned*)needmem((node->len) * sizeof(unsigned));
    memset(node->vals, 0, node->len * sizeof(unsigned));
    memset(node->cnts, 0, node->len * sizeof(unsigned));
    memset(node->rmdupdep, 0, node->len * sizeof(unsigned));
    memset(node->covdep, 0, node->len * sizeof(unsigned));
    return 0;
}

/* delete node and make the node point to the next node */
#define del_node(node) do {			\
	if (node)				\
	{					\
	    struct depnode *tmpnode = node;	\
	    node = node->next;			\
	    freemem(tmpnode->vals);		\
	    freemem(tmpnode->cnts);		\
	    freemem(tmpnode->rmdupdep);		\
	    freemem(tmpnode->covdep);		\
	    freemem(tmpnode);			\
	}					\
    } while(0)

/* construct the bed struct array to a list , return the header node */
static struct depnode * bed_depnode_list(bedreglist_t *bed)
{
    struct depnode *node;
    struct depnode *header = NULL;
    struct depnode *tmpnode = NULL;
    int i;
    for ( i = 0; i < bed->m; ++i)
    {
	node = (struct depnode*)needmem(sizeof(struct depnode));
	node->start = (uint32_t)(bed->a[i] >> 32);
	node->stop = (uint32_t)bed->a[i];
	if (node->start < node->stop) node->start++; // 0-based to 1-based, sametimes inertion variation have same beg and end pos
    
	/* the length of this region should be zero if not allocated memory yet
	 *  Assign the length value when init the vals and cnts */
	node->len = 0;
	if ( isZero(i)) header = node;
	else tmpnode->next = node;
	tmpnode = node;
    }

    INIT_DEBUG(depnode_init(header));
    return header;
}

struct _aux
{
    /* nchr,  total num of chromsome
     * ndata, total num of bam struct
     * maxdepth, the maxmium of depth */
    int nchr, ndata, maxdep;

    /* tgt_len,  length of target region
     * flk_len,  length of flank region 
     * tgt_nreg, total target regions */
    uint64_t tgt_len, flk_len;
    unsigned tgt_nreg;

    /* data, array of bam struct
     * h,  point to bam header */
    bamFile *data;
    bam_header_t *h;

    /* h_tgt, target bed hash
     * h_flk, flank bed hash */
    regHash_t *h_tgt;
    regHash_t *h_flk;

    // count struct of depths, insertsize, flank depths, target regions
    count32_t *c_dep;
    count32_t *c_rmdupdep;
    count32_t *c_isize;
    count32_t *c_flkdep;
    count32_t *c_reg;
};

typedef struct _aux aux_t;

struct _aux  *aux_init()
{
    struct _aux *a;
    a = calloc(1, sizeof(struct _aux));
    a->nchr = a->maxdep = a->ndata = 0;
    a->data = NULL; // bamFile
    a->h = NULL; // bam header
    a->h_tgt = kh_init(reg);
    a->h_flk = kh_init(reg);
    count32_init(a->c_dep);
    count32_init(a->c_rmdupdep);
    count32_init(a->c_isize);
    count32_init(a->c_flkdep);
    count32_init(a->c_reg);
    return a;
}

void destroy_data(void *data)
{
    count32_t *cnt = (count32_t*)data;
    count_destroy(cnt);
}

void aux_destroy(struct _aux *a)
{
    free(a->data);
    bedHand->destroy((void *)a->h_tgt, destroy_data);
    bedHand->destroy((void *)a->h_flk, destroy_void);
    bam_header_destroy(a->h);
    if (a->c_dep->n > 0) count_destroy(a->c_dep);
    if (a->c_rmdupdep->n > 0) count_destroy(a->c_rmdupdep);
    if (a->c_flkdep->n > 0) count_destroy(a->c_flkdep);
    if (a->c_isize->n > 0) count_destroy(a->c_isize);
    if (a->c_reg->n > 0) count_destroy(a->c_reg);
    free(a);
}

typedef struct bamflag
{
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_dup, n_rmdup1, n_rmdup2;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_qcfail;
    uint64_t n_data, n_mdata;
    uint64_t n_qual;
    /* delete n_uniq
     * ref:  https://www.biostars.org/p/59281/ */
    uint64_t n_tgt, n_flk, n_tdata, n_fdata;
    uint64_t n_trmdat; // target rmdup data
} bamflag_t;


/*
 *  flag:
 *      0  qc failed
 *      1  clean read
 *      2  duplicate
 *      3  secondary alignment
 *      -1 normal end
 *      -2 trancate
 *      -3 unmap
 */
// use this macro to stat the flags 
#define flagstat(s, c, ret) do {					\
	++(s)->n_reads;							\
	(s)->n_data += (c)->l_qseq;					\
	if ((c)->flag & BAM_FQCFAIL)					\
	{								\
	    ++(s)->n_qcfail;						\
	    ret = 0;							\
	}								\
        else if ( (c)->flag & BAM_FSECONDARY || (c)->flag & BAM_FSUPPLEMENTARY) {\
            ret = 3;\
        }\
	else								\
	{								\
	    ret = 1;							\
	    if ((c)->flag & BAM_FPAIRED)				\
	    {								\
		++(s)->n_pair_all;					\
		if (((c)->flag & BAM_FPROPER_PAIR) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_pair_good; \
		if ((c)->flag & BAM_FREAD1) ++(s)->n_read1;		\
		if ((c)->flag & BAM_FREAD2) ++(s)->n_read2;		\
		if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) \
		    ++(s)->n_sgltn;					\
		if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP)) \
		{							\
		    ++(s)->n_pair_map;					\
		    if ((c)->mtid != (c)->tid) ++(s)->n_diffchr;	\
		}							\
	    }								\
	    if (!((c)->flag & BAM_FUNMAP))				\
	    {								\
		++(s)->n_mapped;					\
		(s)->n_mdata += (c)->l_qseq;				\
		if ((c)->flag & BAM_FREVERSE) ++(s)->n_mstrand;		\
		else ++(s)->n_pstrand;					\
		if ((c)->flag & BAM_FDUP)				\
		{							\
		    ++(s)->n_dup;					\
		    ret = 2;						\
		}							\
	    }								\
	    else {                                                      \
		ret = -3;						\
	    }								\
	}								\
    } while(0)

static void emit_try_help(void)
{
    fprintf (stderr, "out dir and bed file are mandatory!\n");
    fprintf (stderr, "Try '%s --help' for more information.\n", program_name);
}

void usage(int status)
{
    if (status == 0)
	emit_try_help();
    else
    {
	printf ("\n\
bamdst version: %s\n\
USAGE : %s [OPTION] -p <probe.bed> -o <output_dir> [in1.bam [in2.bam ... ]]\n\
   or : %s [OPTION] -p <probe.bed> -o <output_dir> -\n\
",
	    Version, program_name, program_name);
	puts ("\
Option -o and -p are mandatory:\n\
  -o, --outdir         output dir\n\
  -p, --bed            probe or target regions file, the region file will \n\
                       be merged before calculate depths\n\
");
	puts ("\
Optional parameters:\n\
   -f, --flank [200]   flank n bp of each region\n\
   -q [20]             map quality cutoff value, greater or equal to the value will be count\n\
   --maxdepth [0]      set the max depth to stat the cumu distribution.\n\
   --cutoffdepth [0]   list the coverage of above depths\n\
   --isize [2000]      stat the inferred insert size under this value\n\
   --uncover [5]       region will included in uncover file if below it\n\
   --bamout  BAMFILE   target reads will be exported to this bam file\n\
   -1                  begin position of bed file is 1-based\n\
   -h, --help          print this help info\n\
\n");
/*-d, --rmdup         remove dup reads when calculate depth\n	\*/
	    puts ("\
* Five essential files would be created in the output dir. \n\
* region.tsv.gz and depth.tsv.gz are zipped by bgzip, so you can use tabix \n\
  index these files.\n\n\
 - coverage.report     a report of the coverage information and reads \n\
                       information of whole target regions\n\
 - cumu.plot           distribution data of depth values\n\
 - insert.plot         distribution data of inferred insert size \n\
 - chromosome.report   coverage information for each chromosome\n\
 - region.tsv.gz       mean depth, median depth and coverage of each region\n\
 - depth.tsv.gz        raw depth, rmdup depth, coverage depth of each position\n\
 - uncover.bed         the bad covered or uncovered region in the probe file\n\
\n\
* About depth.tsv.gz:\n\
* There are five columns in this file, including chromosome, position, raw\n\
* depth, rmdep depth, coverage depth\n\
 - chromosome          the chromosome name\n\
 - position            1-based position of each chromosome\n\
 - raw depth           raw depth of position, not filter\n\
 - rmdup depth         remove duplication, and only calculate the reads which\n\
                       are primary mapped and mapQ >= cutoff_mapQ (default 20)\n\
 - coverage depth      calculate the deletions (CIGAR level) into depths,\n\
                       for coverage use.\n\
");
	    puts("============\n");
	    puts(" HOMEPAGE: \n\
      https://github.com/shiquan/bamdst\n");
    }
    exit(EXIT_SUCCESS);
}

#include"ksort.h"

KSORT_INIT_GENERIC(uint32_t)

static float median_cal(const uint32_t * array, int l)
{
    if (l == 0) return 0;
    if (l == 1) return (float)array[0];
    if (l == 2) return (float)(array[0] + array[1]) / 2;
    uint32_t *tmp;
    tmp = (uint32_t*)needmem(l *sizeof(uint32_t));
    memcpy(tmp, array, l * sizeof(uint32_t));
    ks_introsort(uint32_t, l, tmp);
    float med = l & 1 ? tmp[(l >> 1) + 1] :
	(float)(tmp[l >> 1] + tmp[(l >> 1) - 1]) / 2;
    mustfree(tmp);
    return med;
}

static float avg_cal(const uint32_t * array, int l)
{
    if (isNull(l))
	return 0;
    float avg = 0;
    int i;
    for (i = 0; i < l; ++i) avg += (float)array[i];
    avg /= (float)l;
    return avg ;
}

static float coverage_cal(const uint32_t * array, int l)
{
    if (isNull(l)) return 0;
    float cov = 0;
    int i;
    for (i = 0; i < l; ++i)
	if (array[i]) cov++;
    cov /= (float)l;
    return cov * 100;
}


// FIXME: need broken when bed file is truncated
int load_bed_init(char const *fn, aux_t * a)
{
    int ret = 0;
    bedHand->read(fn, a->h_tgt, 0, 0, &ret);
    if (zero_based && ret)
    {
	warnings("This region is not a standard bed format.\n"
		 "Please use parameter \"-1\" if your bed file is 1-based!");
    }
    if (!zero_based) bedHand->base1to0(a->h_tgt);
    bedHand->merge(a->h_tgt);
    int nchr;
    nchr = bedHand->check_length(a->h_tgt, h_chrlen);
    if (nchr != -1)
    {
	errabort(" The bed region is out the range of this chromosome %s. The max length is %u ..",
		 (char*)kh_key(h_chrlen, (khiter_t)nchr), kh_val(h_chrlen, (khiter_t)nchr).length);;
    }
    inf_t *inf1 = bedHand->stat(a->h_tgt);
    a->tgt_len = inf1->length;
    a->tgt_nreg = inf1->total;
  
    bedHand->read(fn, a->h_flk, flank_reg, flank_reg, &ret);
    if (!zero_based) bedHand->base1to0(a->h_flk);
    bedHand->merge(a->h_flk);
    bedHand->check_length(a->h_tgt, h_chrlen); // warp the length accord to the chromosome length
    bedHand->diff(a->h_flk, a->h_tgt);
    inf_t *inf2 = bedHand->stat(a->h_flk);
    a->flk_len = inf2->length;

    mustfree(inf1);
    mustfree(inf2);
    return 0;
}

// this function used to add an region to the bedregion struct
// use this struct to store the uncovered region
int push_bedreg(bedreglist_t *bed, uint32_t begin, uint32_t end)
{
    if (isZero(bed->n))
    {
	bed->n = 2;
	bed->a = (uint64_t*)needmem(bed->n * sizeof(uint64_t));
    }
    else if (bed->m == bed->n)
    {
	bed->n = bed->m << 1;
	bed->a = (uint64_t*)enlarge_empty_mem((void*)bed->a, bed->m * sizeof(uint64_t), bed->n *sizeof(uint64_t));
    }
    bed->a[bed->m++] = (uint64_t)begin << 32 | (uint32_t)end;
    return 0;
}

typedef enum
    {
	CMATCH,
	CDEL,
	CDUP,
        CLOWQ, // low quality, bug report by TOM MORRIS 
	UNKNOWN
    } cntstat_t;

int match_pos(struct depnode * header, uint32_t pos, cntstat_t state)
{
    struct depnode *tmp = header;
    while (tmp && pos > tmp->stop)
    {
	tmp = tmp->next;
    }

    if (isNull(tmp)) return 1; // this chromosome is finished, skip in next loop
    //debug("pos: %u\tstart: %u\tstop: %u", pos, tmp->start, tmp->stop);
    if (pos >= tmp->start)
    {
	if (isZero(tmp->len)) depnode_init(tmp);
	if (state == CMATCH)
	{
	    tmp->vals[pos - tmp->start]++;
	    tmp->rmdupdep[pos - tmp->start]++;
	    tmp->covdep[pos - tmp->start]++;
	}
	else if (state == CDEL)
	{
	    tmp->covdep[pos - tmp->start]++;
	}
	else 
	{
	    tmp->vals[pos - tmp->start]++;
	    tmp->covdep[pos - tmp->start]++;
	}
	return 0;
    }
    return 1; // not reachable
}

/* when deal with a read struct, check the begin of this read and the last position 
 * of this read, if this read is overlap with target region, sum up the depth of
 * each related position */
int readcore(struct depnode *header, bam1_t const * b, cntstat_t state)
{
    struct depnode *tmp = header;
    cntstat_t tmp_state = state;
    int i;
    bam1_core_t const *c = &b->core;
    int pos = c->pos + 1;

    /* while (tmp && pos > tmp->stop) */
    /*   { */
    /*   tmp= tmp->next; */
    /*   } */
    /* header = tmp; */
    if (isNull(tmp)) return 0;
    uint32_t *cigar = bam1_cigar(b);
    uint32_t end = bam_calend(c, cigar);

    if (end >= tmp->start)
    {
	int j = 0, l, s;
	if (pos >= tmp->start)
	{
	    if (isZero(tmp->len)) depnode_init(tmp);
	    tmp->cnts[pos - tmp->start]++;
	}
	for (i = 0; i < c->n_cigar; ++i)
	{
	    s = cigar[i] & 0xf;
	    l = cigar[i] >> BAM_CIGAR_SHIFT;
	    if (s == BAM_CDEL) tmp_state = CDEL;
	    else if (s == BAM_CMATCH) tmp_state = state;
	    else continue;
	    for (j = 0; j < l; ++j)
	    {
		if (pos >= tmp->start) match_pos(tmp, pos, tmp_state);
		pos++;
	    }
	}
	return 1;
    }
    return 0;
}

typedef struct
{
    int tid;
    int lstpos;
    bedreglist_t *tar;
    bedreglist_t *flk;
    bedreglist_t *ucreg;
    count32_t *depvals_of_chr;
    char *name;
    struct depnode *tgt_node;
    struct depnode *flk_node;
    kstring_t *pdepths;
    kstring_t *rcov;
    BGZF *fdep; // write depth to this file
    BGZF *freg; // write coverage of each region to this file
} loopbams_parameters_t;

loopbams_parameters_t * init_loopbams_parameters()
{
    loopbams_parameters_t * para;
    para = (loopbams_parameters_t*)needmem(sizeof(loopbams_parameters_t));
    *para = (loopbams_parameters_t){.tid=-1};
    para->pdepths = (kstring_t*)needmem(sizeof(kstring_t));
    para->rcov = (kstring_t*)needmem(sizeof(kstring_t));
    para->pdepths->l = para->pdepths->m = 0;
    para->rcov->l = para->rcov->m = 0;
    if (outdir) chdir(outdir);
    para->fdep = bgzf_open("depth.tsv.gz", "w");
    if (isNull(para->fdep))
	errabort("failed to open file depth.tsv.gz");
    para->freg = bgzf_open("region.tsv.gz", "w");
    if (isNull(para->freg))
	errabort("failed to open file region.tsv.gz");
    return para;
}

int close_loopbam_parameters(loopbams_parameters_t *para)
{
    if (para->tgt_node) {
	while(para->tgt_node)
	{
	    debug("length: %d\tstart:%u\tend:%u\n", para->tgt_node->len, para->tgt_node->start, para->tgt_node->stop);
	    para->tgt_node = para->tgt_node->next;
	}
	errabort("[close loopbam] target node is still reachable");
    }
    if (para->flk_node) errabort("[close loopbam] flank node is still reachable");
    freemem(para->pdepths->s);
    freemem(para->rcov->s);
    mustfree(para->pdepths);
    mustfree(para->rcov);
    bgzf_close(para->fdep);
    bgzf_close(para->freg);
    mustfree(para);
    return 0;
}

int write_buffer_bgzf(kstring_t *str, BGZF *fp)
{
    int write_size;
    if (str->l)
    {
	bgzf_flush_try(fp, str->l);
	write_size = bgzf_write(fp, str->s, str->l);
	str->l = 0;
	return write_size;
    }
    return 0;
}

int stat_each_region(loopbams_parameters_t *para, aux_t *a)
{
    struct depnode * node = para->tgt_node;
    if (isNull(node)) return 0;
    int j;
    float avg, med, cov1, cov2;
    uint32_t lst_start = 0; // uncover region start
    uint32_t lst_stop = 0; // uncover region stop
  
    if (node->len) {
	avg = avg_cal(node->vals, node->len);
	med = median_cal(node->vals, node->len);
	cov1 = coverage_cal(node->vals, node->len);
	cov2 = coverage_cal(node->covdep, node->len);
	for (j = 0; j < node->len; ++j) {
	    ksprintf(para->pdepths, "%s\t%d\t%u\t%u\t%u\n",
		     para->name, node->start + j, node->vals[j], node->rmdupdep[j], node->covdep[j]);
	    // count_increase will alloc memory space automatically
	    // use covdep to calculate coverage and averge depth
	    count_increase(para->depvals_of_chr, node->covdep[j], uint32_t);
	    count_increase(a->c_dep, node->covdep[j], uint32_t);
	    count_increase(a->c_rmdupdep, node->rmdupdep[j], uint32_t);
	    /* store the uncover region */
	    if (node->covdep[j] < uncover_cutoff) {
		if (isZero(lst_start) && isZero(lst_stop)) {
		    lst_start = node->start+j;
		    lst_stop = lst_start;
		} else {
		    lst_stop = node->start+j;
		}
	    } else if (lst_start > 0) {
		push_bedreg(para->ucreg, lst_start, lst_stop);
		lst_start = 0;
		lst_stop = 0;
	    }
	}
	if (lst_start > 0) {
	    push_bedreg(para->ucreg, lst_start, lst_stop);
	}
	write_buffer_bgzf(para->pdepths, para->fdep);
    } else {
	avg = med = cov1 = cov2 = 0.0;
	for (j = 0; j < node->len; ++j) {
	    ksprintf(para->pdepths, "%s\t%d\t0\t0\t0\n", para->name, node->start+j);
	    write_buffer_bgzf(para->pdepths, para->fdep);
	}
	push_bedreg(para->ucreg, node->start, node->stop); // store uncover region
	count_increaseN(para->depvals_of_chr, 0, node->len, uint32_t);
	count_increaseN(a->c_dep, 0, node->len, uint32_t);
    }
    //ksprintf(para->pdepths,"\n");
    count_increase(a->c_reg, (int)avg, uint32_t);
    ksprintf(para->rcov,"%s\t%u\t%u\t%.2f\t%.1f\t%.2f\t%.2f\n",
	     para->name, node->start-1, node->stop, avg, med, cov1, cov2);

    //if (para->rcov->l > WINDOW_SIZE)
    write_buffer_bgzf(para->rcov, para->freg);
    return 0;
}

// if bam files not contained all chromosomes in the bed file
int check_reachable_regions(loopbams_parameters_t *para, aux_t *a)
{
    khiter_t k;
    khiter_t l;
    int ret = 0;
    for (k = 0; k < kh_end(a->h_tgt); ++k) {
	if (kh_exist(a->h_tgt, k)) {
	    para->name = (char*)kh_key(a->h_tgt, k);
	    para->tar = &kh_val(a->h_tgt, k);
	    if (para->tar->flag == 1) continue; // already reach
	    warnings("%s is not contained in this bam file.", para->name);
	    count32_init(para->depvals_of_chr);
	    para->tar->data = (void*)para->depvals_of_chr;
	    para->flk = &kh_val(a->h_flk, k);
	    para->tgt_node = bed_depnode_list(para->tar);
	    para->flk_node = bed_depnode_list(para->flk);
	    l = kh_put(reg, h_uncov, strdup(para->name), &ret);
	    if (!ret) errabort("this chromosome should be empty! please contact developer to report this bug!");
	    bedreglist_t *ucreg_tmp;
	    ucreg_tmp = (bedreglist_t*)needmem(sizeof(bedreglist_t));
	    kh_val(h_uncov,l) = *ucreg_tmp;
	    para->ucreg = &kh_val(h_uncov, l);
      
	    while (para->tgt_node) {
		int length = para->tgt_node->stop - para->tgt_node->start +1;
		count_increaseN(a->c_dep, 0, length, uint32_t);
		stat_each_region(para, a);
		del_node(para->tgt_node); // no need allocate memory for these nodes
	    }
	    //count_merge(a->c_dep, para->depvals_of_chr, uint32_t);
	    while(para->flk_node) {
		int length = para->flk_node->stop - para->flk_node->start +1;
		count_increaseN(a->c_flkdep, 0, length, uint32_t);
		del_node(para->flk_node); // no need allocate memory for these nodes
	    }
	}
    }
    return 0;
}

int stat_flk_depcnt(loopbams_parameters_t *para, aux_t *a)
{
    int j;
    struct depnode *node = para->flk_node;  
    for (j = 0; j < node->len; ++j)
	count_increase(a->c_flkdep, node->vals[j], uint32_t);
    del_node(para->flk_node);
    depnode_init(para->flk_node);
    return 0;
}

void write_unover_file()
{
    if (outdir) chdir(outdir);
    bedHand->merge(h_uncov);
    bedHand->base1to0(h_uncov);
    bedHand->save("uncover.bed", h_uncov);
    bedHand->destroy(h_uncov, destroy_void);
}

// load bam files and stat the depths
// huge function, need IMPROVE it!!
int load_bamfiles(struct opt_aux *f, aux_t * a, bamflag_t * fs)
{
    // get the chromosome name from header
    bam_header_t *h = a->h;
    loopbams_parameters_t *para = init_loopbams_parameters();
    ksprintf(para->pdepths, "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth\n");
    ksprintf(para->rcov, "#Chr\tStart\tStop\tAvg depth\tMedian\tCoverage\tCoverage(FIX)\n");
    if (outdir) chdir(outdir); // FIXME: if there is no such dir?
    h_uncov_init();
    int i;
    for (i = 0; i < a->ndata; ++i)
    {
	bamFile dat = a->data[i];
	bool goto_next_chromosome = FALSE;
	int ret;
	cntstat_t state;
	// main loop
	while (1)
	{
	    bam1_t *b;
	    b = (bam1_t*)needmem(sizeof(bam1_t));
	    state = CMATCH;
	    ret = bam_read1(dat, b);
	    if (ret == -1) {
		mustfree(b); break; //normal end
	    } if (ret == -2) {
		errabort("%d bam file is truncated!\n", i + 1);
	    }
	    bam1_core_t *c = &b->core;
	    flagstat(fs, c, ret);

            // People usually want to know how many aligned reads with mapping quality greater than or equally to 10 or 20..
	    //if (c->qual > f->mapQ_lim) fs->n_qual++;            
            if (c->qual >= f->mapQ_lim) fs->n_qual++;
            else state = CLOWQ; 
            
	    if (c->tid == -1 || ret == -3) goto endcore; // unmapped~
            // secondary alignment
            if ( ret == 3 ) goto endcore;
            

	    if (ret == 2) {
		state = CDUP;
	    } else {
		if (c->flag & BAM_FREAD1) fs->n_rmdup1++;
		if (c->flag & BAM_FREAD2) fs->n_rmdup2++;
	    }


	    /* stat the insertsize */
	    if (c->isize > 0 && c->isize < f->isize_lim) {
		count_increase(a->c_isize, c->isize, uint32_t);
	    }

	    if (para->tid == c->tid ) {
		if (goto_next_chromosome) goto endcore;
		// Only accepted sorted bam files for effective
		if (para->lstpos > c->pos) {
		    errabort("The bam file is not sorted!");
		}
	    }
	    para->lstpos = c->pos; // 0-base

	    // if tid != c->tid, completed the stat of the last chromosome
	    // init the new chromosome node and clean the memory
	    // FIXME: need multi thread to improve it or NOT?
	    if (para->tid != c->tid) {
		goto_next_chromosome = FALSE; // clean the flag of skiping

		// still have reachable target node
		// release the left nodes
		if (para->tgt_node && para->tid >= 0)
		    while (para->tgt_node) {
			stat_each_region(para, a);
			del_node(para->tgt_node);
			// no need allocate memory space for these node, becase they are all 0s
		    }

		if (para->flk_node && para->tid >= 0) {
		    stat_flk_depcnt(para, a);
		    while (para->flk_node) {
			count_increaseN(a->c_flkdep, 0, para->flk_node->len, uint32_t);
			del_node(para->flk_node);
		    }
		}

		para->tid = c->tid;
		para->name = h->target_name[c->tid];

		// impossible in normal pratices, only happens in the different bam headers	
		if (para->tid > a->nchr) {
		    errabort("chromosome %s is not in bam header!"
			     "It must be use different bam headers"
			     , para->name);
		}
		khiter_t k;
		k = kh_get(reg, a->h_tgt, para->name);
		if (k == kh_end(a->h_tgt)) {
		    para->tgt_node = NULL;
		    para->flk_node = NULL;
		    goto_next_chromosome = TRUE;
		    goto endcore;
		}
		para->tar = &kh_val(a->h_tgt, k);
		para->flk = &kh_val(a->h_flk, k);
		count32_t *tmp;
		count32_init(tmp);
		para->depvals_of_chr = tmp;
		para->tar->data = (void*)para->depvals_of_chr;
		if (para->tar->flag || para->flk->flag) {
		    errabort("bam files are not properly sorted\n");
		}

		para->tgt_node = bed_depnode_list(para->tar);
		para->flk_node = bed_depnode_list(para->flk);
		para->tar->flag = para->flk->flag =1;

		/* the next part is init uncover region hash*/
		k = kh_put(reg, h_uncov, strdup(para->name), &ret);

		bedreglist_t *ucreg_tmp;
		ucreg_tmp = (bedreglist_t*)needmem(sizeof(bedreglist_t));
		kh_val(h_uncov, k) = *ucreg_tmp;
		para->ucreg = &kh_val(h_uncov, k); // FIXME: I don't understand why couldn't use ucreg_tmp directly
		/* finish init */
	    }
	    while (para->flk_node && para->flk_node->stop < para->lstpos+1) {
		stat_flk_depcnt(para, a);
	    }

	    if (para->flk_node && readcore(para->flk_node, b,state)) fs->n_flk++;

	    while (para->tgt_node && para->tgt_node->stop < para->lstpos +1) {
		stat_each_region(para, a);
		del_node(para->tgt_node);
		if (para->tgt_node && isZero(para->tgt_node->len)) depnode_init(para->tgt_node);
	    }
	    if (para->tgt_node && readcore(para->tgt_node, b,state)) {
		if (export_target_bam) bam_write1(bamoutfp, b);	
		fs->n_tgt++;
	    }

	  endcore:
	    bam_destroy1(b);
	}
	bgzf_close(a->data[i]);
    }
    while (para->tgt_node) {
	stat_each_region(para, a);
	del_node(para->tgt_node); // no need allocate memory for these nodes
    }
    while(para->flk_node) {
	count_increaseN(a->c_flkdep, 0, para->flk_node->len, uint32_t);
	del_node(para->flk_node); // no need allocate memory for these nodes
    }
    check_reachable_regions(para, a);
    write_unover_file();
    write_buffer_bgzf(para->pdepths, para->fdep);
    write_buffer_bgzf(para->rcov, para->freg);
    close_loopbam_parameters(para);
    return 0;
}
 
struct regcov
{
    uint64_t cnt, cnt4, cnt10, cnt30, cnt100, cntx;
    float    cov, cov4, cov10, cov30, cov100, covx;
};

struct regcov * regcov_init()
{
    struct regcov *c;
    c = (struct regcov*)needmem(sizeof(struct regcov));
    return c;
}
 
uint64_t cntcov_cal(struct opt_aux *f,
		    struct regcov * cov, count32_t * cnt, uint64_t *data)
{
    uint64_t rawcnt = 0;
    int i;
    *data = 0;
    *cov = (struct regcov){};
    for (i = 0; i < cnt->m; ++i)
    {
	(*data) += cnt->a[i] * i;
	rawcnt += cnt->a[i];
	if (i < 100) cov->cnt100 += cnt->a[i];
	if (i < 30)  cov->cnt30 += cnt->a[i];
	if (i < 10)  cov->cnt10 += cnt->a[i];
	if (i < 4) cov->cnt4 += cnt->a[i];
	if (f->cutoff && i < f->cutoff) cov->cntx += cnt->a[i];
    }
    if (rawcnt == 0) return 0;
    cov->cnt = rawcnt - (uint64_t)cnt->a[0];
    cov->cnt4 = rawcnt - cov->cnt4;
    cov->cnt10 = rawcnt - cov->cnt10;
    cov->cnt30 = rawcnt - cov->cnt30;
    cov->cnt100 = rawcnt - cov->cnt100;
    cov->cov = (float)cov->cnt / rawcnt * 100;
    cov->cov4 = (float)cov->cnt4 / rawcnt *100;
    cov->cov10 = (float)cov->cnt10 / rawcnt * 100;
    cov->cov30 = (float)cov->cnt30 / rawcnt *100;
    cov->cov100 = (float)cov->cnt100 / rawcnt*100;
    if (f->cutoff)
    {
	cov->cntx = rawcnt - cov->cntx;
	cov->covx = (float)cov->cntx / rawcnt*100;
    }
    return rawcnt;
}

// FIXME: need improve soon!!!
float median_cnt(count32_t *cnt)
{
    int i;
    uint64_t sum = 0;
    for (i = 0; i < cnt->m; ++i) sum += (uint64_t)cnt->a[i];
    uint64_t med = sum/2;
    uint64_t num = 0;
    for (i = 0; i < cnt->m; ++i)
    {
	num += (uint64_t)cnt->a[i];
	if (num >= med) return (float)i;
    }
    return 0;
}

int print_report(struct opt_aux *f, aux_t * a, bamflag_t * fs)
{
    int i;
    if (outdir) chdir(outdir);
    FILE *finsert;
    FILE *fdep;
    finsert = open_wfile("insertsize.plot");
    fdep = open_wfile("depth_distribution.plot");
  
    struct regcov *tarcov = regcov_init();
    struct regcov *flkcov = regcov_init();
    struct regcov *regcov = regcov_init();
    uint64_t icnt = 0;
    uint64_t dcnt = 0;
    for (i = 0; i < a->c_isize->m; ++i) icnt += a->c_isize->a[i];
    uint64_t icumu = icnt;
    for (i = 0; i < a->c_isize->m; ++i)
    {
	icumu -= a->c_isize->a[i];
	fprintf(finsert, "%d\t%u\t%f\t%"PRIu64"\t%f\n",
		i, a->c_isize->a[i], (float)a->c_isize->a[i] / icnt, icumu, (float)icumu/ icnt );
    }

    for (i = 0; i < a->c_dep->m; ++i) dcnt += a->c_dep->a[i];
    uint64_t dcumu = dcnt;
    for (i = 0; i < a->c_dep->m; ++i)
    {
	dcumu -= a->c_dep->a[i];
	fprintf(fdep, "%d\t%u\t%f\t%"PRIu64"\t%f\n",
		i, a->c_dep->a[i], (float)a->c_dep->a[i] / dcnt, dcumu, (float)dcumu/dcnt);
    }
    fclose(fdep);
    fclose(finsert);
    cntcov_cal(f, tarcov, a->c_dep, &fs->n_tdata);
    cntcov_cal(f, regcov, a->c_reg, &fs->n_fdata);
    //merge_cnt(a->c_flkdep, a->c_dep);
    cntcov_cal(f, flkcov, a->c_flkdep, &fs->n_fdata);
    for (i = 0; i < a->c_rmdupdep->m; ++i)
	fs->n_trmdat += a->c_rmdupdep->a[i] * i;

    FILE *fchrcov = open_wfile("chromosomes.report");
	{
	    fprintf(fchrcov, "%11s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s",
		    "#Chromosome","DATA(%)","Avg depth","Median","Coverage%","Cov 4x %","Cov 10x %","Cov 30x %","Cov 100x %");
	    if(f->cutoff) fprintf(fchrcov, "\tCov %dx", f->cutoff);
	    fprintf(fchrcov,"\n");
	    khiter_t k;
	    for (k = 0; k != kh_end(a->h_tgt); ++k)
	    {
		if (kh_exist(a->h_tgt, k))
		{
		    char *name = (char*)kh_key(a->h_tgt, k);
		    count32_t *cnt = (count32_t*)kh_val(a->h_tgt,k).data;
		    uint64_t data = 0;
		    struct regcov *chrcov = regcov_init();
		    cntcov_cal(f, chrcov, cnt, &data);
		    uint64_t length = 0;
		    int i;
		    for (i = 0; i < cnt->m; ++i) length += cnt->a[i];
		    float avg, med, per;
		    if (data > 0)
		    {
			avg = (float)data/ length;
			med = median_cnt(cnt);
			per = (float)data/fs->n_tdata*100.0;
			fprintf(fchrcov, "%11s\t%8.2f\t%8.2f\t%9.1f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f",
				name, per, avg, med, chrcov->cov, chrcov->cov4, chrcov->cov10, chrcov->cov30, chrcov->cov100);
			if (f->cutoff) fprintf(fchrcov, "\t%.2f", chrcov->covx);
		    }
		    else
		    {
			fprintf(fchrcov, "%11s\t%8.2f\t%8.2f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f",
				name, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0);
			if (f->cutoff) fprintf(fchrcov, "\t%5.4f", (float)0);
		    }
		    fprintf(fchrcov, "\n");
		}
	    }
	}
	fclose(fchrcov);
	FILE *fc = open_wfile("coverage.report");
	do
	{
	    fprintf(fc, "## The file was created by %s\n", program_name);
	    fprintf(fc, "## Version : %s\n", Version);
	    fprintf(fc, "## Files : ");
	    for (i = 0;  i < f->nfiles; ++i) fprintf(fc, "%s ", f->inputs[i]);
	    fprintf(fc, "\n");
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Raw Reads (All reads)", fs->n_reads); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] QC Fail reads", fs->n_qcfail); 
	    fprintf(fc, "%60s\t%.2f\n", "[Total] Raw Data(Mb)", (float)fs->n_data / 1e6); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Paired Reads", fs->n_pair_all);
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Mapped Reads", fs->n_mapped); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Mapped Reads", (float)fs->n_mapped / fs->n_reads *100); 
	    fprintf(fc, "%60s\t%.2f\n", "[Total] Mapped Data(Mb)", fs->n_mdata / 1e6); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Mapped Data(Mb)", (float)fs->n_mdata / fs->n_data *100); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Properly paired", fs->n_pair_good); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Properly paired", (float)fs->n_pair_good / fs->n_reads *100); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read and mate paired", fs->n_pair_map); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Read and mate paired", (float)fs->n_pair_map / fs->n_reads *100); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Singletons", fs->n_sgltn); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read and mate map to diff chr", fs->n_diffchr); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read1", fs->n_read1); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read2",fs->n_read2); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read1(rmdup)", fs->n_rmdup1); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] Read2(rmdup)",fs->n_rmdup2); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] forward strand reads", fs->n_pstrand); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] backward strand reads", fs->n_mstrand); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] PCR duplicate reads", fs->n_dup); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of PCR duplicate reads", (float)fs->n_dup / fs->n_mapped *100); // change n_reads to n_mapped, 2015/05/25
	    fprintf(fc, "%60s\t%d\n", "[Total] Map quality cutoff value", f->mapQ_lim); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Total] MapQuality above cutoff reads", fs->n_qual); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of MapQ reads in all reads", (float)fs->n_qual / fs->n_reads * 100); 
	    fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of MapQ reads in mapped reads", (float)fs->n_qual / fs->n_mapped *100); 
	    //tgt
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Target] Target Reads", fs->n_tgt);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Reads in all reads", (float)fs->n_tgt / fs->n_reads *100);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Reads in mapped reads", (float)fs->n_tgt / fs->n_mapped*100);
	    fprintf(fc, "%60s\t%.2f\n", "[Target] Target Data(Mb)", (float)fs->n_tdata / 1e6);
	    fprintf(fc, "%60s\t%.2f\n", "[Target] Target Data Rmdup(Mb)", (float)fs->n_trmdat / 1e6);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Data in all data", (float)fs->n_tdata / fs->n_data *100);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Data in mapped data", (float)fs->n_tdata / fs->n_mdata*100);
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Target] Len of region", a->tgt_len);
	    fprintf(fc, "%60s\t%.2f\n", "[Target] Average depth", (float)fs->n_tdata / a->tgt_len);
	    fprintf(fc, "%60s\t%.2f\n", "[Target] Average depth(rmdup)", (float)fs->n_trmdat / a->tgt_len);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>0x)", tarcov->cov);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=4x)", tarcov->cov4);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=10x)", tarcov->cov10);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=30x)", tarcov->cov30);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=100x)", tarcov->cov100);
	    if (f->cutoff)
	    {
		char titles[40];
		sprintf(titles,"[Target] Coverage (>=%ux)", f->cutoff);
		fprintf(fc, "%60s\t%.2f%%\n", titles, tarcov->covx);
	    }
	    //tgt regions
	    fprintf(fc, "%60s\t%u\n", "[Target] Target Region Count", a->tgt_nreg); 
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[Target] Region covered > 0x", regcov->cnt);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered > 0x", regcov->cov);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 4x", regcov->cov4);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 10x", regcov->cov10);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 30x", regcov->cov30);
	    fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 100x", regcov->cov100);
            if (f->cutoff)
	    {
		char titles[60];
                sprintf(titles,"[Target] Fraction Region covered (>=%ux)", f->cutoff);
		fprintf(fc, "%60s\t%.2f%%\n", titles, regcov->covx);
	    }
	    //tgt regions

	    //flk
	    fprintf(fc, "%60s\t%u\n", "[flank] flank size", flank_reg);
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[flank] Len of region (not include target region)", a->flk_len);
	    fprintf(fc, "%60s\t%.2f\n", "[flank] Average depth", (float)fs->n_fdata / a->flk_len);
	    fprintf(fc, "%60s\t%"PRIu64"\n", "[flank] flank Reads", fs->n_flk);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Reads in all reads", (float)fs->n_flk / fs->n_reads *100);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Reads in mapped reads", (float)fs->n_flk / fs->n_mapped*100);
	    fprintf(fc, "%60s\t%.2f\n", "[flank] flank Data(Mb)", (float)fs->n_fdata / 1e6);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Data in all data", (float)fs->n_fdata / fs->n_data *100);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Data in mapped data", (float)fs->n_fdata / fs->n_mdata*100);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>0x)", flkcov->cov);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=4x)", flkcov->cov4);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=10x)", flkcov->cov10);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=30x)", flkcov->cov30);
	    fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=100x)", flkcov->cov100);
            if (f->cutoff)
            {
		char titles[40];
                sprintf(titles,"[flank] Coverage (>=%ux)", f->cutoff);
		fprintf(fc, "%60s\t%.2f%%\n", titles, flkcov->covx);
	    }

	}
	while(0);
  
	mustfree(tarcov);
	mustfree(regcov);
	mustfree(flkcov);
	fclose(fc);
	return 0;
}

enum
    {
	MAXDEPTH,
	CUTOFF,
	INSERTSIZE,
	UNCOVER,
	BAMOUT,
	HELP
    };

static struct option const long_opts[] =
    {
	{"outdir", required_argument, NULL, 'o'},
	{"bedfile", required_argument, NULL, 'p'},
	{"flank", required_argument, NULL, 'f'},
	{"maxdepth", required_argument, NULL, MAXDEPTH},
	{"cutoffdepth", required_argument, NULL, CUTOFF},
	{"isize", required_argument, NULL, INSERTSIZE},
	{"mapthres", required_argument, NULL, 'q'},
	{"uncover", required_argument, NULL, UNCOVER},
	{"bamout", required_argument, NULL, BAMOUT},
	//{"rmdup", no_argument, NULL, 'd'},
	{"help", no_argument, NULL, 'h'},
	{"version", no_argument, NULL, 'v'}
    };

int show_version()
{
    printf("%s\n", Version);
    return 1;
}
int bamdst(int argc, char *argv[])
{
    int n, i;
    char *probe = 0;
  
    struct opt_aux opt = {.inputs=NULL, .isize_lim = 2000, .mapQ_lim = 20};
    while ((n = getopt_long(argc, argv, "o:p:f:q:l:h1v", long_opts, NULL)) >= 0)
    {
	switch (n)
	{
	    //output dir, must have right to write
	    case 'o': outdir = strdup(optarg); break;
		//capture region or just the region you interesting
	    case 'p': probe = strdup(optarg); break;
		//flk the region for more information, default is 200 bp
	    case 'f': flank_reg = atoi(optarg); break;
		//max depth to considered in the cumulative distribution of depths
	    case MAXDEPTH: opt.maxdepth = atoi(optarg); break;
	    case CUTOFF: opt.cutoff = atoi(optarg); break;
		// uncover_cutoff must be greater than 0
	    case UNCOVER: uncover_cutoff = atoi(optarg); assert(uncover_cutoff > 0); break;
	    case INSERTSIZE: opt.isize_lim = atoi(optarg); break;
	    case BAMOUT: export_target_bam = strdup(optarg); break;
	    case 'q': opt.mapQ_lim = atoi(optarg); break;
	    case 'h': usage(1); break;
	    case 'v': return show_version(); 
	    case '1': zero_based = FALSE; break;
		//case 'd': rmdup_mark = TRUE; break;
	    default: usage(0); 
		//more help
	}
		
    }
    if (isNull(outdir) || isNull(probe)) usage(0);
    if (export_target_bam && check_filename_isbam(export_target_bam))
    {
	fprintf(stderr,"--bamout must be a bam file: %s", export_target_bam);
	goto freeall;
    }
  
    n = argc - optind;
    //capable of deals with severl bam files
    aux_t * aux;
    aux = aux_init();
    if (isZero(n)) {
	aux->data = (bamFile*)needmem(sizeof(bamFile));
	aux->data[0] = bgzf_dopen(fileno(stdin), "r");
	aux->h = bam_header_read(aux->data[0]);
	aux->ndata = 1;
	opt.nfiles = 0;
    } else {
	aux->data = (bamFile*)needmem(n * sizeof(bamFile));
	opt.nfiles = n;
	opt.inputs = (char**)needmem(n * sizeof(char*));
	for (i = 0; i < n; ++i)
	{
	    bam_header_t *h_tmp;
	    h_tmp = calloc(1, sizeof(bam_header_t));
	    if ( STREQ(argv[optind + i], "-") )
	    {
		aux->data[i] = bgzf_dopen(fileno(stdin), "r");
		stdin_lock = 1;
	    }
	    else
	    {
		aux->data[i] = bgzf_open(argv[optind + i], "r");
	    }
	    if (aux->data[i] == NULL)
		errabort("%s: %s", argv[optind + i], strerror(errno));
	    h_tmp = bam_header_read(aux->data[i]);
	    if (i == 0) aux->h = h_tmp;
	    else bam_header_destroy(h_tmp);
	    opt.inputs[i] = strdup(argv[optind+i]);
	}
	aux->ndata = n;
    }
    // FIXME: accpet more than one bam files!
    if (export_target_bam)
    {
	bamoutfp = bam_open(export_target_bam, "w");
	if (bamoutfp == NULL)
	    errabort("%s : %s", export_target_bam, strerror(errno));
	bam_header_write(bamoutfp, aux->h);
    }
    h_chrlength_init();
    header2chrhash(aux->h);
    load_bed_init(probe, aux);
    chrhash_destroy();
    freemem(probe);
    aux->c_isize->a = calloc(opt.isize_lim, sizeof(unsigned));
    for (i = 0; i < opt.isize_lim; ++i) aux->c_isize->a[i] = 0;
    aux->c_isize->n = opt.isize_lim;
    aux->nchr = aux->h->n_targets;
    struct bamflag fs = {};
    load_bamfiles(&opt, aux, &fs);
    print_report(&opt, aux, &fs);
    aux_destroy(aux);
    for (i = 0; i < opt.nfiles; ++i) freemem(opt.inputs[i]);
    freemem(opt.inputs);
    if (export_target_bam) bam_close(bamoutfp);
  freeall:
    freemem(export_target_bam);
    freemem(outdir);
    return 0;
}

/* main */
int main(int argc, char *argv[])
{
    return bamdst(argc, argv);
}
