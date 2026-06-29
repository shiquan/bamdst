/* The MIT License

   Copyright (c) 2013-2026 The bamdst Authors

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

/*
 * bamdst — BAM Depth Statistics for Targeted Sequencing
 * ======================================================
 *
 * Architecture overview:
 *
 *   Input:  sorted BAM file(s) + probe BED file
 *   Output: coverage reports, per-base depth TSV, uncover regions
 *
 * Processing pipeline:
 *   1. Parse BED → chromosome-keyed hash of merged target intervals
 *      (h_tgt) and flank-extended intervals (h_flk).
 *   2. For each BAM read, walk its CIGAR to iterate over aligned bases.
 *      For each base within a target/flank interval, increment local
 *      depth counters in a linked list of depnode structs.
 *   3. As reads move past an interval (sorted BAM), finalize that
 *      interval: collect per-base depths into count32_t histograms,
 *      write depth.tsv.gz and region.tsv.gz rows, track uncover regions.
 *   4. After all reads: compute coverage statistics from histograms,
 *      print reports in the requested format (txt/csv/json).
 *
 * Key data types:
 *   depnode      — one target region interval, holds per-base depth arrays
 *   count32_t    — dynamic histogram: count[a[i]] = #positions with depth i
 *   regcov       — summary coverage metrics for a set of regions
 *   bamflag_t    — read-level flag statistics (mapped, paired, dup, etc.)
 *   opt_aux      — user-configurable options (thresholds, cutoffs, ratios)
 */

#include "bedutil.h"
#include "commons.h"
#include "count.h"

// bam.h and sam_header.h are standard header from samtools
#include "bam.h"
#include "sam_header.h"

// khash, kstring and knetfile are standard utils of klib
#include "bgzf.h" // write tabix-able depth.gz file
#include "khash.h"
#include "knetfile.h"
#include "kstring.h"

// GTF support (from PISA)
#include "gtf.h"

#include <sys/stat.h>
#include <pthread.h>

static char const *program_name = "bamdst";
static char const *Version = "2.0.0";

/* flank region will be stat in the coverage report file,
 * this value can be set by -f / --flank */
static int flank_reg = 200;

/* GTF annotation support */
static char *gtf_file = NULL;
static int gtf_level = GTF_LEVEL_GENE;  /* GTF_LEVEL_GENE (default) or GTF_LEVEL_EXON */
static bool gtf_enabled = FALSE;
static bool gtf_parallel = FALSE;       /* TRUE when -p and -g both used */
static bool gtf_stranded = FALSE;       /* TRUE when --gtf-stranded is set */

/* multi-threading */
static int n_threads = 1;               /* -@ / --threads, default 1 */
static bool mt_prepass_done = FALSE;    /* pre-pass counted global reads */

/* extern bedHand from bedutil.c, it is a collection of functions*/
extern bedHandle_t *bedHand;

/* output format: txt (default), csv, json */
typedef enum { FMT_TXT = 0, FMT_CSV, FMT_JSON } output_format_t;
static output_format_t output_fmt = FMT_TXT;

/* only accepted one stdin pipeline */
static bool stdin_lock = FALSE;

/* the bed file is zero based */
static bool zero_based = TRUE;

/* Maximum number of user-defined custom depth cutoffs supported (--cutoffdepth).
 * Each cutoff generates an additional coverage line in the report, e.g. --cutoffdepth 50,200
 * adds "Coverage (>=50x)" and "Coverage (>=200x)" to coverage.report. */
static const int MAX_CUTOFFS = 10;

/* The number of threads after which there are
   diminishing performance gains. */
// enum { DEFAULT_MAX_THREADS = 8 };

/* duplicate will be removed if rmdup_mark is TRUE
 * this option is removed since version 1.0.0
 */
// static bool rmdup_mark = FALSE;

/* export target reads to a specitified bam file */
static char *export_target_bam = NULL;
bamFile bamoutfp;

int check_filename_isbam(char *name)
{
    int length = strlen(name);
    if (strcmp(name + length - 4, ".bam"))
        return 1;
    return 0;
}

// static const int WINDOW_SIZE = 102400;
// force replace the existed files
// static bool is_forced = TRUE;

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
        k = kh_put(chr, h_chrlen, strdup(tbl[4 * i]), &ret);
        kh_val(h_chrlen, k).length = atoi(tbl[4 * i + 1]);
    }
    if (tbl)
        free(tbl);
}

void chrhash_destroy()
{
    khiter_t k;
    if (h_chrlen)
    {
        for (k = 0; k < kh_end(h_chrlen); ++k)
        {
            if (kh_exist(h_chrlen, k))
                freemem((char *)kh_key(h_chrlen, k));
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
    // int cutoff;
    bool cutoff;
    int num_cutoffs;
    int max_cutoffs; // allocated size of cutoffs[] array
    int *cutoffs;    // user-specified depth thresholds for coverage lines
    int mapQ_lim;
    int maxdepth;
    bool depth_ratio;   // enable ratio-based coverage stats (--depthratio)
    int num_ratios;     // number of user-specified ratios
    int max_ratios;     // allocated size of ratios[] array
    float *ratios;      // ratio thresholds, e.g. {0.1, 0.2, 0.5}
};

/*
 * Initialize the options structure with safe defaults.
 * - isize_lim=2000: ignore inferred insert sizes beyond this (visual clarity)
 * - mapQ_lim=20:   only reads with MAPQ >= 20 count as "high quality"
 * - depth_ratio defaults to FALSE for backward compatibility
 *   (user must explicitly enable via --depthratio)
 */
struct opt_aux init_opt_aux()
{
    struct opt_aux opt;
    memset(&opt, 0, sizeof(opt));
    opt.num_cutoffs = 0;
    opt.max_cutoffs = MAX_CUTOFFS;
    opt.cutoff = FALSE;
    opt.cutoffs = malloc(opt.max_cutoffs * sizeof(int));
    if (opt.cutoffs == NULL)
    {
        fprintf(stderr, "Memory allocation failed for cutoffs array\n");
        exit(EXIT_FAILURE);
    }
    opt.inputs = NULL;
    opt.isize_lim = 2000;
    opt.mapQ_lim = 20;
    /* Ratio-based coverage defaults to OFF for backward compatibility.
     * Enable with --depthratio 0.2,0.5 to get coverage at multiples of mean depth. */
    opt.depth_ratio = FALSE;
    opt.num_ratios = 0;
    opt.max_ratios = MAX_CUTOFFS;
    opt.ratios = malloc(opt.max_ratios * sizeof(float));
    if (opt.ratios == NULL)
    {
        fprintf(stderr, "Memory allocation failed for ratios\n");
        exit(EXIT_FAILURE);
    }
    return opt;
}

/*
 * A linked list node representing one target region interval.
 *
 * Each node covers positions [start, stop] on a chromosome.
 * The four per-base arrays are allocated lazily (len == 0 until depnode_init).
 *
 *   vals[]     — raw depth (all reads, no filtering)
 *   cnts[]     — number of distinct reads covering this position
 *   rmdupdep[] — filtered depth: excludes PCR duplicates, secondary alignments,
 *                and reads with MAPQ < mapQ_lim (default 20). Only primary hits.
 *   covdep[]   — coverage depth: like raw depth but also counts deletion (CDEL)
 *                positions. Used for coverage calculations since deletions affect
 *                whether a position is "covered" by a read.
 *
 * Nodes are freed as reads pass beyond their stop position (sorted BAM required).
 */
struct depnode
{
    unsigned len;     // allocated length (stop-start+1), 0 = not yet allocated
    unsigned start;   // 1-based start position of this region
    unsigned stop;    // 1-based end position
    unsigned *vals;     // raw depth per position
    unsigned *cnts;     // distinct read count per position
    unsigned *rmdupdep; // clean depth (no dup, MAPQ-filtered, primary only)
    unsigned *covdep;   // coverage depth (raw + CDEL positions)
    struct depnode *next;
};

// debug the depnode list init
static const char *init_debugmsg[] = {"Success", "Trying to allocated an unempty node",
                                      "Trying to allocated a zero memory", "END of list"};

// debug macro, I think it is a good way to find memeory problem
#define INIT_DEBUG(x)                                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        int _a = x;                                                                                                    \
        if (_a > 0 && _a < 3)                                                                                          \
        {                                                                                                              \
            warnings("%s : %d %s", __FILE__, __LINE__, init_debugmsg[_a]);                                             \
        }                                                                                                              \
    } while (0)

/* node->vals is the depth value of each loc
 * node->cnts is the count of covered reads
 * init the memory before use it */
static int depnode_init(struct depnode *node)
{
    if (isNull(node))
        return 3;
    if (node->len)
        return 1;
    node->len = node->stop - node->start + 1;
    if (isZero(node->len))
        return 2;
    node->vals = (unsigned *)needmem((node->len) * sizeof(unsigned));
    node->cnts = (unsigned *)needmem((node->len) * sizeof(unsigned));
    node->rmdupdep = (unsigned *)needmem((node->len) * sizeof(unsigned));
    node->covdep = (unsigned *)needmem((node->len) * sizeof(unsigned));
    memset(node->vals, 0, node->len * sizeof(unsigned));
    memset(node->cnts, 0, node->len * sizeof(unsigned));
    memset(node->rmdupdep, 0, node->len * sizeof(unsigned));
    memset(node->covdep, 0, node->len * sizeof(unsigned));
    return 0;
}

/* delete node and make the node point to the next node */
#define del_node(node)                                                                                                 \
    do                                                                                                                 \
    {                                                                                                                  \
        if (node)                                                                                                      \
        {                                                                                                              \
            struct depnode *tmpnode = node;                                                                            \
            node = node->next;                                                                                         \
            freemem(tmpnode->vals);                                                                                    \
            freemem(tmpnode->cnts);                                                                                    \
            freemem(tmpnode->rmdupdep);                                                                                \
            freemem(tmpnode->covdep);                                                                                  \
            freemem(tmpnode);                                                                                          \
        }                                                                                                              \
    } while (0)

/* construct the bed struct array to a list , return the header node */
static struct depnode *bed_depnode_list(bedreglist_t *bed)
{
    struct depnode *node;
    struct depnode *header = NULL;
    struct depnode *tmpnode = NULL;
    int i;
    for (i = 0; i < bed->m; ++i)
    {
        node = (struct depnode *)needmem(sizeof(struct depnode));
        node->start = (uint32_t)(bed->a[i] >> 32);
        node->stop = (uint32_t)bed->a[i];
        if (node->start < node->stop)
            node->start++; // 0-based to 1-based, sametimes inertion variation have
                           // same beg and end pos

        /* the length of this region should be zero if not allocated memory yet
         *  Assign the length value when init the vals and cnts */
        node->len = 0;
        if (isZero(i))
            header = node;
        else
            tmpnode->next = node;
        tmpnode = node;
    }

    INIT_DEBUG(depnode_init(header));
    return header;
}

/* Stranded GTF auxiliary data (forward or reverse). */
struct gtf_strand_aux {
    regHash_t *h_tgt;
    regHash_t *h_flk;
    uint64_t len;
    unsigned nreg;
    count32_t *c_dep;
    count32_t *c_rmdupdep;
    count32_t *c_reg;
};

#define GTF_S_INIT {NULL, NULL, 0, 0, NULL, NULL, NULL}

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

    /* GTF regions (when -g is used alongside -p, stored separately) */
    regHash_t *h_gtf;
    regHash_t *h_gtf_flk;
    uint64_t gtf_len;
    unsigned gtf_nreg;

    // count struct of depths, insertsize, flank depths, target regions
    count32_t *c_dep;
    count32_t *c_rmdupdep;
    count32_t *c_isize;
    count32_t *c_flkdep;
    count32_t *c_reg;

    /* GTF depth histograms (non-stranded, or combined) */
    count32_t *c_gtf_dep;
    count32_t *c_gtf_rmdupdep;
    count32_t *c_gtf_reg;

    /* Strand-separated GTF data (when --gtf-stranded is used) */
    struct gtf_strand_aux gtf_fwd;
    struct gtf_strand_aux gtf_rev;
};

typedef struct _aux aux_t;

struct _aux *aux_init()
{
    struct _aux *a;
    a = calloc(1, sizeof(struct _aux));
    a->nchr = a->maxdep = a->ndata = 0;
    a->data = NULL; // bamFile
    a->h = NULL;    // bam header
    a->h_tgt = kh_init(reg);
    a->h_flk = kh_init(reg);
    a->h_gtf = NULL;
    a->h_gtf_flk = NULL;
    count32_init(a->c_dep);
    count32_init(a->c_rmdupdep);
    count32_init(a->c_isize);
    count32_init(a->c_flkdep);
    count32_init(a->c_reg);
    /* GTF histograms are initialised lazily (when -g is used) */
    a->c_gtf_dep = NULL;
    a->c_gtf_rmdupdep = NULL;
    a->c_gtf_reg = NULL;
    a->gtf_fwd = (struct gtf_strand_aux)GTF_S_INIT;
    a->gtf_rev = (struct gtf_strand_aux)GTF_S_INIT;
    return a;
}

void destroy_data(void *data)
{
    count32_t *cnt = (count32_t *)data;
    count_destroy(cnt);
}

void aux_destroy(struct _aux *a)
{
    free(a->data);
    bedHand->destroy((void *)a->h_tgt, destroy_data);
    bedHand->destroy((void *)a->h_flk, destroy_void);
    if (a->h_gtf)      bedHand->destroy((void *)a->h_gtf, destroy_void);
    if (a->h_gtf_flk)  bedHand->destroy((void *)a->h_gtf_flk, destroy_void);
    bam_header_destroy(a->h);
    /* if (a->c_dep->n > 0) count_destroy(a->c_dep); */
    /* if (a->c_rmdupdep->n > 0) count_destroy(a->c_rmdupdep); */
    /* if (a->c_flkdep->n > 0) count_destroy(a->c_flkdep); */
    /* if (a->c_isize->n > 0) count_destroy(a->c_isize); */
    /* if (a->c_reg->n > 0) count_destroy(a->c_reg); */
    count_destroy(a->c_dep);
    count_destroy(a->c_rmdupdep);
    count_destroy(a->c_flkdep);
    count_destroy(a->c_isize);
    count_destroy(a->c_reg);
    if (a->c_gtf_dep)      count_destroy(a->c_gtf_dep);
    if (a->c_gtf_rmdupdep) count_destroy(a->c_gtf_rmdupdep);
    if (a->c_gtf_reg)      count_destroy(a->c_gtf_reg);
    if (a->gtf_fwd.h_tgt)     bedHand->destroy((void *)a->gtf_fwd.h_tgt, destroy_void);
    if (a->gtf_fwd.h_flk)     bedHand->destroy((void *)a->gtf_fwd.h_flk, destroy_void);
    if (a->gtf_rev.h_tgt)     bedHand->destroy((void *)a->gtf_rev.h_tgt, destroy_void);
    if (a->gtf_rev.h_flk)     bedHand->destroy((void *)a->gtf_rev.h_flk, destroy_void);
    if (a->gtf_fwd.c_dep)      count_destroy(a->gtf_fwd.c_dep);
    if (a->gtf_fwd.c_rmdupdep) count_destroy(a->gtf_fwd.c_rmdupdep);
    if (a->gtf_fwd.c_reg)      count_destroy(a->gtf_fwd.c_reg);
    if (a->gtf_rev.c_dep)      count_destroy(a->gtf_rev.c_dep);
    if (a->gtf_rev.c_rmdupdep) count_destroy(a->gtf_rev.c_rmdupdep);
    if (a->gtf_rev.c_reg)      count_destroy(a->gtf_rev.c_reg);

    free(a);
}

/*
 * Read-level statistics counters.
 *
 * These are accumulated across all BAM files and printed in coverage.report.
 * n_data / n_mdata track total / mapped base counts (for Mb conversion).
 * n_tgt / n_flk count reads overlapping target or flank regions.
 * n_trmdat stores total rmdup base count (for rmdup average depth).
 */
typedef struct bamflag
{
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_dup, n_rmdup1, n_rmdup2;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_qcfail;
    uint64_t n_data, n_mdata;
    uint64_t n_qual;
    /* n_uniq was removed — see https://www.biostars.org/p/59281/ */
    uint64_t n_tgt, n_flk, n_tdata, n_fdata;
    uint64_t n_trmdat; // total rmdup bases in target regions
    /* GTF parallel stats */
    uint64_t n_gtf_reads, n_gtf_data;
    uint64_t n_gtf_rmdpdat;
    /* GTF strand-separated stats */
    uint64_t n_gtf_fwd_reads, n_gtf_rev_reads;
    uint64_t n_gtf_fwd_data, n_gtf_rev_data;
    uint64_t n_gtf_fwd_rmdpdat, n_gtf_rev_rmdpdat;
} bamflag_t;

/*
 * Classify a BAM read and update the flag statistics counters.
 *
 * Sets `ret` to one of:
 *   0  — QC failed (BAM_FQCFAIL set)
 *   1  — clean read (passes QC, not duplicate)
 *   2  — PCR/optical duplicate (BAM_FDUP set)
 *  -1  — normal end (not returned by this macro, used in caller)
 *  -2  — truncated file (not returned by this macro)
 *  -3  — unmapped read (BAM_FUNMAP set)
 *
 * The caller uses `ret` to decide the depth-counting behaviour:
 *   ret=1 (clean) → increment both raw and rmdup depths
 *   ret=2 (dup)  → increment raw depth only
 *   ret=0 (QC fail) / ret=-3 (unmapped) → skip depth counting entirely
 */
#define flagstat(s, c, ret)                                                                                            \
    do                                                                                                                 \
    {                                                                                                                  \
        ++(s)->n_reads;                                                                                                \
        (s)->n_data += (c)->l_qseq;                                                                                    \
        if ((c)->flag & BAM_FQCFAIL)                                                                                   \
        {                                                                                                              \
            ++(s)->n_qcfail;                                                                                           \
            ret = 0;                                                                                                   \
        }                                                                                                              \
        else                                                                                                           \
        {                                                                                                              \
            ret = 1;                                                                                                   \
            if ((c)->flag & BAM_FPAIRED)                                                                               \
            {                                                                                                          \
                ++(s)->n_pair_all;                                                                                     \
                if (((c)->flag & BAM_FPROPER_PAIR) && !((c)->flag & BAM_FUNMAP))                                       \
                    ++(s)->n_pair_good;                                                                                \
                if ((c)->flag & BAM_FREAD1)                                                                            \
                    ++(s)->n_read1;                                                                                    \
                if ((c)->flag & BAM_FREAD2)                                                                            \
                    ++(s)->n_read2;                                                                                    \
                if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP))                                            \
                    ++(s)->n_sgltn;                                                                                    \
                if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP))                                           \
                {                                                                                                      \
                    ++(s)->n_pair_map;                                                                                 \
                    if ((c)->mtid != (c)->tid)                                                                         \
                        ++(s)->n_diffchr;                                                                              \
                }                                                                                                      \
            }                                                                                                          \
            if (!((c)->flag & BAM_FUNMAP))                                                                             \
            {                                                                                                          \
                ++(s)->n_mapped;                                                                                       \
                (s)->n_mdata += (c)->l_qseq;                                                                           \
                if ((c)->flag & BAM_FREVERSE)                                                                          \
                    ++(s)->n_mstrand;                                                                                  \
                else                                                                                                   \
                    ++(s)->n_pstrand;                                                                                  \
                if ((c)->flag & BAM_FDUP)                                                                              \
                {                                                                                                      \
                    ++(s)->n_dup;                                                                                      \
                    ret = 2;                                                                                           \
                }                                                                                                      \
            }                                                                                                          \
            else                                                                                                       \
            {                                                                                                          \
                ret = -3;                                                                                              \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

static void emit_try_help(void)
{
    fprintf(stderr, "out dir and bed file are mandatory!\n");
    fprintf(stderr, "Try '%s --help' for more information.\n", program_name);
}

void usage(int status)
{
    if (status == 0)
        emit_try_help();
    else
    {
        printf("\n\
bamdst version: %s\n\
USAGE : %s [OPTION] -p <probe.bed> -o <output_dir> [in1.bam [in2.bam ... ]]\n\
   or : %s [OPTION] -p <probe.bed> -o <output_dir> -\n\
",
               Version, program_name, program_name);
        puts("\
Option -o is mandatory; at least one of -p or -g is required:\n\
  -o, --outdir         output dir\n\
  -p, --bed            probe or target regions file, the region file will \n\
                       be merged before calculate depths\n\
  -g, --gtf            GTF annotation file (gene/exon coverage)\n\
");
        puts("\
Optional parameters:\n\
   -f, --flank [200]   flank n bp of each region\n\
   --gtf-level [gene]  GTF coverage level: gene or exon (default: gene)\n\
   --gtf-stranded      separate forward/reverse strand GTF coverage\n\
   -@, --threads [1]   number of threads (requires .bai index)\n\
   -q [20]             map quality cutoff value, greater or equal to the value will be count\n\
   --maxdepth [0]      set the max depth to stat the cumu distribution.\n\
   --cutoffdepth [0,0] list the coverage of above these depths, allow maximal 10 cutoffs.\n\
   --isize [2000]      stat the inferred insert size under this value\n\
   --uncover [5]       region will included in uncover file if below it\n\
   --bamout  BAMFILE   target reads will be exported to this bam file\n\
   -F, --format [txt]  output format: txt, csv, json (default: txt)\n\
   --depthratio [0,0]  coverage at ratios of avg depth (e.g., 0.2,0.5)\n\
   -1                  begin position of bed file is 1-based\n\
   -h, --help          print this help info\n\
\n");

        puts("\
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
    exit(status == 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}



#include "ksort.h"

KSORT_INIT_GENERIC(uint32_t)

static float median_cal(const uint32_t *array, int l)
{
    if (l == 0)
        return 0;
    if (l == 1)
        return (float)array[0];
    if (l == 2)
        return (float)(array[0] + array[1]) / 2;
    uint32_t *tmp;
    tmp = (uint32_t *)needmem(l * sizeof(uint32_t));
    memcpy(tmp, array, l * sizeof(uint32_t));
    ks_introsort(uint32_t, l, tmp);
    float med = l & 1 ? tmp[l >> 1] : (float)(tmp[l >> 1] + tmp[(l >> 1) - 1]) / 2;
    mustfree(tmp);
    return med;
}

static float avg_cal(const uint32_t *array, int l)
{
    if (isNull(l))
        return 0;
    float avg = 0;
    int i;
    for (i = 0; i < l; ++i)
        avg += (float)array[i];
    avg /= (float)l;
    return avg;
}

static float coverage_cal(const uint32_t *array, int l)
{
    if (isNull(l))
        return 0;
    float cov = 0;
    int i;
    for (i = 0; i < l; ++i)
        if (array[i])
            cov++;
    cov /= (float)l;
    return cov * 100;
}

// FIXME: need broken when bed file is truncated
int load_bed_init(char const *fn, aux_t *a)
{
    int ret = 0;
    bedHand->read(fn, a->h_tgt, 0, 0, &ret);
    if (zero_based && ret)
    {
        warnings("This region is not a standard bed format.\n"
                 "Please use parameter \"-1\" if your bed file is 1-based!");
    }
    if (!zero_based)
        bedHand->base1to0(a->h_tgt);
    bedHand->merge(a->h_tgt);
    int nchr;
    nchr = bedHand->check_length(a->h_tgt, h_chrlen);
    if (nchr != -1)
    {
        errabort(" The bed region is out the range of this chromosome %s. The max "
                 "length is %u ..",
                 (char *)kh_key(h_chrlen, (khiter_t)nchr), kh_val(h_chrlen, (khiter_t)nchr).length);
        ;
    }
    inf_t *inf1 = bedHand->stat(a->h_tgt);
    a->tgt_len = inf1->length;
    a->tgt_nreg = inf1->total;

    bedHand->read(fn, a->h_flk, flank_reg, flank_reg, &ret);
    if (!zero_based)
        bedHand->base1to0(a->h_flk);
    bedHand->merge(a->h_flk);
    bedHand->check_length(a->h_tgt, h_chrlen); // warp the length accord to the chromosome length
    // bedHand->diff(a->h_flk, a->h_tgt);
    inf_t *inf2 = bedHand->stat(a->h_flk);
    a->flk_len = inf2->length;

    mustfree(inf1);
    mustfree(inf2);
    return 0;
}

/*
 * Load GTF annotation and convert to region hashes for coverage analysis.
 *
 * Parses the GTF file, extracts gene or exon intervals (depending on
 * gtf_level), converts them to 0-based regHash_t, and merges overlapping
 * intervals.  The resulting hash is stored in a->h_gtf (and a->h_gtf_flk
 * for flank-extended intervals).
 *
 * When used alongside -p (BED), GTF regions are tracked as a separate
 * parallel target set in load_bamfiles().
 * When used alone (no -p), GTF regions become the primary target.
 */
int load_gtf_init(char const *fn, aux_t *a)
{
    fprintf(stderr, "[bamdst] Loading GTF: %s (level=%s)\n",
              fn, gtf_level == GTF_LEVEL_EXON ? "exon" : "gene");

    struct gtf_spec *G = gtf_read_lite(fn);
    if (G == NULL)
        errabort("Failed to load GTF file: %s", fn);

    a->h_gtf = gtf_to_regHash(G, gtf_level);
    if (a->h_gtf == NULL)
        errabort("No regions extracted from GTF: %s", fn);

    /* merge overlapping intervals (genes may overlap on same strand,
     * and we treat all strands together) */
    bedHand->merge(a->h_gtf);

    /* create flank-extended version */
    a->h_gtf_flk = kh_init(reg);
    {
        khiter_t k;
        for (k = 0; k != kh_end(a->h_gtf); ++k) {
            if (kh_exist(a->h_gtf, k)) {
                const char *chr = kh_key(a->h_gtf, k);
                bedreglist_t *src = &kh_val(a->h_gtf, k);
                int ret2;
                khiter_t k2 = kh_put(reg, a->h_gtf_flk, strdup(chr), &ret2);
                bedreglist_t *dst = &kh_val(a->h_gtf_flk, k2);
                if (ret2) {
                    memset(dst, 0, sizeof(bedreglist_t));
                }
                int i;
                for (i = 0; i < src->m; ++i) {
                    uint32_t beg = (uint32_t)(src->a[i] >> 32);
                    uint32_t end = (uint32_t)src->a[i];
                    /* extend by flank_reg on both sides, clamp to 0 */
                    uint32_t ext_beg = beg > (uint32_t)flank_reg
                                     ? beg - (uint32_t)flank_reg : 0;
                    uint32_t ext_end = end + (uint32_t)flank_reg;
                    /* push extended interval */
                    if (dst->n == 0) {
                        dst->n = 2;
                        dst->a = (uint64_t *)needmem(dst->n * sizeof(uint64_t));
                    } else if (dst->m == dst->n) {
                        dst->n = dst->n + 1024;
                        dst->a = (uint64_t *)
                            enlarge_empty_mem(dst->a,
                                dst->m * sizeof(uint64_t),
                                dst->n * sizeof(uint64_t));
                    }
                    dst->a[dst->m++] = (uint64_t)ext_beg << 32
                                     | (uint32_t)ext_end;
                }
            }
        }
        bedHand->merge(a->h_gtf_flk);
    }

    /* record GTF region statistics (non-stranded combined) */
    {
        inf_t *inf = bedHand->stat(a->h_gtf);
        a->gtf_len = inf->length;
        a->gtf_nreg = inf->total;
        mustfree(inf);
    }

    /* init GTF depth histograms (non-stranded) */
    count32_init(a->c_gtf_dep);
    count32_init(a->c_gtf_rmdupdep);
    count32_init(a->c_gtf_reg);

    /* Strand-separated GTF regions */
    if (gtf_stranded) {
        regHash_t *h_fwd = NULL, *h_rev = NULL;
        gtf_to_regHash_stranded(G, gtf_level, &h_fwd, &h_rev);

        bedHand->merge(h_fwd);
        bedHand->merge(h_rev);

        a->gtf_fwd.h_tgt = h_fwd;
        a->gtf_rev.h_tgt = h_rev;

        inf_t *ifwd = bedHand->stat(h_fwd);
        inf_t *irev = bedHand->stat(h_rev);
        a->gtf_fwd.len = ifwd->length;
        a->gtf_fwd.nreg = ifwd->total;
        a->gtf_rev.len = irev->length;
        a->gtf_rev.nreg = irev->total;
        mustfree(ifwd);
        mustfree(irev);

        count32_init(a->gtf_fwd.c_dep);
        count32_init(a->gtf_fwd.c_rmdupdep);
        count32_init(a->gtf_fwd.c_reg);
        count32_init(a->gtf_rev.c_dep);
        count32_init(a->gtf_rev.c_rmdupdep);
        count32_init(a->gtf_rev.c_reg);

        fprintf(stderr, "[bamdst] GTF stranded: fwd %u regions/%" PRIu64 " bp, "
                "rev %u regions/%" PRIu64 " bp\n",
                a->gtf_fwd.nreg, a->gtf_fwd.len,
                a->gtf_rev.nreg, a->gtf_rev.len);
    }

    fprintf(stderr, "[bamdst] GTF: %u regions, %" PRIu64 " bp\n",
              a->gtf_nreg, a->gtf_len);

    gtf_destroy(G);
    return 0;
}

// this function used to add an region to the bedregion struct
// use this struct to store the uncovered region
int push_bedreg(bedreglist_t *bed, uint32_t begin, uint32_t end)
{
    if (isZero(bed->n))
    {
        bed->n = 2;
        bed->a = (uint64_t *)needmem(bed->n * sizeof(uint64_t));
    }
    else if (bed->m == bed->n)
    {
        bed->n = bed->n + 1024;
        bed->a = (uint64_t *)enlarge_empty_mem((void *)bed->a, bed->m * sizeof(uint64_t), bed->n * sizeof(uint64_t));
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

int match_pos(struct depnode *header, uint32_t pos, cntstat_t state)
{
    struct depnode *tmp = header;
    while (tmp && pos > tmp->stop)
    {
        tmp = tmp->next;
    }

    if (isNull(tmp))
        return 1; // this chromosome is finished, skip in next loop
    // debug("pos: %u\tstart: %u\tstop: %u", pos, tmp->start, tmp->stop);
    if (pos >= tmp->start)
    {
        if (isZero(tmp->len))
            depnode_init(tmp);
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

/* when deal with a read struct, check the begin of this read and the last
 * position of this read, if this read is overlap with target region, sum up the
 * depth of each related position */
int readcore(struct depnode *header, bam1_t const *b, cntstat_t state)
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
    if (isNull(tmp))
        return 0;
    uint32_t *cigar = bam1_cigar(b);
    uint32_t end = bam_calend(c, cigar);

    if (end >= tmp->start)
    {
        int j = 0, l, s;
        if (pos >= tmp->start)
        {
            if (isZero(tmp->len))
                depnode_init(tmp);
            tmp->cnts[pos - tmp->start]++;
        }
        for (i = 0; i < c->n_cigar; ++i)
        {
            s = cigar[i] & 0xf;
            l = cigar[i] >> BAM_CIGAR_SHIFT;
            if (s == BAM_CDEL)
                tmp_state = CDEL;
            // else if (s == BAM_CMATCH) tmp_state = state;
            else if (s == BAM_CMATCH || s == BAM_CEQUAL || s == BAM_CDIFF)
                tmp_state = state;
            else
                continue;
            for (j = 0; j < l; ++j)
            {
                if (pos >= tmp->start)
                    match_pos(tmp, pos, tmp_state);
                pos++;
            }
        }
        return 1;
    }
    return 0;
}

/* Stranded GTF per-loop data. */
struct gtf_strand_para {
    bedreglist_t *tar;
    struct depnode *tgt_node;
    count32_t *depvals_of_chr;
};

#define GTF_SP_INIT {NULL, NULL, NULL}

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
    /* GTF parallel tracking (non-stranded) */
    bedreglist_t *gtf_tar;
    struct depnode *gtf_tgt_node;
    count32_t *gtf_depvals_of_chr;
    /* GTF strand-separated tracking */
    struct gtf_strand_para gtf_fwd;
    struct gtf_strand_para gtf_rev;
    kstring_t *pdepths;
    kstring_t *rcov;
    BGZF *fdep; // write depth to this file
    BGZF *freg; // write coverage of each region to this file
} loopbams_parameters_t;

loopbams_parameters_t *init_loopbams_parameters()
{
    loopbams_parameters_t *para;
    para = (loopbams_parameters_t *)needmem(sizeof(loopbams_parameters_t));
    *para = (loopbams_parameters_t){.tid = -1};
    para->gtf_tar = NULL;
    para->gtf_tgt_node = NULL;
    para->gtf_depvals_of_chr = NULL;
    para->gtf_fwd = (struct gtf_strand_para)GTF_SP_INIT;
    para->gtf_rev = (struct gtf_strand_para)GTF_SP_INIT;
    para->pdepths = (kstring_t *)needmem(sizeof(kstring_t));
    para->rcov = (kstring_t *)needmem(sizeof(kstring_t));
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
    if (para->tgt_node)
    {
        while (para->tgt_node)
        {
            debug("length: %d\tstart:%u\tend:%u\n", para->tgt_node->len, para->tgt_node->start, para->tgt_node->stop);
            para->tgt_node = para->tgt_node->next;
        }
        errabort("[close loopbam] target node is still reachable");
    }
    if (para->flk_node)
        errabort("[close loopbam] flank node is still reachable");
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
    struct depnode *node = para->tgt_node;
    if (isNull(node))
        return 0;
    int j;
    float avg, med, cov1, cov2;
    uint32_t lst_start = 0; // uncover region start
    uint32_t lst_stop = 0;  // uncover region stop

    if (node->len)
    {
        avg = avg_cal(node->vals, node->len);
        med = median_cal(node->vals, node->len);
        cov1 = coverage_cal(node->vals, node->len);
        cov2 = coverage_cal(node->covdep, node->len);
        for (j = 0; j < node->len; ++j)
        {
            ksprintf(para->pdepths, "%s\t%d\t%u\t%u\t%u\n", para->name, node->start + j, node->vals[j],
                     node->rmdupdep[j], node->covdep[j]);
            // count_increase will alloc memory space automatically
            // use covdep to calculate coverage and averge depth
            count_increase(para->depvals_of_chr, node->covdep[j], uint32_t);
            count_increase(a->c_dep, node->covdep[j], uint32_t);
            count_increase(a->c_rmdupdep, node->rmdupdep[j], uint32_t);
            /* store the uncover region */
            if (node->covdep[j] < uncover_cutoff)
            {
                if (isZero(lst_start) && isZero(lst_stop))
                {
                    lst_start = node->start + j;
                    lst_stop = lst_start;
                }
                else
                {
                    lst_stop = node->start + j;
                }
            }
            else if (lst_start > 0)
            {
                push_bedreg(para->ucreg, lst_start, lst_stop);
                lst_start = 0;
                lst_stop = 0;
            }
        }
        if (lst_start > 0)
        {
            push_bedreg(para->ucreg, lst_start, lst_stop);
        }
        write_buffer_bgzf(para->pdepths, para->fdep);
    }
    else
    {
        avg = med = cov1 = cov2 = 0.0;
        for (j = 0; j < node->len; ++j)
        {
            ksprintf(para->pdepths, "%s\t%d\t0\t0\t0\n", para->name, node->start + j);
            write_buffer_bgzf(para->pdepths, para->fdep);
        }
        push_bedreg(para->ucreg, node->start, node->stop); // store uncover region
        count_increaseN(para->depvals_of_chr, 0, node->len, uint32_t);
        count_increaseN(a->c_dep, 0, node->len, uint32_t);
    }
    // ksprintf(para->pdepths,"\n");
    count_increase(a->c_reg, (int)avg, uint32_t);
    ksprintf(para->rcov, "%s\t%u\t%u\t%.2f\t%.1f\t%.2f\t%.2f\n", para->name, node->start - 1, node->stop, avg, med,
             cov1, cov2);

    // if (para->rcov->l > WINDOW_SIZE)
    write_buffer_bgzf(para->rcov, para->freg);
    return 0;
}

// if bam files not contained all chromosomes in the bed file
int check_reachable_regions(loopbams_parameters_t *para, aux_t *a)
{
    khiter_t k;
    khiter_t l;
    int ret = 0;
    for (k = 0; k < kh_end(a->h_tgt); ++k)
    {
        if (kh_exist(a->h_tgt, k))
        {
            para->name = (char *)kh_key(a->h_tgt, k);
            para->tar = &kh_val(a->h_tgt, k);
            if (para->tar->flag == 1)
                continue; // already reach
            warnings("%s is not contained in this bam file.", para->name);
            count32_init(para->depvals_of_chr);
            para->tar->data = (void *)para->depvals_of_chr;
            para->flk = &kh_val(a->h_flk, k);
            para->tgt_node = bed_depnode_list(para->tar);
            para->flk_node = bed_depnode_list(para->flk);
            l = kh_put(reg, h_uncov, strdup(para->name), &ret);
            if (!ret)
                errabort("this chromosome should be empty! please contact developer to "
                         "report this bug!");
            bedreglist_t ucreg_tmp = {0};
            kh_val(h_uncov, l) = ucreg_tmp;
            para->ucreg = &kh_val(h_uncov, l);

            while (para->tgt_node)
            {
                int length = para->tgt_node->stop - para->tgt_node->start + 1;
                count_increaseN(a->c_dep, 0, length, uint32_t);
                stat_each_region(para, a);
                del_node(para->tgt_node); // no need allocate memory for these nodes
            }
            // count_merge(a->c_dep, para->depvals_of_chr, uint32_t);
            while (para->flk_node)
            {
                int length = para->flk_node->stop - para->flk_node->start + 1;
                count_increaseN(a->c_flkdep, 0, length, uint32_t);
                del_node(para->flk_node); // no need allocate memory for these nodes
            }
        }
    }
    return 0;
}

/*
 * Like stat_each_region, but for GTF regions tracked in parallel.
 * Only accumulates depth histograms (no per-base TSV for GTF).
 */
int stat_each_region_gtf(loopbams_parameters_t *para, aux_t *a)
{
    struct depnode *node = para->gtf_tgt_node;
    if (isNull(node) || isNull(a->c_gtf_dep))
        return 0;
    int j;
    if (node->len) {
        for (j = 0; j < node->len; ++j) {
            count_increase(a->c_gtf_dep, node->covdep[j], uint32_t);
            count_increase(a->c_gtf_rmdupdep, node->rmdupdep[j], uint32_t);
            if (para->gtf_depvals_of_chr)
                count_increase(para->gtf_depvals_of_chr, node->covdep[j], uint32_t);
        }
    } else {
        int length = node->stop - node->start + 1;
        count_increaseN(a->c_gtf_dep, 0, length, uint32_t);
        if (para->gtf_depvals_of_chr)
            count_increaseN(para->gtf_depvals_of_chr, 0, length, uint32_t);
    }
    /* record region-level mean depth for GTF */
    if (a->c_gtf_reg) {
        float avg = node->len ? avg_cal(node->vals, node->len) : 0.0f;
        count_increase(a->c_gtf_reg, (int)avg, uint32_t);
    }
    return 0;
}

/*
 * Accumulate depth from a depnode into a gtf_strand_aux histogram set.
 * Node is passed explicitly (fwd or rev), stats go into 'gs'.
 */
static int gtf_s_depnode_stat(struct depnode *node, struct gtf_strand_aux *gs)
{
    if (isNull(node) || isNull(gs->c_dep))
        return 0;
    int j;
    if (node->len) {
        for (j = 0; j < node->len; ++j) {
            count_increase(gs->c_dep, node->covdep[j], uint32_t);
            count_increase(gs->c_rmdupdep, node->rmdupdep[j], uint32_t);
        }
    } else {
        int length = node->stop - node->start + 1;
        count_increaseN(gs->c_dep, 0, length, uint32_t);
    }
    if (gs->c_reg && node->len) {
        float avg = avg_cal(node->vals, node->len);
        count_increase(gs->c_reg, (int)avg, uint32_t);
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
    return 0;
}

void write_unover_file()
{
    if (outdir) chdir(outdir);
    bedHand->merge(h_uncov);
    bedHand->base1to0(h_uncov);
    bedHand->save("uncover.bed", h_uncov);
    bedHand->destroy(h_uncov, destroy_data);
}

/* ================================================================
 * Multi-threading support
 * ================================================================ */

/* Merge count32_t histogram. */
static void mt_count_merge(count32_t *dst, count32_t *src)
{
    if (!src || src->m == 0) return;
    int i, max_m = dst->m > src->m ? dst->m : src->m;
    count_resize(dst, max_m, uint32_t);
    for (i = 0; i < src->m; ++i) dst->a[i] += src->a[i];
    if (src->m > dst->m) dst->m = src->m;
}

/* Merge bamflag_t: only target-specific fields (global counts come from pre-pass). */
static void mt_flag_merge(bamflag_t *dst, const bamflag_t *src)
{
    if (mt_prepass_done) {
        dst->n_tgt        += src->n_tgt;
        dst->n_flk        += src->n_flk;
        dst->n_tdata      += src->n_tdata;
        dst->n_fdata      += src->n_fdata;
        dst->n_trmdat     += src->n_trmdat;
        dst->n_gtf_reads  += src->n_gtf_reads;
        dst->n_gtf_data   += src->n_gtf_data;
        dst->n_gtf_rmdpdat += src->n_gtf_rmdpdat;
        return;
    }
    dst->n_reads      += src->n_reads;
    dst->n_mapped     += src->n_mapped;
    dst->n_pair_map   += src->n_pair_map;
    dst->n_pair_all   += src->n_pair_all;
    dst->n_pair_good  += src->n_pair_good;
    dst->n_sgltn      += src->n_sgltn;
    dst->n_read1      += src->n_read1;
    dst->n_read2      += src->n_read2;
    dst->n_dup        += src->n_dup;
    dst->n_rmdup1     += src->n_rmdup1;
    dst->n_rmdup2     += src->n_rmdup2;
    dst->n_diffchr    += src->n_diffchr;
    dst->n_pstrand    += src->n_pstrand;
    dst->n_mstrand    += src->n_mstrand;
    dst->n_qcfail     += src->n_qcfail;
    dst->n_data       += src->n_data;
    dst->n_mdata      += src->n_mdata;
    dst->n_qual       += src->n_qual;
    dst->n_tgt        += src->n_tgt;
    dst->n_flk        += src->n_flk;
    dst->n_tdata      += src->n_tdata;
    dst->n_fdata      += src->n_fdata;
    dst->n_trmdat     += src->n_trmdat;
    dst->n_gtf_reads  += src->n_gtf_reads;
    dst->n_gtf_data   += src->n_gtf_data;
    dst->n_gtf_rmdpdat += src->n_gtf_rmdpdat;
}

/* ---- per-thread context ---- */
struct mt_ctx {
    pthread_t handle;
    int tid;
    const char *bam_path;
    const bam_header_t *hdr;
    regHash_t *h_tgt, *h_flk, *h_gtf;
    int8_t *chrom_mask;
    int n_chrom_total, mapQ_lim, isize_lim;

    /* results */
    count32_t *c_dep, *c_rmdupdep, *c_isize, *c_flkdep, *c_reg;
    bamflag_t fs;
    count32_t *c_gtf_dep, *c_gtf_rmdupdep;
};

/* ---- simplified depth stat (no BGZF output, just histogram) ---- */
static void mts_stat_region(struct depnode *node,
    count32_t *c_dep, count32_t *c_rmdupdep,
    count32_t *c_reg, count32_t *depvals_chr)
{
    int j;
    if (node->len) {
        for (j = 0; j < node->len; ++j) {
            count_increase(c_dep, node->covdep[j], uint32_t);
            count_increase(c_rmdupdep, node->rmdupdep[j], uint32_t);
            if (depvals_chr)
                count_increase(depvals_chr, node->covdep[j], uint32_t);
        }
        count_increase(c_reg, (int)avg_cal(node->vals, node->len), uint32_t);
    } else {
        int len = node->stop - node->start + 1;
        count_increaseN(c_dep, 0, len, uint32_t);
        if (depvals_chr)
            count_increaseN(depvals_chr, 0, len, uint32_t);
    }
}
static void mts_stat_flk(struct depnode *node, count32_t *c_flkdep)
{
    int j;
    for (j = 0; j < node->len; ++j)
        count_increase(c_flkdep, node->vals[j], uint32_t);
}
static void mts_stat_gtf(struct depnode *node,
    count32_t *c_dep, count32_t *c_rmdupdep)
{
    int j;
    if (node->len) {
        for (j = 0; j < node->len; ++j) {
            count_increase(c_dep, node->covdep[j], uint32_t);
            count_increase(c_rmdupdep, node->rmdupdep[j], uint32_t);
        }
    } else {
        int len = node->stop - node->start + 1;
        count_increaseN(c_dep, 0, len, uint32_t);
    }
}

/* ---- thread worker: sequential BAM scan with chromosome mask ---- */
static void *mt_worker(void *arg)
{
    struct mt_ctx *ctx = (struct mt_ctx *)arg;
    if (!ctx->bam_path) return NULL;

    bamFile fp = bam_open(ctx->bam_path, "r");
    if (!fp) {
        warnings("[thread %d] cannot open %s", ctx->tid, ctx->bam_path);
        return NULL;
    }
    { bam_header_t *wh = bam_header_read(fp); bam_header_destroy(wh); }

    count32_init(ctx->c_dep); count32_init(ctx->c_rmdupdep);
    count32_init(ctx->c_isize); count32_init(ctx->c_flkdep);
    count32_init(ctx->c_reg);
    memset(&ctx->fs, 0, sizeof(bamflag_t));
    if (ctx->h_gtf) {
        count32_init(ctx->c_gtf_dep);
        count32_init(ctx->c_gtf_rmdupdep);
    }

    struct depnode *tgt_node = NULL, *flk_node = NULL, *gtf_node = NULL;
    count32_t *depvals_chr = NULL;
    bedreglist_t *tar = NULL;
    int cur_tid = -1;

    bam1_t *b = (bam1_t *)needmem(sizeof(bam1_t));

    while (bam_read1(fp, b) >= 0) {
        bam1_core_t *c = &b->core;
        if (c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY)
            continue;

        int sf;
        flagstat(&ctx->fs, c, sf);
        if (c->qual >= ctx->mapQ_lim) ctx->fs.n_qual++;
        if (c->tid == -1 || sf == -3) continue;

        cntstat_t state = CMATCH;
        if (sf == 2) state = CDUP;
        else if (c->qual < ctx->mapQ_lim) state = CLOWQ;
        if (sf == 1) {
            if (c->flag & BAM_FREAD1) ctx->fs.n_rmdup1++;
            if (c->flag & BAM_FREAD2) ctx->fs.n_rmdup2++;
        }
        if (c->isize > 0 && c->isize < ctx->isize_lim)
            count_increase(ctx->c_isize, c->isize, uint32_t);

        /* Chromosome transition */
        if (c->tid != cur_tid) {
            while (tgt_node) {
                mts_stat_region(tgt_node, ctx->c_dep, ctx->c_rmdupdep,
                                ctx->c_reg, depvals_chr);
                del_node(tgt_node);
            }
            while (flk_node) { mts_stat_flk(flk_node, ctx->c_flkdep); del_node(flk_node); }
            while (gtf_node) { mts_stat_gtf(gtf_node, ctx->c_gtf_dep, ctx->c_gtf_rmdupdep); del_node(gtf_node); }

            cur_tid = c->tid;
            tgt_node = NULL; flk_node = NULL; gtf_node = NULL;
            tar = NULL; depvals_chr = NULL;

            if (cur_tid >= ctx->n_chrom_total || !ctx->chrom_mask[cur_tid])
                continue;

            const char *chr_name = ctx->hdr->target_name[cur_tid];
            khiter_t k = kh_get(reg, ctx->h_tgt, chr_name);
            if (k == kh_end(ctx->h_tgt)) continue;

            tar = &kh_val(ctx->h_tgt, k);
            tgt_node = bed_depnode_list(tar);
            count32_init(depvals_chr);
            tar->data = (void *)depvals_chr;

            k = kh_get(reg, ctx->h_flk, chr_name);
            if (k != kh_end(ctx->h_flk))
                flk_node = bed_depnode_list(&kh_val(ctx->h_flk, k));

            if (ctx->h_gtf) {
                k = kh_get(reg, ctx->h_gtf, chr_name);
                if (k != kh_end(ctx->h_gtf))
                    gtf_node = bed_depnode_list(&kh_val(ctx->h_gtf, k));
            }
        }

        if (!tgt_node) continue;

        /* flank */
        while (flk_node && flk_node->stop < c->pos + 1) {
            mts_stat_flk(flk_node, ctx->c_flkdep); del_node(flk_node);
        }
        if (flk_node && readcore(flk_node, b, state)) ctx->fs.n_flk++;

        /* target */
        while (tgt_node && tgt_node->stop < c->pos + 1) {
            mts_stat_region(tgt_node, ctx->c_dep, ctx->c_rmdupdep,
                            ctx->c_reg, depvals_chr);
            del_node(tgt_node);
            if (tgt_node && isZero(tgt_node->len)) depnode_init(tgt_node);
        }
        if (tgt_node && readcore(tgt_node, b, state)) ctx->fs.n_tgt++;

        /* GTF */
        while (gtf_node && gtf_node->stop < c->pos + 1) {
            mts_stat_gtf(gtf_node, ctx->c_gtf_dep, ctx->c_gtf_rmdupdep);
            del_node(gtf_node);
            if (gtf_node && isZero(gtf_node->len)) depnode_init(gtf_node);
        }
        if (gtf_node && readcore(gtf_node, b, state)) ctx->fs.n_gtf_reads++;
    }

    /* finalize last chromosome */
    while (tgt_node) {
        mts_stat_region(tgt_node, ctx->c_dep, ctx->c_rmdupdep,
                        ctx->c_reg, depvals_chr);
        del_node(tgt_node);
    }
    while (flk_node) { mts_stat_flk(flk_node, ctx->c_flkdep); del_node(flk_node); }
    while (gtf_node) { mts_stat_gtf(gtf_node, ctx->c_gtf_dep, ctx->c_gtf_rmdupdep); del_node(gtf_node); }

    bam_destroy1(b);
    bgzf_close(fp);
    return NULL;
}

/* ---- chromosome partitioning (largest-first bin packing) ---- */
static int8_t **mt_partition(const bam_header_t *hdr, regHash_t *h_tgt)
{
    int n_chr = hdr->n_targets, i, t;
    int8_t **mask = malloc(n_threads * sizeof(int8_t *));
    for (t = 0; t < n_threads; ++t) mask[t] = calloc(n_chr, sizeof(int8_t));

    uint64_t *sz = calloc(n_chr, sizeof(uint64_t));
    for (i = 0; i < n_chr; ++i) {
        khiter_t k = kh_get(reg, h_tgt, hdr->target_name[i]);
        if (k != kh_end(h_tgt)) {
            bedreglist_t *ta = &kh_val(h_tgt, k);
            int j;
            for (j = 0; j < ta->m; ++j)
                sz[i] += ((uint32_t)ta->a[j]) - (ta->a[j] >> 32);
        }
    }
    uint64_t *ld = calloc(n_threads, sizeof(uint64_t));
    for (i = 0; i < n_chr; ++i) {
        if (!sz[i]) continue;
        int best = 0;
        for (t = 1; t < n_threads; ++t) if (ld[t] < ld[best]) best = t;
        mask[best][i] = 1; ld[best] += sz[i];
    }
    free(sz); free(ld);
    return mask;
}

/*
 * load_bamfiles_mt — multi-threaded dispatch.
 * Partitions chromosomes, spawns workers, merges results.
 * Pre-pass MUST have been called first (sets mt_prepass_done).
 */
static int load_bamfiles_mt(struct opt_aux *f, aux_t *a, bamflag_t *fs,
                             const bam_index_t *bai)
{
    (void)bai; /* BAI presence validated by caller; workers scan sequentially */
    int t;
    int8_t **mask = mt_partition(a->h, a->h_tgt);
    struct mt_ctx *ctx = calloc(n_threads, sizeof(struct mt_ctx));

    for (t = 0; t < n_threads; ++t) {
        ctx[t].tid = t;
        ctx[t].bam_path = f->nfiles > 0 ? f->inputs[0] : NULL;
        ctx[t].hdr = a->h;
        ctx[t].h_tgt = a->h_tgt; ctx[t].h_flk = a->h_flk;
        ctx[t].h_gtf = a->h_gtf;
        ctx[t].chrom_mask = mask[t];
        ctx[t].n_chrom_total = a->h->n_targets;
        ctx[t].mapQ_lim = f->mapQ_lim;
        ctx[t].isize_lim = f->isize_lim;
        pthread_create(&ctx[t].handle, NULL, mt_worker, &ctx[t]);
    }
    for (t = 0; t < n_threads; ++t) {
        pthread_join(ctx[t].handle, NULL);
        if (ctx[t].c_dep)      { mt_count_merge(a->c_dep, ctx[t].c_dep);
                                 count_destroy(ctx[t].c_dep); }
        if (ctx[t].c_rmdupdep) { mt_count_merge(a->c_rmdupdep, ctx[t].c_rmdupdep);
                                 count_destroy(ctx[t].c_rmdupdep); }
        if (ctx[t].c_isize)    { mt_count_merge(a->c_isize, ctx[t].c_isize);
                                 count_destroy(ctx[t].c_isize); }
        if (ctx[t].c_flkdep)   { mt_count_merge(a->c_flkdep, ctx[t].c_flkdep);
                                 count_destroy(ctx[t].c_flkdep); }
        if (ctx[t].c_reg)      { mt_count_merge(a->c_reg, ctx[t].c_reg);
                                 count_destroy(ctx[t].c_reg); }
        mt_flag_merge(fs, &ctx[t].fs);
        if (ctx[t].c_gtf_dep && a->c_gtf_dep) {
            mt_count_merge(a->c_gtf_dep, ctx[t].c_gtf_dep);
            count_destroy(ctx[t].c_gtf_dep);
        }
        if (ctx[t].c_gtf_rmdupdep && a->c_gtf_rmdupdep) {
            mt_count_merge(a->c_gtf_rmdupdep, ctx[t].c_gtf_rmdupdep);
            count_destroy(ctx[t].c_gtf_rmdupdep);
        }
    }
    /* Fill unprocessed chromosomes with zero histograms
     * (mirrors check_reachable_regions in single-threaded path). */
    {
        int ci;
        for (ci = 0; ci < a->h->n_targets; ++ci) {
            khiter_t k = kh_get(reg, a->h_tgt, a->h->target_name[ci]);
            if (k != kh_end(a->h_tgt)) {
                bedreglist_t *t = &kh_val(a->h_tgt, k);
                if (!t->data) {
                    count32_t *z;
                    count32_init(z);
                    t->data = (void *)z;
                }
            }
        }
    }

    for (t = 0; t < n_threads; ++t) free(mask[t]);
    free(mask); free(ctx);
    return 0;
}

/* ---- end multi-threading support ---- */


/*
 * Core processing loop: iterate through all BAM reads and accumulate depth
 * statistics for target and flank regions.
 *
 * Algorithm:
 *   For each read:
 *     1. Classify via flagstat() → determine if clean/dup/qc_fail/unmapped.
 *     2. If unmapped or QC-fail → skip (goto endcore).
 *     3. If chromosome changes → finalize previous chromosome's regions,
 *        initialize depnode list for the new chromosome.
 *     4. Walk CIGAR: for each aligned base, call match_pos() which checks
 *        if the position falls in any target/flank depnode and increments
 *        the appropriate depth counters.
 *     5. If the read has passed a depnode's stop position, finalize
 *        the node (stat_each_region) and free it.
 *
 * Assumptions:
 *   - BAM files MUST be sorted by coordinate (checked at runtime).
 *   - Multiple BAM files share the same header/sort order (first file's
 *     header is used for chromosome name lookups).
 *
 * Output (streamed during processing):
 *   - depth.tsv.gz: per-base depths for every target position
 *   - region.tsv.gz: per-region summary (mean depth, median, coverage)
 *   - uncover.bed:  regions with depth below uncover_cutoff
 */
int load_bamfiles(struct opt_aux *f, aux_t *a, bamflag_t *fs)
{
    // get the chromosome name from header
    bam_header_t *h = a->h;
    loopbams_parameters_t *para = init_loopbams_parameters();
    ksprintf(para->pdepths, "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth\n");
    ksprintf(para->rcov, "#Chr\tStart\tStop\tAvg depth\tMedian\tCoverage\tCoverage(FIX)\n");
    if (outdir) chdir(outdir);
    h_uncov_init();
    int i;
    for (i = 0; i < a->ndata; ++i)
    {
        bamFile dat = a->data[i];
        bool goto_next_chromosome = FALSE;
        int ret;
        cntstat_t state;
        // main loop
        bam1_t *b;
        b = (bam1_t *)needmem(sizeof(bam1_t));

        while (1)
        {
            state = CMATCH;
            ret = bam_read1(dat, b);
            if (ret == -1)
            {
                break; // normal end
            }
            if (ret == -2)
            {
                errabort("%d bam file is truncated!\n", i + 1);
            }
            bam1_core_t *c = &b->core;
            // skip secondary and supplement alignment first
            if (c->flag & BAM_FSECONDARY || c->flag & BAM_FSUPPLEMENTARY)
                continue;

            flagstat(fs, c, ret);

            // People usually want to know how many aligned reads with mapping quality
            // greater than or equally to 10 or 20.. if (c->qual > f->mapQ_lim)
            // fs->n_qual++;
            if (c->qual >= f->mapQ_lim)
                fs->n_qual++;
            else
                state = CLOWQ;

            if (c->tid == -1 || ret == -3)
                continue; // goto endcore; // unmapped~

            if (ret == 2)
            {
                state = CDUP;
            }
            else
            {
                if (c->flag & BAM_FREAD1)
                    fs->n_rmdup1++;
                if (c->flag & BAM_FREAD2)
                    fs->n_rmdup2++;
            }

            /* stat the insertsize */
            if (c->isize > 0 && c->isize < f->isize_lim)
            {
                count_increase(a->c_isize, c->isize, uint32_t);
            }

            if (para->tid == c->tid)
            {
                if (goto_next_chromosome)
                    continue; // goto endcore;
                // Only accepted sorted bam files for effective
                if (para->lstpos > c->pos)
                {
                    errabort("The bam file is not sorted!");
                }
            }
            para->lstpos = c->pos; // 0-base

            // if tid != c->tid, completed the stat of the last chromosome
            // init the new chromosome node and clean the memory
            // FIXME: need multi thread to improve it or NOT?
            if (para->tid != c->tid)
            {
                goto_next_chromosome = FALSE; // clean the flag of skiping

                // still have reachable target node
                // release the left nodes
                if (para->tgt_node && para->tid >= 0)
                    while (para->tgt_node)
                    {
                        stat_each_region(para, a);
                        del_node(para->tgt_node);
                        // no need allocate memory space for these node, becase they are all
                        // 0s
                    }

                if (para->flk_node && para->tid >= 0)
                {
                    stat_flk_depcnt(para, a);
                    while (para->flk_node)
                    {
                        count_increaseN(a->c_flkdep, 0, para->flk_node->len, uint32_t);
                        del_node(para->flk_node);
                    }
                }

                /* GTF parallel cleanup */
                if (para->gtf_tgt_node && para->tid >= 0 && a->h_gtf)
                {
                    while (para->gtf_tgt_node)
                    {
                        stat_each_region_gtf(para, a);
                        del_node(para->gtf_tgt_node);
                    }
                }

                /* GTF stranded cleanup */
                if (gtf_stranded && para->tid >= 0) {
                    while (para->gtf_fwd.tgt_node) {
                        gtf_s_depnode_stat(para->gtf_fwd.tgt_node,
                                           &a->gtf_fwd);
                        del_node(para->gtf_fwd.tgt_node);
                    }
                    while (para->gtf_rev.tgt_node) {
                        gtf_s_depnode_stat(para->gtf_rev.tgt_node,
                                           &a->gtf_rev);
                        del_node(para->gtf_rev.tgt_node);
                    }
                }

                para->tid = c->tid;
                para->name = h->target_name[c->tid];

                // impossible in normal pratices, only happens in the different bam
                // headers
                if (para->tid > a->nchr)
                {
                    errabort("chromosome %s is not in bam header!"
                             "It must be use different bam headers",
                             para->name);
                }
                khiter_t k;
                k = kh_get(reg, a->h_tgt, para->name);
                if (k == kh_end(a->h_tgt))
                {
                    para->tgt_node = NULL;
                    para->flk_node = NULL;
                    goto_next_chromosome = TRUE;
                    continue;
                    // goto endcore;
                }
                para->tar = &kh_val(a->h_tgt, k);
                para->flk = &kh_val(a->h_flk, k);
                count32_t *tmp;
                count32_init(tmp);
                para->depvals_of_chr = tmp;
                para->tar->data = (void *)para->depvals_of_chr;
                if (para->tar->flag || para->flk->flag)
                {
                    errabort("bam files are not properly sorted\n");
                }

                para->tgt_node = bed_depnode_list(para->tar);
                para->flk_node = bed_depnode_list(para->flk);
                para->tar->flag = para->flk->flag = 1;

                /* GTF parallel: init GTF depnodes for this chromosome */
                if (a->h_gtf) {
                    k = kh_get(reg, a->h_gtf, para->name);
                    if (k != kh_end(a->h_gtf)) {
                        para->gtf_tar = &kh_val(a->h_gtf, k);
                        count32_t *gtf_tmp;
                        count32_init(gtf_tmp);
                        para->gtf_depvals_of_chr = gtf_tmp;
                        para->gtf_tar->data = (void *)para->gtf_depvals_of_chr;
                        para->gtf_tgt_node = bed_depnode_list(para->gtf_tar);
                    } else {
                        para->gtf_tar = NULL;
                        para->gtf_tgt_node = NULL;
                        para->gtf_depvals_of_chr = NULL;
                    }
                }

                /* GTF stranded: init fwd/rev depnodes */
                if (gtf_stranded) {
                    khiter_t kf, kr;
                    kf = kh_get(reg, a->gtf_fwd.h_tgt, para->name);
                    if (kf != kh_end(a->gtf_fwd.h_tgt)) {
                        para->gtf_fwd.tar = &kh_val(a->gtf_fwd.h_tgt, kf);
                        para->gtf_fwd.tgt_node = bed_depnode_list(para->gtf_fwd.tar);
                    } else {
                        para->gtf_fwd.tar = NULL;
                        para->gtf_fwd.tgt_node = NULL;
                    }
                    kr = kh_get(reg, a->gtf_rev.h_tgt, para->name);
                    if (kr != kh_end(a->gtf_rev.h_tgt)) {
                        para->gtf_rev.tar = &kh_val(a->gtf_rev.h_tgt, kr);
                        para->gtf_rev.tgt_node = bed_depnode_list(para->gtf_rev.tar);
                    } else {
                        para->gtf_rev.tar = NULL;
                        para->gtf_rev.tgt_node = NULL;
                    }
                }

                /* the next part is init uncover region hash*/
                k = kh_put(reg, h_uncov, strdup(para->name), &ret);

                bedreglist_t ucreg_tmp = {0};
                kh_val(h_uncov, k) = ucreg_tmp;
                para->ucreg = &kh_val(h_uncov, k);
                                                   /* finish init */
            }
            while (para->flk_node && para->flk_node->stop < para->lstpos + 1)
            {
                stat_flk_depcnt(para, a);
            }

            if (para->flk_node && readcore(para->flk_node, b, state))
                fs->n_flk++;

            while (para->tgt_node && para->tgt_node->stop < para->lstpos + 1)
            {
                stat_each_region(para, a);
                del_node(para->tgt_node);
                if (para->tgt_node && isZero(para->tgt_node->len))
                    depnode_init(para->tgt_node);
            }
            if (para->tgt_node && readcore(para->tgt_node, b, state))
            {
                if (export_target_bam)
                    bam_write1(bamoutfp, b);
                fs->n_tgt++;
            }

            /* GTF parallel: advance and count GTF depnodes */
            while (para->gtf_tgt_node &&
                   para->gtf_tgt_node->stop < para->lstpos + 1)
            {
                stat_each_region_gtf(para, a);
                del_node(para->gtf_tgt_node);
                if (para->gtf_tgt_node &&
                    isZero(para->gtf_tgt_node->len))
                    depnode_init(para->gtf_tgt_node);
            }
            if (para->gtf_tgt_node &&
                readcore(para->gtf_tgt_node, b, state))
            {
                fs->n_gtf_reads++;
            }

            /* GTF stranded: count by read strand */
            if (gtf_stranded) {
                int read_strand = (c->flag & BAM_FREVERSE) ? 1 : 0;
                struct depnode **gtf_node = read_strand
                    ? &para->gtf_rev.tgt_node
                    : &para->gtf_fwd.tgt_node;
                struct gtf_strand_aux *gs = read_strand
                    ? &a->gtf_rev : &a->gtf_fwd;

                while (*gtf_node &&
                       (*gtf_node)->stop < para->lstpos + 1)
                {
                    struct depnode *cur = *gtf_node;
                    gtf_s_depnode_stat(cur, gs);
                    *gtf_node = cur->next;
                    freemem(cur->vals);
                    freemem(cur->cnts);
                    freemem(cur->rmdupdep);
                    freemem(cur->covdep);
                    freemem(cur);
                    if (*gtf_node && isZero((*gtf_node)->len))
                        depnode_init(*gtf_node);
                }
                if (*gtf_node && readcore(*gtf_node, b, state)) {
                    if (read_strand)
                        fs->n_gtf_rev_reads++;
                    else
                        fs->n_gtf_fwd_reads++;
                }
            }

            // endcore:
        }
        bam_destroy1(b);
        bgzf_close(a->data[i]);
    }
    while (para->tgt_node)
    {
        stat_each_region(para, a);
        del_node(para->tgt_node); // no need allocate memory for these nodes
    }
    while (para->flk_node)
    {
        count_increaseN(a->c_flkdep, 0, para->flk_node->len, uint32_t);
        del_node(para->flk_node); // no need allocate memory for these nodes
    }
    /* GTF parallel cleanup */
    while (para->gtf_tgt_node)
    {
        stat_each_region_gtf(para, a);
        del_node(para->gtf_tgt_node);
    }
    /* GTF stranded cleanup */
    while (para->gtf_fwd.tgt_node) {
        gtf_s_depnode_stat(para->gtf_fwd.tgt_node, &a->gtf_fwd);
        del_node(para->gtf_fwd.tgt_node);
    }
    while (para->gtf_rev.tgt_node) {
        gtf_s_depnode_stat(para->gtf_rev.tgt_node, &a->gtf_rev);
        del_node(para->gtf_rev.tgt_node);
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
    uint64_t cnt, cnt4, cnt10, cnt30, cnt100, cnt02x, cnt05x;
    float cov, cov4, cov10, cov30, cov100, cov02x, cov05x;
    /* custom cutoff coverage (--cutoffdepth) */
    uint64_t cnt_array[10];
    float cov_array[10];
    /* Coverage after duplicate removal (rmdup depth).
     * These use the rmdup depth distribution (c_rmdupdep), which excludes
     * PCR duplicates, secondary alignments, and low-MAPQ reads. */
    uint64_t cnt_rmdup, cnt4_rmdup, cnt10_rmdup, cnt30_rmdup, cnt100_rmdup;
    float cov_rmdup, cov4_rmdup, cov10_rmdup, cov30_rmdup, cov100_rmdup;
    uint64_t cnt_array_rmdup[10];
    float cov_array_rmdup[10];
    /* Coverage at multiples of average depth (--depthratio).
     * e.g. ratio=0.5 means "fraction of positions with depth >= 0.5 * avg_depth". */
    uint64_t cnt_ratio[10];
    float cov_ratio[10];
};

struct regcov *regcov_init()
{
    struct regcov *c;
    c = (struct regcov *)needmem(sizeof(struct regcov));
    return c;
}

uint64_t cntcov_cal(struct opt_aux *f, struct regcov *cov, count32_t *cnt, uint64_t *data)
{
    uint64_t rawcnt = 0;
    int i;
    *data = 0;
    *cov = (struct regcov){};
    for (i = 0; i < cnt->m; ++i)
    {
        (*data) += cnt->a[i] * i;
        rawcnt += cnt->a[i];
        if (i < 100)
            cov->cnt100 += cnt->a[i];
        if (i < 30)
            cov->cnt30 += cnt->a[i];
        if (i < 10)
            cov->cnt10 += cnt->a[i];
        if (i < 4)
            cov->cnt4 += cnt->a[i];
        if (f->cutoff)
        {
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                if (i < f->cutoffs[x])
                    cov->cnt_array[x] += cnt->a[i];
            }
        }
    }
    if (rawcnt == 0)
        return 0;
    cov->cnt = rawcnt - (uint64_t)cnt->a[0];
    cov->cnt4 = rawcnt - cov->cnt4;
    cov->cnt10 = rawcnt - cov->cnt10;
    cov->cnt30 = rawcnt - cov->cnt30;
    cov->cnt100 = rawcnt - cov->cnt100;
    cov->cov = (float)cov->cnt / rawcnt * 100;
    cov->cov4 = (float)cov->cnt4 / rawcnt * 100;
    cov->cov10 = (float)cov->cnt10 / rawcnt * 100;
    cov->cov30 = (float)cov->cnt30 / rawcnt * 100;
    cov->cov100 = (float)cov->cnt100 / rawcnt * 100;
    if (f->cutoff)
    {
        for (int x = 0; x < f->num_cutoffs; x++)
        {
            cov->cnt_array[x] = rawcnt - cov->cnt_array[x];
            cov->cov_array[x] = (float)cov->cnt_array[x] / rawcnt * 100;
        }
    }
    return rawcnt;
}

uint64_t cntcov_cal2(struct opt_aux *f, struct regcov *cov, count32_t *cnt, uint64_t *data, uint64_t tgt_len)
{
    uint64_t rawcnt = 0;
    int i;
    *data = 0;
    *cov = (struct regcov){};
    // average depth
    for (i = 0; i < cnt->m; ++i)
    {
        (*data) += cnt->a[i] * i;
    }
    uint64_t avg = (*data) / tgt_len;
    uint64_t avg_02 = 0.2 * avg;
    uint64_t avg_05 = 0.5 * avg;
    for (i = 0; i < cnt->m; ++i)
    {
        // (*data) += cnt->a[i] * i;
        rawcnt += cnt->a[i];
        if (i < 100)
            cov->cnt100 += cnt->a[i];
        if (i < 30)
            cov->cnt30 += cnt->a[i];
        if (i < 10)
            cov->cnt10 += cnt->a[i];
        if (i < 4)
            cov->cnt4 += cnt->a[i];
        if (i < avg_02)
            cov->cnt02x += cnt->a[i];
        if (i < avg_05)
            cov->cnt05x += cnt->a[i];
        if (f->cutoff)
        {
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                if (i < f->cutoffs[x])
                    cov->cnt_array[x] += cnt->a[i];
            }
        }
    }
    if (rawcnt == 0)
        return 0;
    cov->cnt = rawcnt - (uint64_t)cnt->a[0];
    cov->cnt4 = rawcnt - cov->cnt4;
    cov->cnt10 = rawcnt - cov->cnt10;
    cov->cnt30 = rawcnt - cov->cnt30;
    cov->cnt100 = rawcnt - cov->cnt100;
    cov->cnt02x = rawcnt - cov->cnt02x;
    cov->cnt05x = rawcnt - cov->cnt05x;
    cov->cov = (float)cov->cnt / rawcnt * 100;
    cov->cov4 = (float)cov->cnt4 / rawcnt * 100;
    cov->cov10 = (float)cov->cnt10 / rawcnt * 100;
    cov->cov30 = (float)cov->cnt30 / rawcnt * 100;
    cov->cov100 = (float)cov->cnt100 / rawcnt * 100;
    cov->cov02x = (float)cov->cnt02x / rawcnt * 100;
    cov->cov05x = (float)cov->cnt05x / rawcnt * 100;
    if (f->cutoff)
    {
        for (int x = 0; x < f->num_cutoffs; x++)
        {
            cov->cnt_array[x] = rawcnt - cov->cnt_array[x];
            cov->cov_array[x] = (float)cov->cnt_array[x] / rawcnt * 100;
        }
    }
    return rawcnt;
}

/*
 * Compute coverage statistics using the duplicate-removed (rmdup) depth distribution.
 *
 * The rmdup depth (c_rmdupdep) counts only reads that are:
 *   - primary alignment (not secondary/supplementary)
 *   - not marked as PCR duplicate (FLAG 0x400)
 *   - MAPQ >= mapQ_lim (default 20)
 *
 * This gives a conservative estimate of unique molecule coverage,
 * excluding PCR artefacts and low-quality mappings.
 *
 * Populates the regcov->*_rmdup fields with counts and percentages
 * for thresholds: >0x, >=4x, >=10x, >=30x, >=100x, plus any --cutoffdepth values.
 */
uint64_t cntcov_cal_rmdup(struct opt_aux *f, struct regcov *cov, count32_t *cnt)
{
    uint64_t rawcnt = 0;
    uint64_t data = 0;
    int i;

    for (i = 0; i < cnt->m; ++i)
    {
        data += cnt->a[i] * i;
        rawcnt += cnt->a[i];
    }

    if (rawcnt == 0)
        return 0;

    uint64_t below_100 = 0, below_30 = 0, below_10 = 0, below_4 = 0;
    for (i = 0; i < cnt->m; ++i)
    {
        if (i < 100) below_100 += cnt->a[i];
        if (i < 30)  below_30  += cnt->a[i];
        if (i < 10)  below_10  += cnt->a[i];
        if (i < 4)   below_4   += cnt->a[i];
        if (f->cutoff)
        {
            for (int x = 0; x < f->num_cutoffs; x++)
                if (i < f->cutoffs[x])
                    cov->cnt_array_rmdup[x] += cnt->a[i];
        }
    }

    cov->cnt_rmdup  = rawcnt - (uint64_t)cnt->a[0];
    cov->cnt4_rmdup  = rawcnt - below_4;
    cov->cnt10_rmdup = rawcnt - below_10;
    cov->cnt30_rmdup = rawcnt - below_30;
    cov->cnt100_rmdup = rawcnt - below_100;

    cov->cov_rmdup  = (float)cov->cnt_rmdup  / rawcnt * 100;
    cov->cov4_rmdup  = (float)cov->cnt4_rmdup  / rawcnt * 100;
    cov->cov10_rmdup = (float)cov->cnt10_rmdup / rawcnt * 100;
    cov->cov30_rmdup = (float)cov->cnt30_rmdup / rawcnt * 100;
    cov->cov100_rmdup = (float)cov->cnt100_rmdup / rawcnt * 100;

    if (f->cutoff)
    {
        for (int x = 0; x < f->num_cutoffs; x++)
        {
            cov->cnt_array_rmdup[x] = rawcnt - cov->cnt_array_rmdup[x];
            cov->cov_array_rmdup[x] = (float)cov->cnt_array_rmdup[x] / rawcnt * 100;
        }
    }

    return rawcnt;
}

/*
 * Compute coverage at multiples of the average target depth.
 *
 * For each user-specified ratio R (e.g. 0.2, 0.5), calculates:
 *   threshold = R * avg_depth
 *   coverage  = fraction of positions with depth >= threshold
 *
 * Example: avg_depth = 100x, ratio = 0.2  → coverage at >=20x
 *          avg_depth = 100x, ratio = 0.5  → coverage at >=50x
 *
 * This is useful for assessing uniformity: a well-covered sample should have
 * high coverage at 0.2x and 0.5x of the mean depth.
 *
 * Uses a single-pass prefix sum to compute all ratios in O(n) time,
 * rather than scanning the distribution once per ratio.
 */
void cntcov_cal_ratio(struct opt_aux *f, struct regcov *cov, count32_t *cnt, uint64_t tgt_len)
{
    if (!f->depth_ratio || f->num_ratios == 0)
        return;

    /* First pass: compute total bases and total positions */
    uint64_t total_data = 0, rawcnt = 0;
    int i;

    for (i = 0; i < cnt->m; ++i)
    {
        total_data += cnt->a[i] * i;
        rawcnt += cnt->a[i];
    }

    if (rawcnt == 0 || tgt_len == 0)
        return;

    /* Mean depth over the entire target region */
    float avg_depth = (float)total_data / tgt_len;

    /*
     * Single scan with cumulative sum: as we walk depth values from 0 upward,
     * cumsum tracks positions with depth <= current i. When i reaches the
     * threshold for a ratio, we know all positions with depth below that
     * threshold have been counted, and can compute coverage = (rawcnt - cumsum).
     */
    uint64_t cumsum = 0;
    int r = 0;
    for (i = 0; i < cnt->m && r < f->num_ratios; ++i)
    {
        cumsum += cnt->a[i];
        while (r < f->num_ratios)
        {
            uint64_t threshold = (uint64_t)(f->ratios[r] * avg_depth);
            if (i < threshold)
                break;
            cov->cnt_ratio[r] = rawcnt - cumsum;
            cov->cov_ratio[r] = (float)cov->cnt_ratio[r] / rawcnt * 100;
            r++;
        }
    }
    /* Remaining ratios with threshold 0 or beyond the depth range → 0% coverage */
    for (; r < f->num_ratios; r++)
    {
        cov->cnt_ratio[r] = 0;
        cov->cov_ratio[r] = 0;
    }
}

/*
 * Compute the median from a count distribution using cumulative-sum lookup.
 *
 * The input cnt->a[i] holds the number of positions with depth == i.
 * We walk forward accumulating counts until we reach the midpoint of the
 * total, then return the depth value at that position. O(n) time, O(1) space.
 */
float median_cnt(count32_t *cnt)
{
    int i;
    uint64_t sum = 0;
    for (i = 0; i < cnt->m; ++i)
        sum += (uint64_t)cnt->a[i];
    uint64_t med = sum / 2;
    uint64_t num = 0;
    for (i = 0; i < cnt->m; ++i)
    {
        num += (uint64_t)cnt->a[i];
        if (num >= med)
            return (float)i;
    }
    return 0;
}

float average_cnt(count32_t *cnt)
{
    int i;
    uint64_t sum = 0;
    uint64_t num = 0;
    for (i = 0; i < cnt->m; ++i)
    {
        sum += (uint64_t)cnt->a[i] * i;
        num += (uint64_t)cnt->a[i];
    }
    return (float)sum / num;
}

/* ====== format-aware output helpers ====== */

/* key-value: integer (uint64_t) */
#define FMT_KV_U64(fp, key, val) do { \
    if (output_fmt == FMT_CSV) fprintf(fp, "%s,%" PRIu64 "\n", key, (uint64_t)(val)); \
    else fprintf(fp, "%60s\t%" PRIu64 "\n", key, (uint64_t)(val)); \
} while(0)

/* key-value: integer (int) */
#define FMT_KV_INT(fp, key, val) do { \
    if (output_fmt == FMT_CSV) fprintf(fp, "%s,%d\n", key, (int)(val)); \
    else fprintf(fp, "%60s\t%d\n", key, (int)(val)); \
} while(0)

/* key-value: float with custom format (NaN-safe) */
#define FMT_KV_FLOAT(fp, key, fmt, val) do { \
    if (isnan((float)(val))) { \
        if (output_fmt == FMT_CSV) fprintf(fp, "%s,N/A\n", key); \
        else fprintf(fp, "%60s\tN/A\n", key); \
    } else { \
        if (output_fmt == FMT_CSV) fprintf(fp, "%s," fmt "\n", key, val); \
        else fprintf(fp, "%60s\t" fmt "\n", key, val); \
    } \
} while(0)

/* key-value: string */
#define FMT_KV_STR(fp, key, val) do { \
    if (output_fmt == FMT_CSV) fprintf(fp, "%s,%s\n", key, val); \
    else fprintf(fp, "%60s\t%s\n", key, val); \
} while(0)

/* plot/distribution header */
static void plot_header(FILE *fp)
{
    if (output_fmt == FMT_CSV)
        fprintf(fp, "depth,count,fraction,cumulative_count,cumulative_fraction\n");
    /* TXT: no header for plot files */
}

/* insert size distribution header */
static void isize_header(FILE *fp)
{
    if (output_fmt == FMT_CSV)
        fprintf(fp, "size,count,fraction,cumulative_count,cumulative_fraction\n");
}

/* chromosome table header */
static void chrom_header(FILE *fp, struct opt_aux *f)
{
    if (output_fmt == FMT_CSV) {
        fprintf(fp, "Chromosome,DATA(%%),Avg depth,Median,Coverage%%,Cov 4x %%,Cov 10x %%,Cov 30x %%,Cov 100x %%");
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, ",Cov %dx %%", f->cutoffs[x]);
        }
        fprintf(fp, "\n");
    } else {
        fprintf(fp, "%11s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s",
                "#Chromosome", "DATA(%)", "Avg depth", "Median",
                "Coverage%", "Cov 4x %", "Cov 10x %", "Cov 30x %", "Cov 100x %");
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, "\tCov %dx %%", f->cutoffs[x]);
        }
        fprintf(fp, "\n");
    }
}

/* chromosome table row */
static void chrom_row(FILE *fp, struct opt_aux *f, const char *name,
                      float per, float avg, float med, struct regcov *c)
{
    if (output_fmt == FMT_CSV) {
        fprintf(fp, "%s,%.2f,%.2f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f",
                name, per, avg, med, c->cov, c->cov4, c->cov10, c->cov30, c->cov100);
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, ",%.2f", c->cov_array[x]);
        }
        fprintf(fp, "\n");
    } else {
        fprintf(fp, "%11s\t%8.2f\t%8.2f\t%9.1f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f",
                name, per, avg, med, c->cov, c->cov4, c->cov10, c->cov30, c->cov100);
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, "\t%.2f", c->cov_array[x]);
        }
        fprintf(fp, "\n");
    }
}

/* chromosome table row for empty data */
static void chrom_row_empty(FILE *fp, struct opt_aux *f, const char *name)
{
    if (output_fmt == FMT_CSV) {
        fprintf(fp, "%s,0.00,0.00,0.0,0.0,0.0,0.0,0.0,0.0", name);
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, ",0.0000");
        }
        fprintf(fp, "\n");
    } else {
        fprintf(fp, "%11s\t%8.2f\t%8.2f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f",
                name, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0);
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++)
                fprintf(fp, "\t%5.4f", (float)0);
        }
        fprintf(fp, "\n");
    }
}

/* distribution data row (depth or insert size) */
static void dist_row(FILE *fp, int idx, unsigned count, float frac,
                     uint64_t cumu_cnt, float cumu_frac)
{
    if (output_fmt == FMT_CSV)
        fprintf(fp, "%d,%u,%.6f,%" PRIu64 ",%.6f\n", idx, count, frac, cumu_cnt, cumu_frac);
    else
        fprintf(fp, "%d\t%u\t%f\t%" PRIu64 "\t%f\n", idx, count, frac, cumu_cnt, cumu_frac);
}

/* ====== JSON report writer ====== */
static void write_report_json(struct opt_aux *f, aux_t *a, bamflag_t *fs,
                              struct regcov *tarcov, struct regcov *flkcov,
                              struct regcov *regcov, float iavg, uint64_t imed)
{
    int i;
    FILE *fp = open_wfile("report.json");
    bool first;

    fprintf(fp, "{\n");

    /* header */
    fprintf(fp, "  \"program\": \"%s\",\n", program_name);
    fprintf(fp, "  \"version\": \"%s\",\n", Version);
    fprintf(fp, "  \"input_files\": [");
    first = TRUE;
    for (i = 0; i < f->nfiles; ++i) {
        if (!first) fprintf(fp, ", ");
        fprintf(fp, "\"%s\"", f->inputs[i]);
        first = FALSE;
    }
    fprintf(fp, "],\n");

    /* reads section */
    fprintf(fp, "  \"reads\": {\n");
    fprintf(fp, "    \"total\": %" PRIu64 ",\n", fs->n_reads);
    fprintf(fp, "    \"qc_fail\": %" PRIu64 ",\n", fs->n_qcfail);
    fprintf(fp, "    \"raw_data_mb\": %.2f,\n", (float)fs->n_data / 1e6);
    fprintf(fp, "    \"paired\": %" PRIu64 ",\n", fs->n_pair_all);
    fprintf(fp, "    \"mapped\": %" PRIu64 ",\n", fs->n_mapped);
    fprintf(fp, "    \"mapped_pct\": %.2f,\n", (float)fs->n_mapped / fs->n_reads * 100);
    fprintf(fp, "    \"mapped_data_mb\": %.2f,\n", fs->n_mdata / 1e6);
    fprintf(fp, "    \"mapped_data_pct\": %.2f,\n", (float)fs->n_mdata / fs->n_data * 100);
    fprintf(fp, "    \"properly_paired\": %" PRIu64 ",\n", fs->n_pair_good);
    fprintf(fp, "    \"properly_paired_pct\": %.2f,\n", (float)fs->n_pair_good / fs->n_reads * 100);
    fprintf(fp, "    \"read_mate_paired\": %" PRIu64 ",\n", fs->n_pair_map);
    fprintf(fp, "    \"read_mate_paired_pct\": %.2f,\n", (float)fs->n_pair_map / fs->n_reads * 100);
    fprintf(fp, "    \"singletons\": %" PRIu64 ",\n", fs->n_sgltn);
    fprintf(fp, "    \"diff_chr\": %" PRIu64 ",\n", fs->n_diffchr);
    fprintf(fp, "    \"read1\": %" PRIu64 ",\n", fs->n_read1);
    fprintf(fp, "    \"read2\": %" PRIu64 ",\n", fs->n_read2);
    fprintf(fp, "    \"read1_rmdup\": %" PRIu64 ",\n", fs->n_rmdup1);
    fprintf(fp, "    \"read2_rmdup\": %" PRIu64 ",\n", fs->n_rmdup2);
    fprintf(fp, "    \"forward_strand\": %" PRIu64 ",\n", fs->n_pstrand);
    fprintf(fp, "    \"backward_strand\": %" PRIu64 ",\n", fs->n_mstrand);
    fprintf(fp, "    \"pcr_duplicate\": %" PRIu64 ",\n", fs->n_dup);
    fprintf(fp, "    \"pcr_duplicate_pct\": %.2f,\n", (float)fs->n_dup / fs->n_mapped * 100);
    fprintf(fp, "    \"mapq_cutoff\": %d,\n", f->mapQ_lim);
    fprintf(fp, "    \"mapq_above_cutoff\": %" PRIu64 ",\n", fs->n_qual);
    fprintf(fp, "    \"mapq_above_cutoff_in_all_pct\": %.2f,\n", (float)fs->n_qual / fs->n_reads * 100);
    fprintf(fp, "    \"mapq_above_cutoff_in_mapped_pct\": %.2f\n", (float)fs->n_qual / fs->n_mapped * 100);
    fprintf(fp, "  },\n");

    /* insert size */
    fprintf(fp, "  \"insert_size\": {\n");
    if (isnan(iavg))
        fprintf(fp, "    \"average\": null,\n");
    else
        fprintf(fp, "    \"average\": %.2f,\n", iavg);
    fprintf(fp, "    \"median\": %" PRIu64 "\n", imed);
    fprintf(fp, "  },\n");

    /* target section */
    fprintf(fp, "  \"target\": {\n");
    fprintf(fp, "    \"reads\": %" PRIu64 ",\n", fs->n_tgt);
    fprintf(fp, "    \"reads_in_all_pct\": %.2f,\n", (float)fs->n_tgt / fs->n_reads * 100);
    fprintf(fp, "    \"reads_in_mapped_pct\": %.2f,\n", (float)fs->n_tgt / fs->n_mapped * 100);
    fprintf(fp, "    \"data_mb\": %.2f,\n", (float)fs->n_tdata / 1e6);
    fprintf(fp, "    \"data_rmdup_mb\": %.2f,\n", (float)fs->n_trmdat / 1e6);
    fprintf(fp, "    \"data_in_all_pct\": %.2f,\n", (float)fs->n_tdata / fs->n_data * 100);
    fprintf(fp, "    \"data_in_mapped_pct\": %.2f,\n", (float)fs->n_tdata / fs->n_mdata * 100);
    fprintf(fp, "    \"region_len\": %" PRIu64 ",\n", a->tgt_len);
    fprintf(fp, "    \"avg_depth\": %.2f,\n", (float)fs->n_tdata / a->tgt_len);
    fprintf(fp, "    \"avg_depth_rmdup\": %.2f,\n", (float)fs->n_trmdat / a->tgt_len);
    fprintf(fp, "    \"coverage\": {\n");
    fprintf(fp, "      \"gt_02x_avg\": %.2f,\n", tarcov->cov02x);
    fprintf(fp, "      \"gt_05x_avg\": %.2f,\n", tarcov->cov05x);
    fprintf(fp, "      \"gt_0x\": %.2f,\n", tarcov->cov);
    fprintf(fp, "      \"ge_4x\": %.2f,\n", tarcov->cov4);
    fprintf(fp, "      \"ge_10x\": %.2f,\n", tarcov->cov10);
    fprintf(fp, "      \"ge_30x\": %.2f,\n", tarcov->cov30);
    fprintf(fp, "      \"ge_100x\": %.2f", tarcov->cov100);
    if (f->cutoff) {
        for (int x = 0; x < f->num_cutoffs; x++)
            fprintf(fp, ",\n      \"ge_%dx\": %.2f", f->cutoffs[x], tarcov->cov_array[x]);
    }
    fprintf(fp, "\n    },\n");
    /* rmdup coverage */
    fprintf(fp, "    \"coverage_rmdup\": {\n");
    fprintf(fp, "      \"gt_0x\": %.2f,\n", tarcov->cov_rmdup);
    fprintf(fp, "      \"ge_4x\": %.2f,\n", tarcov->cov4_rmdup);
    fprintf(fp, "      \"ge_10x\": %.2f,\n", tarcov->cov10_rmdup);
    fprintf(fp, "      \"ge_30x\": %.2f,\n", tarcov->cov30_rmdup);
    fprintf(fp, "      \"ge_100x\": %.2f", tarcov->cov100_rmdup);
    if (f->cutoff) {
        for (int x = 0; x < f->num_cutoffs; x++)
            fprintf(fp, ",\n      \"ge_%dx\": %.2f", f->cutoffs[x], tarcov->cov_array_rmdup[x]);
    }
    fprintf(fp, "\n    },\n");
    /* ratio-based coverage */
    if (f->depth_ratio && f->num_ratios > 0) {
        fprintf(fp, "    \"coverage_ratio\": {\n");
        for (int r = 0; r < f->num_ratios; r++) {
            fprintf(fp, "      \"gt_%.2f_avg\": %.2f", f->ratios[r], tarcov->cov_ratio[r]);
            if (r + 1 < f->num_ratios) fprintf(fp, ",");
            fprintf(fp, "\n");
        }
        fprintf(fp, "    },\n");
    }
    fprintf(fp, "    \"region_count\": %u,\n", a->tgt_nreg);
    fprintf(fp, "    \"regions_covered_gt_0x\": %" PRIu64 ",\n", regcov->cnt);
    fprintf(fp, "    \"regions_fraction\": {\n");
    fprintf(fp, "      \"gt_0x\": %.2f,\n", regcov->cov);
    fprintf(fp, "      \"ge_4x\": %.2f,\n", regcov->cov4);
    fprintf(fp, "      \"ge_10x\": %.2f,\n", regcov->cov10);
    fprintf(fp, "      \"ge_30x\": %.2f,\n", regcov->cov30);
    fprintf(fp, "      \"ge_100x\": %.2f", regcov->cov100);
    if (f->cutoff) {
        for (int x = 0; x < f->num_cutoffs; x++)
            fprintf(fp, ",\n      \"ge_%dx\": %.2f", f->cutoffs[x], regcov->cov_array[x]);
    }
    fprintf(fp, "\n    }\n");
    fprintf(fp, "  },\n");

    /* flank section */
    fprintf(fp, "  \"flank\": {\n");
    fprintf(fp, "    \"size\": %u,\n", flank_reg);
    fprintf(fp, "    \"region_len\": %" PRIu64 ",\n", a->flk_len);
    fprintf(fp, "    \"avg_depth\": %.2f,\n", (float)fs->n_fdata / a->flk_len);
    fprintf(fp, "    \"reads\": %" PRIu64 ",\n", fs->n_flk);
    fprintf(fp, "    \"reads_in_all_pct\": %.2f,\n", (float)fs->n_flk / fs->n_reads * 100);
    fprintf(fp, "    \"reads_in_mapped_pct\": %.2f,\n", (float)fs->n_flk / fs->n_mapped * 100);
    fprintf(fp, "    \"data_mb\": %.2f,\n", (float)fs->n_fdata / 1e6);
    fprintf(fp, "    \"data_in_all_pct\": %.2f,\n", (float)fs->n_fdata / fs->n_data * 100);
    fprintf(fp, "    \"data_in_mapped_pct\": %.2f,\n", (float)fs->n_fdata / fs->n_mdata * 100);
    fprintf(fp, "    \"coverage\": {\n");
    fprintf(fp, "      \"gt_0x\": %.2f,\n", flkcov->cov);
    fprintf(fp, "      \"ge_4x\": %.2f,\n", flkcov->cov4);
    fprintf(fp, "      \"ge_10x\": %.2f,\n", flkcov->cov10);
    fprintf(fp, "      \"ge_30x\": %.2f,\n", flkcov->cov30);
    fprintf(fp, "      \"ge_100x\": %.2f", flkcov->cov100);
    if (f->cutoff) {
        for (int x = 0; x < f->num_cutoffs; x++)
            fprintf(fp, ",\n      \"ge_%dx\": %.2f", f->cutoffs[x], flkcov->cov_array[x]);
    }
    fprintf(fp, "\n    }\n");
    fprintf(fp, "  },\n");

    /* chromosomes */
    fprintf(fp, "  \"chromosomes\": [\n");
    first = TRUE;
    khiter_t k;
    for (k = 0; k != kh_end(a->h_tgt); ++k) {
        if (kh_exist(a->h_tgt, k)) {
            char *name = (char *)kh_key(a->h_tgt, k);
            count32_t *cnt = (count32_t *)kh_val(a->h_tgt, k).data;
            uint64_t data = 0;
            struct regcov *chrcov = regcov_init();
            cntcov_cal(f, chrcov, cnt, &data);
            uint64_t length = 0;
            for (i = 0; i < cnt->m; ++i) length += cnt->a[i];
            if (!first) fprintf(fp, ",\n");
            fprintf(fp, "    {\"name\": \"%s\"", name);
            if (data > 0) {
                float avg = (float)data / length;
                float med = median_cnt(cnt);
                float per = (float)data / fs->n_tdata * 100.0;
                fprintf(fp, ", \"data_pct\": %.2f, \"avg_depth\": %.2f, \"median\": %.1f"
                        ", \"coverage\": %.2f, \"cov_4x\": %.2f, \"cov_10x\": %.2f"
                        ", \"cov_30x\": %.2f, \"cov_100x\": %.2f",
                        per, avg, med, chrcov->cov, chrcov->cov4, chrcov->cov10,
                        chrcov->cov30, chrcov->cov100);
                if (f->cutoff) {
                    for (int x = 0; x < f->num_cutoffs; x++)
                        fprintf(fp, ", \"cov_%dx\": %.2f", f->cutoffs[x], chrcov->cov_array[x]);
                }
            } else {
                fprintf(fp, ", \"data_pct\": 0, \"avg_depth\": 0, \"median\": 0"
                        ", \"coverage\": 0, \"cov_4x\": 0, \"cov_10x\": 0"
                        ", \"cov_30x\": 0, \"cov_100x\": 0");
                if (f->cutoff) {
                    for (int x = 0; x < f->num_cutoffs; x++)
                        fprintf(fp, ", \"cov_%dx\": 0", f->cutoffs[x]);
                }
            }
            fprintf(fp, "}");
            first = FALSE;
            free(chrcov);
        }
    }
    fprintf(fp, "\n  ],\n");

    /* depth distribution */
    fprintf(fp, "  \"depth_distribution\": [\n");
    first = TRUE;
    uint64_t dcnt = 0, dcumu;
    for (i = 0; i < a->c_dep->m; ++i) dcnt += a->c_dep->a[i];
    dcumu = dcnt;
    for (i = 0; i < a->c_dep->m && dcnt > 0; ++i) {
        dcumu -= a->c_dep->a[i];
        if (!first) fprintf(fp, ",\n");
        fprintf(fp, "    {\"depth\": %d, \"count\": %u, \"fraction\": %.6f"
                ", \"cumulative_count\": %" PRIu64 ", \"cumulative_fraction\": %.6f}",
                i, a->c_dep->a[i], (float)a->c_dep->a[i] / dcnt,
                dcumu, (float)dcumu / dcnt);
        first = FALSE;
    }
    fprintf(fp, "\n  ],\n");

    /* insert size distribution */
    fprintf(fp, "  \"insert_size_distribution\": [\n");
    first = TRUE;
    uint64_t icnt = 0, icumu;
    for (i = 0; i < a->c_isize->m; ++i) icnt += a->c_isize->a[i];
    icumu = icnt;
    for (i = 0; i < a->c_isize->m && i < f->isize_lim && icnt > 0; ++i) {
        icumu -= a->c_isize->a[i];
        if (!first) fprintf(fp, ",\n");
        fprintf(fp, "    {\"size\": %d, \"count\": %u, \"fraction\": %.6f"
                ", \"cumulative_count\": %" PRIu64 ", \"cumulative_fraction\": %.6f}",
                i, a->c_isize->a[i], (float)a->c_isize->a[i] / icnt,
                icumu, (float)icumu / icnt);
        first = FALSE;
    }
    fprintf(fp, "\n  ]\n");

    fprintf(fp, "}\n");
    fclose(fp);
}

/* helper: pick file extension by format */
static const char *fmt_ext(const char *txt_name, const char *csv_name)
{
    if (output_fmt == FMT_CSV) return csv_name;
    return txt_name;
}

/* ====== main report function ====== */
int print_report(struct opt_aux *f, aux_t *a, bamflag_t *fs)
{
    int i;
    if (outdir) chdir(outdir);

    /* ---- compute shared stats ---- */
    /* insert size distribution */
    uint64_t icnt = 0;
    for (i = 0; i < a->c_isize->m; ++i)
        icnt += a->c_isize->a[i];

    /* depth distribution */
    uint64_t dcnt = 0;
    for (i = 0; i < a->c_dep->m; ++i)
        dcnt += a->c_dep->a[i];

    float iavg = average_cnt(a->c_isize);
    uint64_t imed = median_cnt(a->c_isize);

    struct regcov *tarcov = regcov_init();
    struct regcov *flkcov = regcov_init();
    struct regcov *regcov = regcov_init();
    cntcov_cal2(f, tarcov, a->c_dep, &fs->n_tdata, a->tgt_len);
    cntcov_cal(f, regcov, a->c_reg, &fs->n_fdata);
    cntcov_cal(f, flkcov, a->c_flkdep, &fs->n_fdata);
    for (i = 0; i < a->c_rmdupdep->m; ++i)
        fs->n_trmdat += a->c_rmdupdep->a[i] * i;

    /* Compute duplicate-removed coverage statistics (rmdup depth distribution) */
    cntcov_cal_rmdup(f, tarcov, a->c_rmdupdep);
    /* Compute ratio-based coverage (coverage at multiples of mean depth) */
    cntcov_cal_ratio(f, tarcov, a->c_dep, a->tgt_len);

    /* ---- JSON: single file, then return ---- */
    if (output_fmt == FMT_JSON) {
        write_report_json(f, a, fs, tarcov, flkcov, regcov, iavg, imed);
        mustfree(tarcov);
        mustfree(regcov);
        mustfree(flkcov);
        return 0;
    }

    /* ---- TXT / CSV: separate files ---- */

    /* insert size distribution */
    FILE *finsert = open_wfile(fmt_ext("insertsize.plot", "insertsize.csv"));
    isize_header(finsert);
    uint64_t icumu = icnt;
    for (i = 0; i < a->c_isize->m && i < f->isize_lim && icnt > 0; ++i) {
        icumu -= a->c_isize->a[i];
        dist_row(finsert, i, a->c_isize->a[i], (float)a->c_isize->a[i] / icnt,
                 icumu, (float)icumu / icnt);
    }
    fclose(finsert);

    /* depth distribution */
    FILE *fdep = open_wfile(fmt_ext("depth_distribution.plot", "depth_distribution.csv"));
    plot_header(fdep);
    uint64_t dcumu = dcnt;
    for (i = 0; i < a->c_dep->m && dcnt > 0; ++i) {
        dcumu -= a->c_dep->a[i];
        dist_row(fdep, i, a->c_dep->a[i], (float)a->c_dep->a[i] / dcnt,
                 dcumu, (float)dcumu / dcnt);
    }
    fclose(fdep);

    /* chromosomes report */
    FILE *fchrcov = open_wfile(fmt_ext("chromosomes.report", "chromosomes.csv"));
    chrom_header(fchrcov, f);
    {
        khiter_t k;
        for (k = 0; k != kh_end(a->h_tgt); ++k) {
            if (kh_exist(a->h_tgt, k)) {
                char *name = (char *)kh_key(a->h_tgt, k);
                count32_t *cnt = (count32_t *)kh_val(a->h_tgt, k).data;
                uint64_t data = 0;
                struct regcov *chrcov = regcov_init();
                cntcov_cal(f, chrcov, cnt, &data);
                uint64_t length = 0;
                int j;
                for (j = 0; j < cnt->m; ++j) length += cnt->a[j];
                if (data > 0) {
                    float avg = (float)data / length;
                    float med = median_cnt(cnt);
                    float per = (float)data / fs->n_tdata * 100.0;
                    chrom_row(fchrcov, f, name, per, avg, med, chrcov);
                } else {
                    chrom_row_empty(fchrcov, f, name);
                }
                free(chrcov);
            }
        }
    }
    fclose(fchrcov);

    /* coverage report */
    FILE *fc = open_wfile(fmt_ext("coverage.report", "coverage.csv"));
    if (output_fmt == FMT_TXT) {
        fprintf(fc, "## The file was created by %s\n", program_name);
        fprintf(fc, "## Version : %s\n", Version);
        fprintf(fc, "## Files : ");
        for (i = 0; i < f->nfiles; ++i)
            fprintf(fc, "%s ", f->inputs[i]);
        fprintf(fc, "\n");
    } else {
        fprintf(fc, "# program: %s\n", program_name);
        fprintf(fc, "# version: %s\n", Version);
        fprintf(fc, "# files:");
        for (i = 0; i < f->nfiles; ++i)
            fprintf(fc, " %s", f->inputs[i]);
        fprintf(fc, "\n");
    }

    FMT_KV_U64(fc, "[Total] Raw Reads (All reads)", fs->n_reads);
    FMT_KV_U64(fc, "[Total] QC Fail reads", fs->n_qcfail);
    FMT_KV_FLOAT(fc, "[Total] Raw Data(Mb)", "%.2f", (float)fs->n_data / 1e6);
    FMT_KV_U64(fc, "[Total] Paired Reads", fs->n_pair_all);
    FMT_KV_U64(fc, "[Total] Mapped Reads", fs->n_mapped);
    FMT_KV_FLOAT(fc, "[Total] Fraction of Mapped Reads", "%.2f%%", (float)fs->n_mapped / fs->n_reads * 100);
    FMT_KV_FLOAT(fc, "[Total] Mapped Data(Mb)", "%.2f", fs->n_mdata / 1e6);
    FMT_KV_FLOAT(fc, "[Total] Fraction of Mapped Data(Mb)", "%.2f%%", (float)fs->n_mdata / fs->n_data * 100);
    FMT_KV_U64(fc, "[Total] Properly paired", fs->n_pair_good);
    FMT_KV_FLOAT(fc, "[Total] Fraction of Properly paired", "%.2f%%", (float)fs->n_pair_good / fs->n_reads * 100);
    FMT_KV_U64(fc, "[Total] Read and mate paired", fs->n_pair_map);
    FMT_KV_FLOAT(fc, "[Total] Fraction of Read and mate paired", "%.2f%%", (float)fs->n_pair_map / fs->n_reads * 100);
    FMT_KV_U64(fc, "[Total] Singletons", fs->n_sgltn);
    FMT_KV_U64(fc, "[Total] Read and mate map to diff chr", fs->n_diffchr);
    FMT_KV_U64(fc, "[Total] Read1", fs->n_read1);
    FMT_KV_U64(fc, "[Total] Read2", fs->n_read2);
    FMT_KV_U64(fc, "[Total] Read1(rmdup)", fs->n_rmdup1);
    FMT_KV_U64(fc, "[Total] Read2(rmdup)", fs->n_rmdup2);
    FMT_KV_U64(fc, "[Total] forward strand reads", fs->n_pstrand);
    FMT_KV_U64(fc, "[Total] backward strand reads", fs->n_mstrand);
    FMT_KV_U64(fc, "[Total] PCR duplicate reads", fs->n_dup);
    FMT_KV_FLOAT(fc, "[Total] Fraction of PCR duplicate reads", "%.2f%%", (float)fs->n_dup / fs->n_mapped * 100);
    FMT_KV_INT(fc, "[Total] Map quality cutoff value", f->mapQ_lim);
    FMT_KV_U64(fc, "[Total] MapQuality above cutoff reads", fs->n_qual);
    FMT_KV_FLOAT(fc, "[Total] Fraction of MapQ reads in all reads", "%.2f%%", (float)fs->n_qual / fs->n_reads * 100);
    FMT_KV_FLOAT(fc, "[Total] Fraction of MapQ reads in mapped reads", "%.2f%%", (float)fs->n_qual / fs->n_mapped * 100);
    FMT_KV_FLOAT(fc, "[Insert size] Average", "%.2f", iavg);
    FMT_KV_U64(fc, "[Insert size] Median", imed);

    /* target */
    FMT_KV_U64(fc, "[Target] Target Reads", fs->n_tgt);
    FMT_KV_FLOAT(fc, "[Target] Fraction of Target Reads in all reads", "%.2f%%", (float)fs->n_tgt / fs->n_reads * 100);
    FMT_KV_FLOAT(fc, "[Target] Fraction of Target Reads in mapped reads", "%.2f%%", (float)fs->n_tgt / fs->n_mapped * 100);
    FMT_KV_FLOAT(fc, "[Target] Target Data(Mb)", "%.2f", (float)fs->n_tdata / 1e6);
    FMT_KV_FLOAT(fc, "[Target] Target Data Rmdup(Mb)", "%.2f", (float)fs->n_trmdat / 1e6);
    FMT_KV_FLOAT(fc, "[Target] Fraction of Target Data in all data", "%.2f%%", (float)fs->n_tdata / fs->n_data * 100);
    FMT_KV_FLOAT(fc, "[Target] Fraction of Target Data in mapped data", "%.2f%%", (float)fs->n_tdata / fs->n_mdata * 100);
    FMT_KV_U64(fc, "[Target] Len of region", a->tgt_len);
    FMT_KV_FLOAT(fc, "[Target] Average depth", "%.2f", (float)fs->n_tdata / a->tgt_len);
    FMT_KV_FLOAT(fc, "[Target] Average depth(rmdup)", "%.2f", (float)fs->n_trmdat / a->tgt_len);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>0.2*(Average depth)x)", "%.2f%%", tarcov->cov02x);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>0.5*(Average depth)x)", "%.2f%%", tarcov->cov05x);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>0x)", "%.2f%%", tarcov->cov);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>=4x)", "%.2f%%", tarcov->cov4);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>=10x)", "%.2f%%", tarcov->cov10);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>=30x)", "%.2f%%", tarcov->cov30);
    FMT_KV_FLOAT(fc, "[Target] Coverage (>=100x)", "%.2f%%", tarcov->cov100);
    if (f->cutoff) {
        char titles[40];
        for (int x = 0; x < f->num_cutoffs; x++) {
            sprintf(titles, "[Target] Coverage (>=%ux)", f->cutoffs[x]);
            FMT_KV_FLOAT(fc, titles, "%.2f%%", tarcov->cov_array[x]);
        }
    }
    /* rmdup depth coverage */
    FMT_KV_FLOAT(fc, "[Target] Coverage(rmdup) (>0x)", "%.2f%%", tarcov->cov_rmdup);
    FMT_KV_FLOAT(fc, "[Target] Coverage(rmdup) (>=4x)", "%.2f%%", tarcov->cov4_rmdup);
    FMT_KV_FLOAT(fc, "[Target] Coverage(rmdup) (>=10x)", "%.2f%%", tarcov->cov10_rmdup);
    FMT_KV_FLOAT(fc, "[Target] Coverage(rmdup) (>=30x)", "%.2f%%", tarcov->cov30_rmdup);
    FMT_KV_FLOAT(fc, "[Target] Coverage(rmdup) (>=100x)", "%.2f%%", tarcov->cov100_rmdup);
    if (f->cutoff) {
        char titles[60];
        for (int x = 0; x < f->num_cutoffs; x++) {
            sprintf(titles, "[Target] Coverage(rmdup) (>=%ux)", f->cutoffs[x]);
            FMT_KV_FLOAT(fc, titles, "%.2f%%", tarcov->cov_array_rmdup[x]);
        }
    }
    /* ratio-based coverage */
    if (f->depth_ratio && f->num_ratios > 0) {
        char titles[60];
        for (int r = 0; r < f->num_ratios; r++) {
            sprintf(titles, "[Target] Coverage (>%.2f*Avg)", f->ratios[r]);
            FMT_KV_FLOAT(fc, titles, "%.2f%%", tarcov->cov_ratio[r]);
        }
    }
    /* target regions */
    FMT_KV_INT(fc, "[Target] Target Region Count", a->tgt_nreg);
    FMT_KV_U64(fc, "[Target] Region covered > 0x", regcov->cnt);
    FMT_KV_FLOAT(fc, "[Target] Fraction Region covered > 0x", "%.2f%%", regcov->cov);
    FMT_KV_FLOAT(fc, "[Target] Fraction Region covered >= 4x", "%.2f%%", regcov->cov4);
    FMT_KV_FLOAT(fc, "[Target] Fraction Region covered >= 10x", "%.2f%%", regcov->cov10);
    FMT_KV_FLOAT(fc, "[Target] Fraction Region covered >= 30x", "%.2f%%", regcov->cov30);
    FMT_KV_FLOAT(fc, "[Target] Fraction Region covered >= 100x", "%.2f%%", regcov->cov100);
    if (f->cutoff) {
        char titles[60];
        for (int x = 0; x < f->num_cutoffs; x++) {
            sprintf(titles, "[Target] Fraction Region covered (>=%ux)", f->cutoffs[x]);
            FMT_KV_FLOAT(fc, titles, "%.2f%%", regcov->cov_array[x]);
        }
    }

    /* flank */
    FMT_KV_INT(fc, "[flank] flank size", flank_reg);
    FMT_KV_U64(fc, "[flank] Len of region (not include target region)", a->flk_len);
    FMT_KV_FLOAT(fc, "[flank] Average depth", "%.2f",
                 a->flk_len > 0 ? (float)fs->n_fdata / a->flk_len : 0.0f);
    FMT_KV_U64(fc, "[flank] flank Reads", fs->n_flk);
    FMT_KV_FLOAT(fc, "[flank] Fraction of flank Reads in all reads", "%.2f%%", (float)fs->n_flk / fs->n_reads * 100);
    FMT_KV_FLOAT(fc, "[flank] Fraction of flank Reads in mapped reads", "%.2f%%", (float)fs->n_flk / fs->n_mapped * 100);
    FMT_KV_FLOAT(fc, "[flank] flank Data(Mb)", "%.2f", (float)fs->n_fdata / 1e6);
    FMT_KV_FLOAT(fc, "[flank] Fraction of flank Data in all data", "%.2f%%", (float)fs->n_fdata / fs->n_data * 100);
    FMT_KV_FLOAT(fc, "[flank] Fraction of flank Data in mapped data", "%.2f%%", (float)fs->n_fdata / fs->n_mdata * 100);
    FMT_KV_FLOAT(fc, "[flank] Coverage (>0x)", "%.2f%%", flkcov->cov);
    FMT_KV_FLOAT(fc, "[flank] Coverage (>=4x)", "%.2f%%", flkcov->cov4);
    FMT_KV_FLOAT(fc, "[flank] Coverage (>=10x)", "%.2f%%", flkcov->cov10);
    FMT_KV_FLOAT(fc, "[flank] Coverage (>=30x)", "%.2f%%", flkcov->cov30);
    FMT_KV_FLOAT(fc, "[flank] Coverage (>=100x)", "%.2f%%", flkcov->cov100);
    if (f->cutoff) {
        char titles[40];
        for (int x = 0; x < f->num_cutoffs; x++) {
            sprintf(titles, "[flank] Coverage (>=%ux)", f->cutoffs[x]);
            FMT_KV_FLOAT(fc, titles, "%.2f%%", flkcov->cov_array[x]);
        }
    }

    /* ====== GTF section (only when -p and -g used together) ====== */
    if (gtf_parallel && a->c_gtf_dep && a->c_gtf_dep->m > 0) {
        const char *gtf_tag = gtf_level == GTF_LEVEL_EXON
                            ? "[GTF Exon]" : "[GTF Gene]";

        /* Compute GTF coverage stats */
        uint64_t gtf_tdata = 0;
        struct regcov *gtfcov = regcov_init();
        cntcov_cal2(f, gtfcov, a->c_gtf_dep, &gtf_tdata, a->gtf_len);
        fs->n_gtf_data = gtf_tdata;

        /* GTF rmdup */
        for (int i = 0; i < a->c_gtf_rmdupdep->m; ++i)
            fs->n_gtf_rmdpdat += a->c_gtf_rmdupdep->a[i] * i;
        cntcov_cal_rmdup(f, gtfcov, a->c_gtf_rmdupdep);

        /* GTF ratio coverage */
        cntcov_cal_ratio(f, gtfcov, a->c_gtf_dep, a->gtf_len);

        /* GTF region-level stats */
        uint64_t gtf_rdat = 0;
        struct regcov *gtfregcov = regcov_init();
        if (a->c_gtf_reg)
            cntcov_cal(f, gtfregcov, a->c_gtf_reg, &gtf_rdat);

        char buf[60];
        sprintf(buf, "%s Reads", gtf_tag);
        FMT_KV_U64(fc, buf, fs->n_gtf_reads);
        sprintf(buf, "%s Fraction of Reads in all reads", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%",
                     (float)fs->n_gtf_reads / fs->n_reads * 100);
        sprintf(buf, "%s Fraction of Reads in mapped reads", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%",
                     (float)fs->n_gtf_reads / fs->n_mapped * 100);
        sprintf(buf, "%s Data(Mb)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f",
                     (float)fs->n_gtf_data / 1e6);
        sprintf(buf, "%s Len of region", gtf_tag);
        FMT_KV_U64(fc, buf, a->gtf_len);
        sprintf(buf, "%s Average depth", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f",
                     a->gtf_len > 0
                     ? (float)fs->n_gtf_data / a->gtf_len : 0.0f);
        sprintf(buf, "%s Average depth(rmdup)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f",
                     a->gtf_len > 0
                     ? (float)fs->n_gtf_rmdpdat / a->gtf_len : 0.0f);
        sprintf(buf, "%s Coverage (>0.2*(Avg)x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov02x);
        sprintf(buf, "%s Coverage (>0.5*(Avg)x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov05x);
        sprintf(buf, "%s Coverage (>0x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov);
        sprintf(buf, "%s Coverage (>=4x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov4);
        sprintf(buf, "%s Coverage (>=10x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov10);
        sprintf(buf, "%s Coverage (>=30x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov30);
        sprintf(buf, "%s Coverage (>=100x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov100);
        if (f->cutoff) {
            for (int x = 0; x < f->num_cutoffs; x++) {
                sprintf(buf, "%s Coverage (>=%ux)", gtf_tag, f->cutoffs[x]);
                FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov_array[x]);
            }
        }
        /* GTF rmdup coverage */
        sprintf(buf, "%s Coverage(rmdup) (>0x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov_rmdup);
        sprintf(buf, "%s Coverage(rmdup) (>=4x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov4_rmdup);
        sprintf(buf, "%s Coverage(rmdup) (>=10x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov10_rmdup);
        sprintf(buf, "%s Coverage(rmdup) (>=30x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov30_rmdup);
        sprintf(buf, "%s Coverage(rmdup) (>=100x)", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov100_rmdup);
        /* GTF ratio coverage */
        if (f->depth_ratio && f->num_ratios > 0) {
            for (int r = 0; r < f->num_ratios; r++) {
                sprintf(buf, "%s Coverage (>%.2f*Avg)", gtf_tag, f->ratios[r]);
                FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfcov->cov_ratio[r]);
            }
        }
        /* GTF region-level */
        sprintf(buf, "%s Region Count", gtf_tag);
        FMT_KV_INT(fc, buf, a->gtf_nreg);
        sprintf(buf, "%s Region covered > 0x", gtf_tag);
        FMT_KV_U64(fc, buf, gtfregcov->cnt);
        sprintf(buf, "%s Fraction Region covered > 0x", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfregcov->cov);
        sprintf(buf, "%s Fraction Region covered >= 4x", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfregcov->cov4);
        sprintf(buf, "%s Fraction Region covered >= 10x", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfregcov->cov10);
        sprintf(buf, "%s Fraction Region covered >= 30x", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfregcov->cov30);
        sprintf(buf, "%s Fraction Region covered >= 100x", gtf_tag);
        FMT_KV_FLOAT(fc, buf, "%.2f%%", gtfregcov->cov100);

        mustfree(gtfcov);
        mustfree(gtfregcov);
    }

    /* ====== GTF Stranded sections (--gtf-stranded) ====== */
    if (gtf_stranded && a->gtf_fwd.c_dep && a->gtf_fwd.c_dep->m > 0) {

        const char *strand_tags[2] = {"[GTF Gene Forward]",
                                       "[GTF Gene Reverse]"};
        struct gtf_strand_aux *gs[2] = {&a->gtf_fwd, &a->gtf_rev};
        uint64_t *n_reads[2]  = {&fs->n_gtf_fwd_reads, &fs->n_gtf_rev_reads};

        for (int si = 0; si < 2; si++) {
            const char *tag = strand_tags[si];
            struct gtf_strand_aux *g = gs[si];

            if (!g->c_dep || g->c_dep->m == 0) continue;

            uint64_t gdata = 0;
            struct regcov *gc = regcov_init();
            cntcov_cal2(f, gc, g->c_dep, &gdata, g->len);
            *n_reads[si] = (*n_reads[si] > 0) ? *n_reads[si] : 0;

            cntcov_cal_rmdup(f, gc, g->c_rmdupdep);
            cntcov_cal_ratio(f, gc, g->c_dep, g->len);

            uint64_t grdat = 0;
            struct regcov *grc = regcov_init();
            if (g->c_reg) cntcov_cal(f, grc, g->c_reg, &grdat);

            char buf[80];
            sprintf(buf, "%s Reads", tag);
            FMT_KV_U64(fc, buf, *n_reads[si]);
            sprintf(buf, "%s Len of region", tag);
            FMT_KV_U64(fc, buf, g->len);
            sprintf(buf, "%s Average depth", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f",
                         g->len > 0 ? (float)gdata / g->len : 0.0f);
            sprintf(buf, "%s Coverage (>0x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov);
            sprintf(buf, "%s Coverage (>=4x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov4);
            sprintf(buf, "%s Coverage (>=10x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov10);
            sprintf(buf, "%s Coverage (>=30x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov30);
            sprintf(buf, "%s Coverage (>=100x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov100);
            sprintf(buf, "%s Coverage(rmdup) (>0x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov_rmdup);
            sprintf(buf, "%s Coverage(rmdup) (>=10x)", tag);
            FMT_KV_FLOAT(fc, buf, "%.2f%%", gc->cov10_rmdup);
            sprintf(buf, "%s Region Count", tag);
            FMT_KV_INT(fc, buf, g->nreg);

            mustfree(gc);
            mustfree(grc);
        }
    }

    fclose(fc);
    mustfree(tarcov);
    mustfree(regcov);
    mustfree(flkcov);
    return 0;
}

enum
{
    MAXDEPTH,
    CUTOFF,
    INSERTSIZE,
    UNCOVER,
    BAMOUT,
    FORMAT_OPT,
    DEPTHRATIO,
    GTFLEVEL,
    GTFSTRANDED,
    NTHREADS,
    HELP
};

static struct option const long_opts[] = {{"outdir", required_argument, NULL, 'o'},
                                          {"bed", required_argument, NULL, 'p'},
                                          {"gtf", required_argument, NULL, 'g'},
                                          {"gtf-level", required_argument, NULL, GTFLEVEL},
                                          {"gtf-stranded", no_argument, NULL, GTFSTRANDED},
                                          {"flank", required_argument, NULL, 'f'},
                                          {"maxdepth", required_argument, NULL, MAXDEPTH},
                                          {"cutoffdepth", required_argument, NULL, CUTOFF},
                                          {"isize", required_argument, NULL, INSERTSIZE},
                                          {"mapthres", required_argument, NULL, 'q'},
                                          {"uncover", required_argument, NULL, UNCOVER},
                                          {"bamout", required_argument, NULL, BAMOUT},
                                          {"format", required_argument, NULL, FORMAT_OPT},
                                          {"depthratio", required_argument, NULL, DEPTHRATIO},
                                          //{"rmdup", no_argument, NULL, 'd'},
                                          {"threads", required_argument, NULL, NTHREADS},
                                          {"help", no_argument, NULL, 'h'},
                                          {"version", no_argument, NULL, 'v'}};

int show_version()
{
    printf("%s\n", Version);
    return 1;
}
int bamdst(int argc, char *argv[])
{
    int n, i;
    char *probe = 0;

    // struct opt_aux opt = {.inputs = NULL, .isize_lim = 2000, .mapQ_lim = 20};
    struct opt_aux opt = init_opt_aux();
    while ((n = getopt_long(argc, argv, "o:p:g:f:q:F:h1v@:", long_opts, NULL)) >= 0)
    {
        switch (n)
        {
        // output dir, must have right to write
        case 'o':
            outdir = strdup(optarg);
            break;
        // capture region or just the region you interesting
        case 'p':
            probe = strdup(optarg);
            break;
        // GTF annotation file for gene/exon coverage
        case 'g':
            gtf_file = strdup(optarg);
            gtf_enabled = TRUE;
            break;
        case GTFLEVEL:
            if (strcmp(optarg, "exon") == 0)
                gtf_level = GTF_LEVEL_EXON;
            else if (strcmp(optarg, "gene") == 0)
                gtf_level = GTF_LEVEL_GENE;
            else {
                fprintf(stderr, "Unknown --gtf-level '%s'. Valid: gene, exon\n", optarg);
                usage(0);
            }
            break;
        case GTFSTRANDED:
            gtf_stranded = TRUE;
            break;
        // flk the region for more information, default is 200 bp
        case 'f':
            flank_reg = atoi(optarg);
            break;
        // max depth to considered in the cumulative distribution of depths
        case MAXDEPTH:
            opt.maxdepth = atoi(optarg);
            break;
        case CUTOFF:
            opt.cutoff = TRUE;
            char *token = strtok(optarg, ",");
            int tmp;
            while (token != NULL)
            {
                if (opt.num_cutoffs >= opt.max_cutoffs)
                {
                    fprintf(stderr, "Too many cutoff depths specified\n");
                    exit(EXIT_FAILURE);
                }
                tmp = atoi(token);
                // opt.cutoffs[opt.num_cutoffs++] = atoi(token);
                opt.cutoffs[opt.num_cutoffs++] = tmp;
                token = strtok(NULL, ",");
            }
            break;
        // uncover_cutoff must be greater than 0
        case UNCOVER:
            uncover_cutoff = atoi(optarg);
            assert(uncover_cutoff > 0);
            break;
        case INSERTSIZE:
            opt.isize_lim = atoi(optarg);
            break;
        case BAMOUT:
            export_target_bam = strdup(optarg);
            break;
        case 'F':
        case FORMAT_OPT:
            if (strcmp(optarg, "csv") == 0)
                output_fmt = FMT_CSV;
            else if (strcmp(optarg, "json") == 0)
                output_fmt = FMT_JSON;
            else if (strcmp(optarg, "txt") == 0)
                output_fmt = FMT_TXT;
            else {
                fprintf(stderr, "Unknown format '%s'. Valid: txt, csv, json\n", optarg);
                usage(0);
            }
            break;
        case DEPTHRATIO:
            {
                opt.depth_ratio = TRUE;
                char *ratio_token = strtok(optarg, ",");
                opt.num_ratios = 0;
                while (ratio_token != NULL)
                {
                    if (opt.num_ratios >= opt.max_ratios)
                    {
                        fprintf(stderr, "Too many depth ratios specified (max %d)\n", opt.max_ratios);
                        exit(EXIT_FAILURE);
                    }
                    float ratio = atof(ratio_token);
                    if (ratio <= 0.0f || ratio > 1.0f)
                    {
                        fprintf(stderr, "Invalid depth ratio: %s (must be 0 < x <= 1)\n", ratio_token);
                        exit(EXIT_FAILURE);
                    }
                    opt.ratios[opt.num_ratios++] = ratio;
                    ratio_token = strtok(NULL, ",");
                }
            }
            break;
        case 'q':
            opt.mapQ_lim = atoi(optarg);
            break;
        case 'h':
            usage(1);
            break;
        case NTHREADS:
        case '@':
            n_threads = atoi(optarg);
            if (n_threads < 1) n_threads = 1;
            break;
        case 'v':
            return show_version();
        case '1':
            zero_based = FALSE;
            break;
        // case 'd': rmdup_mark = TRUE; break;
        default:
            usage(0);
            // more help
        }
    }
    if (isNull(outdir))
        usage(0);
    if (isNull(probe) && !gtf_enabled)
        usage(0);
    if (export_target_bam && check_filename_isbam(export_target_bam))
    {
        fprintf(stderr, "--bamout must be a bam file: %s", export_target_bam);
        goto freeall;
    }

    n = argc - optind;
    /* Verify output directory exists and is writable.
     * We do not create it — the user is responsible for that. */
    {
        struct stat st;
        if (stat(outdir, &st) != 0) {
            fprintf(stderr, "Output directory does not exist: %s\n", outdir);
            exit(EXIT_FAILURE);
        }
        if (!S_ISDIR(st.st_mode)) {
            fprintf(stderr, "Output path is not a directory: %s\n", outdir);
            exit(EXIT_FAILURE);
        }
        if (access(outdir, W_OK) != 0) {
            fprintf(stderr, "Output directory is not writable: %s\n", outdir);
            exit(EXIT_FAILURE);
        }
    }
    // capable of deals with severl bam files
    aux_t *aux;
    aux = aux_init();
    if (isZero(n))
    {
        aux->data = (bamFile *)needmem(sizeof(bamFile));
        aux->data[0] = bgzf_dopen(fileno(stdin), "r");
        aux->h = bam_header_read(aux->data[0]);
        aux->ndata = 1;
        opt.nfiles = 0;
    }
    else
    {
        aux->data = (bamFile *)needmem(n * sizeof(bamFile));
        opt.nfiles = n;
        opt.inputs = (char **)needmem(n * sizeof(char *));
        for (i = 0; i < n; ++i)
        {
            bam_header_t *h_tmp;
            // h_tmp = calloc(1, sizeof(bam_header_t));
            if (STREQ(argv[optind + i], "-"))
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
            if (i == 0)
                aux->h = h_tmp;
            else
                bam_header_destroy(h_tmp);
            opt.inputs[i] = realpath(argv[optind + i], NULL);
            if (!opt.inputs[i])
                opt.inputs[i] = strdup(argv[optind + i]);
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
    bool has_bed = (probe != NULL);
    if (has_bed)
        load_bed_init(probe, aux);
    chrhash_destroy();
    if (probe) { freemem(probe); probe = NULL; }

    /* Load GTF if enabled.  When -p is absent, GTF becomes the primary
     * target (stored in h_tgt/h_flk).  When both -p and -g are present,
     * GTF is stored as a parallel set in h_gtf/h_gtf_flk. */
    if (gtf_enabled) {
        load_gtf_init(gtf_file, aux);
        if (has_bed) {
            gtf_parallel = TRUE;  /* -p + -g together */
        } else {
            /* GTF-only mode: use GTF as primary target */
            bedHand->destroy((void *)aux->h_tgt, destroy_void);
            bedHand->destroy((void *)aux->h_flk, destroy_void);
            aux->h_tgt = aux->h_gtf;
            aux->h_flk = aux->h_gtf_flk;
            aux->tgt_len = aux->gtf_len;
            aux->tgt_nreg = aux->gtf_nreg;
            /* move GTF depth histograms to primary slots */
            if (aux->c_gtf_dep) {
                count_destroy(aux->c_dep);
                aux->c_dep = aux->c_gtf_dep;
                aux->c_gtf_dep = NULL;
            }
            if (aux->c_gtf_rmdupdep) {
                count_destroy(aux->c_rmdupdep);
                aux->c_rmdupdep = aux->c_gtf_rmdupdep;
                aux->c_gtf_rmdupdep = NULL;
            }
            if (aux->c_gtf_reg) {
                count_destroy(aux->c_reg);
                aux->c_reg = aux->c_gtf_reg;
                aux->c_gtf_reg = NULL;
            }
            aux->h_gtf = NULL;
            aux->h_gtf_flk = NULL;
        }
    }
    if (aux->c_isize->n < opt.isize_lim)
    {
        unsigned *new_a = realloc(aux->c_isize->a, opt.isize_lim * sizeof(unsigned));
        if (isNull(new_a))
            errabort("Failed to reallocate isize counter");
        aux->c_isize->a = new_a;
        aux->c_isize->n = opt.isize_lim;
    }
    aux->c_isize->m = 0;
    memset(aux->c_isize->a, 0, aux->c_isize->n * sizeof(unsigned));
    aux->nchr = aux->h->n_targets;
    struct bamflag fs = {};

    /* multi-threaded: check BAI index */
    if (n_threads > 1) {
        const bam_index_t *bai = NULL;
        const char *bam_path = (opt.nfiles > 0) ? opt.inputs[0] : NULL;
        if (bam_path) {
            bai = bam_index_load(bam_path);
            if (bai == NULL) {
                int64_t *blk = NULL;
                int64_t nb = bgzf_scan_blocks(bam_path, &blk);
                if (nb > 0) {
                    fprintf(stderr, "[bamdst] No BAI; %" PRId64
                            " blocks, %d threads.\n", nb, n_threads);
                    free(blk);
                } else {
                    errabort("'%s' is not a valid BGZF/BAM file.", bam_path);
                }
            }
            /* Light pre-pass: count total reads (no CIGAR walk, no depth). */
            {
                bamFile pp = bam_open(bam_path, "r");
                if (!pp) errabort("pre-pass: cannot re-open %s", bam_path);
                bam_header_t *pp_hdr = bam_header_read(pp); /* skip header */
                bam_header_destroy(pp_hdr);
                bam1_t *bb = (bam1_t *)needmem(sizeof(bam1_t));
                int pp_ret;
                while ((pp_ret = bam_read1(pp, bb)) >= 0) {
                    bam1_core_t *cc = &bb->core;
                    if (cc->flag & BAM_FSECONDARY ||
                        cc->flag & BAM_FSUPPLEMENTARY)
                        continue;
                    int dummy;
                    flagstat(&fs, cc, dummy);
                    if (cc->qual >= opt.mapQ_lim) fs.n_qual++;
                }
                bam_destroy1(bb);
                bgzf_close(pp);
                mt_prepass_done = TRUE;
            }
            load_bamfiles_mt(&opt, aux, &fs, bai);
            if (bai) bam_index_destroy((bam_index_t *)bai);
        } else {
            warnings("stdin input: cannot seek; falling back to single-thread");
            load_bamfiles(&opt, aux, &fs);
        }
    } else {
        load_bamfiles(&opt, aux, &fs);
    }
    print_report(&opt, aux, &fs);
    aux_destroy(aux);
    for (i = 0; i < opt.nfiles; ++i)
        freemem(opt.inputs[i]);
    freemem(opt.inputs);
    if (export_target_bam)
        bam_close(bamoutfp);
    freemem(opt.cutoffs);
    freemem(opt.ratios);
freeall:
    freemem(export_target_bam);
    if (gtf_file) freemem(gtf_file);
    freemem(outdir);
    return 0;
}

/* main */
int main(int argc, char *argv[])
{
    return bamdst(argc, argv);
}
