/* The MIT License

   Copyright (c) 2026 pzweuj - Ported to HTSlib
   Copyright (c) 2022, 2023, 2024 Authors
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


#include "bedutil.h"
#include "commons.h"
#include "count.h"

// bgzf for writing tabix-able depth.gz file
#include <htslib/bgzf.h>

// htslib for BAM/CRAM/SAM file reading
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

// khash, kstring and knetfile are standard utils of klib
#include "khash.h"
#include "knetfile.h"
#include "kstring.h"
#include <sys/stat.h>

static char const *program_name = "xamdst";
static char const *Version = "2.0.0";

/* flank region will be stat in the coverage report file,
 * this value can be set by -f / --flank */
static int flank_reg = 200;

/* extern bedHand from bedutil.c, it is a collection of functions*/
extern bedHandle_t *bedHand;

/* only accepted one stdin pipeline */
static bool stdin_lock = FALSE;

/* the bed file is zero based */
static bool zero_based = TRUE;

/* 定义 max cutoff 最大值为10 */
static const int MAX_CUTOFFS = 10;

/* I/O 缓冲区大小优化 - 64KB 缓冲区减少系统调用次数 */
static const int WRITE_BUFFER_SIZE = 65536;

/* The number of threads after which there are
   diminishing performance gains. */
// enum { DEFAULT_MAX_THREADS = 8 };

/* duplicate will be removed if rmdup_mark is TRUE
 * this option is removed since version 1.0.0
 */
// static bool rmdup_mark = FALSE;

/* export target reads to a specitified bam file */
static char *export_target_bam = NULL;
htsFile *bamoutfp = NULL;

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

// warp sam header to retrieve chromosome length hash
void header2chrhash(sam_hdr_t *h)
{
    int i, ret;
    khiter_t k;
    int n_targets = sam_hdr_nref(h);
    for (i = 0; i < n_targets; i++)
    {
        const char *name = sam_hdr_tid2name(h, i);
        hts_pos_t len = sam_hdr_tid2len(h, i);
        k = kh_put(chr, h_chrlen, strdup(name), &ret);
        kh_val(h_chrlen, k).length = (uint32_t)len;
    }
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
    int max_cutoffs; // 最大允许的 cutoff 数量
    int *cutoffs;    // 指向 cutoff 数组的指针
    int mapQ_lim;
    int maxdepth;
    // 新增：深度比例参数
    bool depth_ratio;   // 是否启用深度比例统计
    int num_ratios;     // 比例数量
    int max_ratios;     // 最大比例数量
    float *ratios;      // 比例数组 (如 0.1, 0.2, 0.5)
    // 新增：参考基因组路径（用于 CRAM 文件）
    char *reference;    // 参考基因组 FASTA 文件路径
    // 新增：线程数
    int nthreads;       // htslib 多线程压缩/解压缩线程数
};

// 一个函数来初始化 opt_aux 结构体
struct opt_aux init_opt_aux()
{
    struct opt_aux opt;
    memset(&opt, 0, sizeof(opt)); // 初始化所有成员为0
    opt.num_cutoffs = 0;
    opt.max_cutoffs = MAX_CUTOFFS; // 设置最大 cutoff 数量
    opt.cutoff = FALSE;
    opt.cutoffs = malloc(opt.max_cutoffs * sizeof(int)); // 分配内存
    if (opt.cutoffs == NULL)
    {
        // 处理内存分配失败
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    // 其他成员的初始化...
    opt.inputs = NULL;
    opt.isize_lim = 2000;
    opt.mapQ_lim = 20;
    
    // 初始化深度比例参数
    opt.depth_ratio = TRUE;  // 默认启用
    opt.num_ratios = 2;      // 默认2个比例
    opt.max_ratios = MAX_CUTOFFS;
    opt.ratios = malloc(opt.max_ratios * sizeof(float));
    if (opt.ratios == NULL)
    {
        fprintf(stderr, "Memory allocation failed for ratios\n");
        exit(EXIT_FAILURE);
    }
    // 设置默认比例值 0.2 和 0.5
    opt.ratios[0] = 0.2f;
    opt.ratios[1] = 0.5f;
    
    // 初始化参考基因组路径
    opt.reference = NULL;

    // 初始化线程数（0 = 单线程模式，不使用多线程）
    opt.nthreads = 0;

    return opt;
}

struct depnode
{
    unsigned len; // len will be set to 0, if not allocated memory
    unsigned start;
    unsigned stop;
    unsigned *vals;     // raw depth
    unsigned *cnts;     // reads count in this position
    unsigned *rmdupdep; // clean depth, rmdup, mapQ > 20, primary hit
    unsigned *covdep;   // coverage depth
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

    /* data, array of htsFile pointers for BAM/CRAM/SAM
     * h,  point to sam header */
    htsFile **data;
    sam_hdr_t *h;

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

struct _aux *aux_init()
{
    struct _aux *a;
    a = calloc(1, sizeof(struct _aux));
    a->nchr = a->maxdep = a->ndata = 0;
    a->data = NULL; // htsFile array
    a->h = NULL;    // sam header
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
    count32_t *cnt = (count32_t *)data;
    count_destroy(cnt);
}

void aux_destroy(struct _aux *a)
{
    free(a->data);
    bedHand->destroy((void *)a->h_tgt, destroy_data);
    bedHand->destroy((void *)a->h_flk, destroy_void);
    sam_hdr_destroy(a->h);
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
   or : %s [OPTION] -p <probe.bed> -o <output_dir> [in1.cram] -T <ref.fa>\n\
   or : %s [OPTION] -p <probe.bed> -o <output_dir> -\n\
",
               Version, program_name, program_name, program_name);
        puts("\
Option -o and -p are mandatory:\n\
  -o, --outdir         output dir (will be created if not exists)\n\
  -p, --bed            probe or target regions file, the region file will \n\
                       be merged before calculate depths\n\
");
        puts("\
Optional parameters:\n\
   -T, --reference FILE  reference genome FASTA file (required for CRAM input)\n\
   -f, --flank [200]   flank n bp of each region\n\
   -q [20]             map quality cutoff value, greater or equal to the value will be count\n\
   --maxdepth [0]      set the max depth to stat the cumu distribution.\n\
   --cutoffdepth [0,0] list the coverage of above these depths, allow maximal 10 cutoffs.\n\
   --depthratio [0.2,0.5] coverage at ratios of average depth (e.g., 0.1,0.2,0.5)\n\
   --isize [2000]      stat the inferred insert size under this value\n\
   --uncover [5]       region will included in uncover file if below it\n\
   --bamout  BAMFILE   target reads will be exported to this bam file\n\
   --threads [0]       number of threads for BAM/CRAM I/O (0 = single-threaded)\n\
   -1                  begin position of bed file is 1-based\n\
   -h, --help          print this help info\n\
\n");

        puts("\
* Supported input formats: BAM, CRAM, SAM\n\
* CRAM files require a reference genome (-T/--reference option)\n\
* The output directory will be created automatically if it doesn't exist\n\
\n\
* Five essential files would be created in the output dir. \n\
* region.tsv.gz and depth.tsv.gz are zipped by bgzip, so you can use tabix \n\
  index these files.\n\n\
 - coverage.report     a report of the coverage information and reads \n\
                       information of whole target regions\n\
 - coverage.report.json  JSON format of coverage report for easy parsing\n\
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
      https://github.com/pzweuj/xamdst\n");
    }
    exit(EXIT_SUCCESS);
}

/**
 * @brief 如果文件夹路径不存在则创建该路径
 * @param path
 * @param mode
 * @return
 */
int mkdirp(const char *path, mode_t mode)
{
    char tmp[256];
    char *p;

    // 确保路径不以 '/' 结束，除非它是根目录
    snprintf(tmp, sizeof(tmp), "%.*s", (int)(strlen(path) < sizeof(tmp) - 1) ? strlen(path) : sizeof(tmp) - 1, path);

    for (p = &tmp[1]; *p; p++)
    {
        if (*p == '/')
        {
            *p = '\0'; // 分割点
            if (mkdir(tmp, mode) != 0 && errno != EEXIST)
            {
                fprintf(stderr, "Failed to create directory '%s': %s\n", tmp, strerror(errno));
                return -1;
            }
            *p = '/'; // 恢复路径
        }
    }
    // 创建最后的目录
    if (mkdir(tmp, mode) != 0 && errno != EEXIST)
    {
        fprintf(stderr, "Failed to create directory '%s': %s\n", tmp, strerror(errno));
        return -1;
    }
    return 0;
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
    float med = l & 1 ? tmp[(l >> 1) + 1] : (float)(tmp[l >> 1] + tmp[(l >> 1) - 1]) / 2;
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

// 加载并初始化 BED 文件，包含截断检测
int load_bed_init(char const *fn, aux_t *a)
{
    int ret = 0;
    bedHand->read(fn, a->h_tgt, 0, 0, &ret);
    
    // 检查 BED 文件是否为空
    if (kh_size(a->h_tgt) == 0)
    {
        errabort("Failed to read BED file: %s (file may be empty or corrupted)", fn);
    }
    
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
    uint32_t *cigar = bam_get_cigar(b);
    uint32_t end = bam_endpos(b);

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

loopbams_parameters_t *init_loopbams_parameters()
{
    loopbams_parameters_t *para;
    para = (loopbams_parameters_t *)needmem(sizeof(loopbams_parameters_t));
    *para = (loopbams_parameters_t){.tid = -1};
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
        // 优化：只有当缓冲区超过阈值时才写入
        if (para->pdepths->l > WRITE_BUFFER_SIZE)
        {
            write_buffer_bgzf(para->pdepths, para->fdep);
        }
    }
    else
    {
        avg = med = cov1 = cov2 = 0.0;
        for (j = 0; j < node->len; ++j)
        {
            ksprintf(para->pdepths, "%s\t%d\t0\t0\t0\n", para->name, node->start + j);
        }
        // 优化：批量写入而不是每行都写
        if (para->pdepths->l > WRITE_BUFFER_SIZE)
        {
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

    // 优化：只有当缓冲区超过阈值时才写入
    if (para->rcov->l > WRITE_BUFFER_SIZE)
    {
        write_buffer_bgzf(para->rcov, para->freg);
    }
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
            bedreglist_t *ucreg_tmp;
            ucreg_tmp = (bedreglist_t *)needmem(sizeof(bedreglist_t));
            kh_val(h_uncov, l) = *ucreg_tmp;
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
    bedHand->destroy(h_uncov, destroy_data);
}

// load bam files and stat the depths
// huge function, need IMPROVE it!!
int load_bamfiles(struct opt_aux *f, aux_t *a, bamflag_t *fs)
{
    // get the chromosome name from header
    sam_hdr_t *h = a->h;
    loopbams_parameters_t *para = init_loopbams_parameters();
    ksprintf(para->pdepths, "#Chr\tPos\tRaw Depth\tRmdup depth\tCover depth\n");
    ksprintf(para->rcov, "#Chr\tStart\tStop\tAvg depth\tMedian\tCoverage\tCoverage(FIX)\n");
    if (outdir) chdir(outdir);
    h_uncov_init();
    int i;
    for (i = 0; i < a->ndata; ++i)
    {
        htsFile *dat = a->data[i];
        bool goto_next_chromosome = FALSE;
        int ret;
        cntstat_t state;
        // main loop
        bam1_t *b;
        b = bam_init1();

        while (1)
        {
            state = CMATCH;
            ret = sam_read1(dat, h, b);
            if (ret == -1)
            {
                break; // normal end
            }
            if (ret < -1)
            {
                errabort("%d bam/cram file is truncated or corrupted!\n", i + 1);
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
            // 
            // 多线程优化说明：
            // 当前实现是单线程顺序处理，可以考虑以下多线程方案：
            // 1. 按染色体并行：每个线程处理一个染色体的深度统计
            //    - 优点：实现简单，线程间无数据竞争
            //    - 缺点：需要预先读取整个 BAM 文件或使用索引
            // 2. 生产者-消费者模式：一个线程读取 BAM，多个线程统计
            //    - 优点：可以流式处理
            //    - 缺点：需要线程安全的队列和同步机制
            // 目前保持单线程实现，因为 I/O 通常是瓶颈
            //
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

                para->tid = c->tid;
                para->name = (char *)sam_hdr_tid2name(h, c->tid);

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

                /* 初始化 uncover region hash
                 * 注意：这里不能直接使用 ucreg_tmp 指针，因为 kh_val 返回的是值的拷贝
                 * 我们需要获取 hash 表中实际存储位置的地址
                 * 这样后续对 para->ucreg 的修改才能反映到 hash 表中
                 */
                k = kh_put(reg, h_uncov, strdup(para->name), &ret);

                bedreglist_t *ucreg_tmp;
                ucreg_tmp = (bedreglist_t *)needmem(sizeof(bedreglist_t));
                memset(ucreg_tmp, 0, sizeof(bedreglist_t));  // 确保初始化为零
                kh_val(h_uncov, k) = *ucreg_tmp;
                mustfree(ucreg_tmp);  // 释放临时变量，数据已拷贝到 hash 表
                para->ucreg = &kh_val(h_uncov, k);  // 获取 hash 表中的实际地址
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
                    sam_write1(bamoutfp, a->h, b);
                fs->n_tgt++;
            }

            // endcore:
        }
        bam_destroy1(b);
        hts_close(a->data[i]);
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
    // MAX_CUTOFFS
    uint64_t cnt_array[10];
    float cov_array[10];
    
    // 新增：rmdup depth 覆盖度统计
    uint64_t cnt_rmdup, cnt4_rmdup, cnt10_rmdup, cnt30_rmdup, cnt100_rmdup;
    float cov_rmdup, cov4_rmdup, cov10_rmdup, cov30_rmdup, cov100_rmdup;
    uint64_t cnt_array_rmdup[10];
    float cov_array_rmdup[10];
    
    // 新增：基于比例的覆盖度统计
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

// 计算 rmdup depth 的覆盖度统计
uint64_t cntcov_cal_rmdup(struct opt_aux *f, struct regcov *cov, count32_t *cnt, uint64_t *data, uint64_t tgt_len)
{
    uint64_t rawcnt = 0;
    int i;
    *data = 0;
    
    // 计算总数据量和平均深度
    for (i = 0; i < cnt->m; ++i)
    {
        (*data) += cnt->a[i] * i;
        rawcnt += cnt->a[i];
    }
    
    if (rawcnt == 0)
        return 0;
    
    // 计算各阈值下的覆盖度
    uint64_t below_100 = 0, below_30 = 0, below_10 = 0, below_4 = 0;
    
    for (i = 0; i < cnt->m; ++i)
    {
        if (i < 100)
            below_100 += cnt->a[i];
        if (i < 30)
            below_30 += cnt->a[i];
        if (i < 10)
            below_10 += cnt->a[i];
        if (i < 4)
            below_4 += cnt->a[i];
        if (f->cutoff)
        {
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                if (i < f->cutoffs[x])
                    cov->cnt_array_rmdup[x] += cnt->a[i];
            }
        }
    }
    
    // 计算 rmdup 覆盖度
    cov->cnt_rmdup = rawcnt - (uint64_t)cnt->a[0];
    cov->cnt4_rmdup = rawcnt - below_4;
    cov->cnt10_rmdup = rawcnt - below_10;
    cov->cnt30_rmdup = rawcnt - below_30;
    cov->cnt100_rmdup = rawcnt - below_100;
    
    cov->cov_rmdup = (float)cov->cnt_rmdup / rawcnt * 100;
    cov->cov4_rmdup = (float)cov->cnt4_rmdup / rawcnt * 100;
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

// 计算基于比例的覆盖度统计
void cntcov_cal_ratio(struct opt_aux *f, struct regcov *cov, count32_t *cnt, uint64_t tgt_len)
{
    if (!f->depth_ratio || f->num_ratios == 0)
        return;
    
    // 计算总数据量
    uint64_t total_data = 0;
    uint64_t rawcnt = 0;
    int i;
    
    for (i = 0; i < cnt->m; ++i)
    {
        total_data += cnt->a[i] * i;
        rawcnt += cnt->a[i];
    }
    
    if (rawcnt == 0 || tgt_len == 0)
        return;
    
    // 计算平均深度
    float avg_depth = (float)total_data / tgt_len;
    
    // 对每个比例计算覆盖度
    for (int r = 0; r < f->num_ratios; r++)
    {
        uint64_t threshold = (uint64_t)(f->ratios[r] * avg_depth);
        uint64_t below_threshold = 0;
        
        for (i = 0; i < cnt->m && i <= threshold; ++i)
        {
            below_threshold += cnt->a[i];
        }
        
        cov->cnt_ratio[r] = rawcnt - below_threshold;
        cov->cov_ratio[r] = (float)cov->cnt_ratio[r] / rawcnt * 100;
    }
}

// 计算深度分布的中位数
// 使用累积分布查找，时间复杂度 O(n)
float median_cnt(count32_t *cnt)
{
    if (cnt->m == 0)
        return 0;
    
    int i;
    uint64_t sum = 0;
    
    // 计算总数
    for (i = 0; i < cnt->m; ++i)
        sum += (uint64_t)cnt->a[i];
    
    if (sum == 0)
        return 0;
    
    uint64_t med = sum / 2;
    uint64_t num = 0;
    
    // 查找中位数位置
    for (i = 0; i < cnt->m; ++i)
    {
        num += (uint64_t)cnt->a[i];
        if (num >= med)
        {
            // 对于偶数个元素，返回两个中间值的平均
            if (sum % 2 == 0 && num == med && i + 1 < cnt->m)
            {
                // 找下一个非零位置
                int j = i + 1;
                while (j < cnt->m && cnt->a[j] == 0)
                    j++;
                if (j < cnt->m)
                    return (float)(i + j) / 2.0f;
            }
            return (float)i;
        }
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

// JSON 格式化输出 - 缩进级别跟踪
static int json_indent_level = 0;
static const char *JSON_INDENT = "  ";  // 2空格缩进

static void json_newline_indent(kstring_t *s)
{
    kputc('\n', s);
    for (int i = 0; i < json_indent_level; i++)
        kputs(JSON_INDENT, s);
}

static void json_start_object(kstring_t *s)
{
    kputc('{', s);
    json_indent_level++;
}

static void json_end_object(kstring_t *s)
{
    // Remove trailing comma if present
    if (s->l > 0 && s->s[s->l - 1] == ',')
        s->l--;
    json_indent_level--;
    json_newline_indent(s);
    kputc('}', s);
}

static void json_start_array(kstring_t *s)
{
    kputc('[', s);
}

static void json_end_array(kstring_t *s)
{
    if (s->l > 0 && s->s[s->l - 1] == ',')
        s->l--;
    kputc(']', s);
}

static void json_key(kstring_t *s, const char *key)
{
    json_newline_indent(s);
    ksprintf(s, "\"%s\": ", key);
}

static void json_string(kstring_t *s, const char *key, const char *value)
{
    json_newline_indent(s);
    ksprintf(s, "\"%s\": \"%s\",", key, value);
}

static void json_int(kstring_t *s, const char *key, int64_t value)
{
    json_newline_indent(s);
    ksprintf(s, "\"%s\": %" PRId64 ",", key, value);
}

static void json_uint(kstring_t *s, const char *key, uint64_t value)
{
    json_newline_indent(s);
    ksprintf(s, "\"%s\": %" PRIu64 ",", key, value);
}

static void json_float(kstring_t *s, const char *key, float value)
{
    json_newline_indent(s);
    ksprintf(s, "\"%s\": %.2f,", key, value);
}

/* Generate JSON format coverage report */
int print_report_json(struct opt_aux *f, aux_t *a, bamflag_t *fs,
                      struct regcov *tarcov, struct regcov *flkcov, 
                      struct regcov *regcov, float iavg, uint64_t imed)
{
    if (outdir) chdir(outdir);
    FILE *fjson = open_wfile("coverage.report.json");
    if (!fjson) {
        warnings("Failed to create coverage.report.json");
        return -1;
    }
    // 设置更大的文件缓冲区
    setvbuf(fjson, NULL, _IOFBF, WRITE_BUFFER_SIZE);
    
    // 重置缩进级别
    json_indent_level = 0;
    
    kstring_t json = {0, 0, NULL};
    
    json_start_object(&json);
    
    // Version and files
    json_string(&json, "version", Version);
    json_key(&json, "files");
    json_start_array(&json);
    for (int i = 0; i < f->nfiles; ++i) {
        ksprintf(&json, "\"%s\",", f->inputs[i]);
    }
    json_end_array(&json);
    kputc(',', &json);
    
    // Total section
    json_key(&json, "total");
    json_start_object(&json);
    json_uint(&json, "raw_reads", fs->n_reads);
    json_uint(&json, "qc_fail_reads", fs->n_qcfail);
    json_float(&json, "raw_data_mb", (float)fs->n_data / 1e6);
    json_uint(&json, "paired_reads", fs->n_pair_all);
    json_uint(&json, "mapped_reads", fs->n_mapped);
    json_float(&json, "mapped_reads_fraction", (float)fs->n_mapped / fs->n_reads * 100);
    json_float(&json, "mapped_data_mb", (float)fs->n_mdata / 1e6);
    json_float(&json, "mapped_data_fraction", (float)fs->n_mdata / fs->n_data * 100);
    json_uint(&json, "properly_paired", fs->n_pair_good);
    json_float(&json, "properly_paired_fraction", (float)fs->n_pair_good / fs->n_reads * 100);
    json_uint(&json, "read_mate_paired", fs->n_pair_map);
    json_float(&json, "read_mate_paired_fraction", (float)fs->n_pair_map / fs->n_reads * 100);
    json_uint(&json, "singletons", fs->n_sgltn);
    json_uint(&json, "diff_chr", fs->n_diffchr);
    json_uint(&json, "read1", fs->n_read1);
    json_uint(&json, "read2", fs->n_read2);
    json_uint(&json, "read1_rmdup", fs->n_rmdup1);
    json_uint(&json, "read2_rmdup", fs->n_rmdup2);
    json_uint(&json, "forward_strand", fs->n_pstrand);
    json_uint(&json, "backward_strand", fs->n_mstrand);
    json_uint(&json, "pcr_duplicates", fs->n_dup);
    json_float(&json, "pcr_duplicates_fraction", (float)fs->n_dup / fs->n_mapped * 100);
    json_int(&json, "mapq_cutoff", f->mapQ_lim);
    json_uint(&json, "mapq_reads", fs->n_qual);
    json_float(&json, "mapq_reads_fraction_all", (float)fs->n_qual / fs->n_reads * 100);
    json_float(&json, "mapq_reads_fraction_mapped", (float)fs->n_qual / fs->n_mapped * 100);
    json_end_object(&json);
    kputc(',', &json);
    
    // Insert size section
    json_key(&json, "insert_size");
    json_start_object(&json);
    json_float(&json, "average", iavg);
    json_uint(&json, "median", imed);
    json_end_object(&json);
    kputc(',', &json);
    
    // Target section
    json_key(&json, "target");
    json_start_object(&json);
    json_uint(&json, "target_reads", fs->n_tgt);
    json_float(&json, "target_reads_fraction_all", (float)fs->n_tgt / fs->n_reads * 100);
    json_float(&json, "target_reads_fraction_mapped", (float)fs->n_tgt / fs->n_mapped * 100);
    json_float(&json, "target_data_mb", (float)fs->n_tdata / 1e6);
    json_float(&json, "target_data_rmdup_mb", (float)fs->n_trmdat / 1e6);
    json_float(&json, "target_data_fraction_all", (float)fs->n_tdata / fs->n_data * 100);
    json_float(&json, "target_data_fraction_mapped", (float)fs->n_tdata / fs->n_mdata * 100);
    json_uint(&json, "region_length", a->tgt_len);
    json_float(&json, "average_depth", (float)fs->n_tdata / a->tgt_len);
    json_float(&json, "average_depth_rmdup", (float)fs->n_trmdat / a->tgt_len);
    
    // Coverage sub-object
    json_key(&json, "coverage");
    json_start_object(&json);
    json_float(&json, "gt_0x", tarcov->cov);
    json_float(&json, "gte_4x", tarcov->cov4);
    json_float(&json, "gte_10x", tarcov->cov10);
    json_float(&json, "gte_30x", tarcov->cov30);
    json_float(&json, "gte_100x", tarcov->cov100);
    json_float(&json, "gt_0_2_avg", tarcov->cov02x);
    json_float(&json, "gt_0_5_avg", tarcov->cov05x);
    if (f->cutoff) {
        json_key(&json, "custom");
        json_start_object(&json);
        for (int x = 0; x < f->num_cutoffs; x++) {
            char key[32];
            sprintf(key, "gte_%dx", f->cutoffs[x]);
            json_float(&json, key, tarcov->cov_array[x]);
        }
        json_end_object(&json);
        kputc(',', &json);
    }
    if (f->depth_ratio && f->num_ratios > 0) {
        json_key(&json, "ratio_based");
        json_start_object(&json);
        for (int r = 0; r < f->num_ratios; r++) {
            char key[32];
            sprintf(key, "gt_%.1f_avg", f->ratios[r]);
            json_float(&json, key, tarcov->cov_ratio[r]);
        }
        json_end_object(&json);
        kputc(',', &json);
    }
    json_end_object(&json);
    kputc(',', &json);
    
    // Coverage rmdup sub-object
    json_key(&json, "coverage_rmdup");
    json_start_object(&json);
    json_float(&json, "gt_0x", tarcov->cov_rmdup);
    json_float(&json, "gte_4x", tarcov->cov4_rmdup);
    json_float(&json, "gte_10x", tarcov->cov10_rmdup);
    json_float(&json, "gte_30x", tarcov->cov30_rmdup);
    json_float(&json, "gte_100x", tarcov->cov100_rmdup);
    if (f->cutoff) {
        json_key(&json, "custom");
        json_start_object(&json);
        for (int x = 0; x < f->num_cutoffs; x++) {
            char key[32];
            sprintf(key, "gte_%dx", f->cutoffs[x]);
            json_float(&json, key, tarcov->cov_array_rmdup[x]);
        }
        json_end_object(&json);
        kputc(',', &json);
    }
    json_end_object(&json);
    kputc(',', &json);
    
    // Region coverage
    json_uint(&json, "region_count", a->tgt_nreg);
    json_key(&json, "region_coverage");
    json_start_object(&json);
    json_uint(&json, "gt_0x_count", regcov->cnt);
    json_float(&json, "gt_0x_fraction", regcov->cov);
    json_float(&json, "gte_4x_fraction", regcov->cov4);
    json_float(&json, "gte_10x_fraction", regcov->cov10);
    json_float(&json, "gte_30x_fraction", regcov->cov30);
    json_float(&json, "gte_100x_fraction", regcov->cov100);
    json_end_object(&json);
    
    json_end_object(&json);  // end target
    kputc(',', &json);
    
    // Flank section
    json_key(&json, "flank");
    json_start_object(&json);
    json_int(&json, "flank_size", flank_reg);
    json_uint(&json, "region_length", a->flk_len);
    json_float(&json, "average_depth", (float)fs->n_fdata / a->flk_len);
    json_uint(&json, "flank_reads", fs->n_flk);
    json_float(&json, "flank_reads_fraction_all", (float)fs->n_flk / fs->n_reads * 100);
    json_float(&json, "flank_reads_fraction_mapped", (float)fs->n_flk / fs->n_mapped * 100);
    json_float(&json, "flank_data_mb", (float)fs->n_fdata / 1e6);
    json_float(&json, "flank_data_fraction_all", (float)fs->n_fdata / fs->n_data * 100);
    json_float(&json, "flank_data_fraction_mapped", (float)fs->n_fdata / fs->n_mdata * 100);
    json_key(&json, "coverage");
    json_start_object(&json);
    json_float(&json, "gt_0x", flkcov->cov);
    json_float(&json, "gte_4x", flkcov->cov4);
    json_float(&json, "gte_10x", flkcov->cov10);
    json_float(&json, "gte_30x", flkcov->cov30);
    json_float(&json, "gte_100x", flkcov->cov100);
    json_end_object(&json);
    
    json_end_object(&json);  // end flank
    
    json_end_object(&json);  // end root
    
    // Write to file
    fprintf(fjson, "%s\n", json.s);
    fclose(fjson);
    free(json.s);
    
    return 0;
}

int print_report(struct opt_aux *f, aux_t *a, bamflag_t *fs)
{
    int i;
    if (outdir) chdir(outdir);
    FILE *finsert;
    FILE *fdep;
    finsert = open_wfile("insertsize.plot");
    fdep = open_wfile("depth_distribution.plot");
    // 设置更大的文件缓冲区以减少系统调用
    setvbuf(finsert, NULL, _IOFBF, WRITE_BUFFER_SIZE);
    setvbuf(fdep, NULL, _IOFBF, WRITE_BUFFER_SIZE);

    struct regcov *tarcov = regcov_init();
    struct regcov *flkcov = regcov_init();
    struct regcov *regcov = regcov_init();
    uint64_t icnt = 0;
    uint64_t dcnt = 0;
    for (i = 0; i < a->c_isize->m; ++i)
        icnt += a->c_isize->a[i];
    uint64_t icumu = icnt;
    for (i = 0; i < a->c_isize->m && i < f->isize_lim; ++i)
    {
        icumu -= a->c_isize->a[i];
        fprintf(finsert, "%d\t%u\t%f\t%" PRIu64 "\t%f\n", i, a->c_isize->a[i], (float)a->c_isize->a[i] / icnt, icumu,
                (float)icumu / icnt);
    }

    for (i = 0; i < a->c_dep->m; ++i)
        dcnt += a->c_dep->a[i];
    uint64_t dcumu = dcnt;
    for (i = 0; i < a->c_dep->m; ++i)
    {
        dcumu -= a->c_dep->a[i];
        fprintf(fdep, "%d\t%u\t%f\t%" PRIu64 "\t%f\n", i, a->c_dep->a[i], (float)a->c_dep->a[i] / dcnt, dcumu,
                (float)dcumu / dcnt);
    }
    // 计算 avg 和 med
    float iavg = average_cnt(a->c_isize);
    uint64_t imed = median_cnt(a->c_isize);

    fclose(fdep);
    fclose(finsert);
    // cntcov_cal(f, tarcov, a->c_dep, &fs->n_tdata);
    cntcov_cal2(f, tarcov, a->c_dep, &fs->n_tdata, a->tgt_len);
    cntcov_cal(f, regcov, a->c_reg, &fs->n_fdata);
    // merge_cnt(a->c_flkdep, a->c_dep);
    cntcov_cal(f, flkcov, a->c_flkdep, &fs->n_fdata);
    for (i = 0; i < a->c_rmdupdep->m; ++i)
        fs->n_trmdat += a->c_rmdupdep->a[i] * i;
    
    // 计算 rmdup depth 覆盖度统计
    uint64_t rmdup_data;
    cntcov_cal_rmdup(f, tarcov, a->c_rmdupdep, &rmdup_data, a->tgt_len);
    
    // 计算基于比例的覆盖度统计
    cntcov_cal_ratio(f, tarcov, a->c_dep, a->tgt_len);

    FILE *fchrcov = open_wfile("chromosomes.report");
    // 设置更大的文件缓冲区
    setvbuf(fchrcov, NULL, _IOFBF, WRITE_BUFFER_SIZE);
    {
        fprintf(fchrcov, "%11s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s", "#Chromosome", "DATA(%)", "Avg depth",
                "Median", "Coverage%", "Cov 4x %", "Cov 10x %", "Cov 30x %", "Cov 100x %");
        if (f->cutoff)
        {
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                fprintf(fchrcov, "\tCov %dx %%", f->cutoffs[x]);
            }
        }
        fprintf(fchrcov, "\n");
        khiter_t k;
        for (k = 0; k != kh_end(a->h_tgt); ++k)
        {
            if (kh_exist(a->h_tgt, k))
            {
                char *name = (char *)kh_key(a->h_tgt, k);
                count32_t *cnt = (count32_t *)kh_val(a->h_tgt, k).data;
                uint64_t data = 0;
                struct regcov *chrcov = regcov_init();
                cntcov_cal(f, chrcov, cnt, &data);
                uint64_t length = 0;
                int i;
                for (i = 0; i < cnt->m; ++i)
                    length += cnt->a[i];
                float avg, med, per;
                if (data > 0)
                {
                    avg = (float)data / length;
                    med = median_cnt(cnt);
                    per = (float)data / fs->n_tdata * 100.0;
                    fprintf(fchrcov, "%11s\t%8.2f\t%8.2f\t%9.1f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f", name, per, avg,
                            med, chrcov->cov, chrcov->cov4, chrcov->cov10, chrcov->cov30, chrcov->cov100);
                    if (f->cutoff)
                    {
                        for (int x = 0; x < f->num_cutoffs; x++)
                        {
                            fprintf(fchrcov, "\t%.2f", chrcov->cov_array[x]);
                        }
                    }
                }
                else
                {
                    fprintf(fchrcov, "%11s\t%8.2f\t%8.2f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f", name, (float)0,
                            (float)0, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0);
                    if (f->cutoff)
                    {
                        for (int x = 0; x < f->num_cutoffs; x++)
                        {
                            fprintf(fchrcov, "\t%5.4f", (float)0);
                        }
                    }
                }
                fprintf(fchrcov, "\n");
                free(chrcov);
            }
        }
    }
    fclose(fchrcov);
    FILE *fc = open_wfile("coverage.report");
    // 设置更大的文件缓冲区以减少系统调用
    setvbuf(fc, NULL, _IOFBF, WRITE_BUFFER_SIZE);
    do
    {
        fprintf(fc, "## The file was created by %s\n", program_name);
        fprintf(fc, "## Version : %s\n", Version);
        fprintf(fc, "## Files : ");
        for (i = 0; i < f->nfiles; ++i)
            fprintf(fc, "%s ", f->inputs[i]);
        fprintf(fc, "\n");
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Raw Reads (All reads)", fs->n_reads);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] QC Fail reads", fs->n_qcfail);
        fprintf(fc, "%60s\t%.2f\n", "[Total] Raw Data(Mb)", (float)fs->n_data / 1e6);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Paired Reads", fs->n_pair_all);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Mapped Reads", fs->n_mapped);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Mapped Reads", (float)fs->n_mapped / fs->n_reads * 100);
        fprintf(fc, "%60s\t%.2f\n", "[Total] Mapped Data(Mb)", fs->n_mdata / 1e6);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Mapped Data(Mb)", (float)fs->n_mdata / fs->n_data * 100);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Properly paired", fs->n_pair_good);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Properly paired",
                (float)fs->n_pair_good / fs->n_reads * 100);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read and mate paired", fs->n_pair_map);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of Read and mate paired",
                (float)fs->n_pair_map / fs->n_reads * 100);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Singletons", fs->n_sgltn);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read and mate map to diff chr", fs->n_diffchr);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read1", fs->n_read1);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read2", fs->n_read2);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read1(rmdup)", fs->n_rmdup1);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] Read2(rmdup)", fs->n_rmdup2);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] forward strand reads", fs->n_pstrand);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] backward strand reads", fs->n_mstrand);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] PCR duplicate reads", fs->n_dup);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of PCR duplicate reads",
                (float)fs->n_dup / fs->n_mapped * 100); // change n_reads to n_mapped, 2015/05/25
        fprintf(fc, "%60s\t%d\n", "[Total] Map quality cutoff value", f->mapQ_lim);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Total] MapQuality above cutoff reads", fs->n_qual);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of MapQ reads in all reads",
                (float)fs->n_qual / fs->n_reads * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[Total] Fraction of MapQ reads in mapped reads",
                (float)fs->n_qual / fs->n_mapped * 100);
        // insert
        fprintf(fc, "%60s\t%.2f\n", "[Insert size] Average", iavg);
        fprintf(fc, "%60s\t%.ld\n", "[Insert size] Median", imed);
        // tgt
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Target] Target Reads", fs->n_tgt);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Reads in all reads",
                (float)fs->n_tgt / fs->n_reads * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Reads in mapped reads",
                (float)fs->n_tgt / fs->n_mapped * 100);
        fprintf(fc, "%60s\t%.2f\n", "[Target] Target Data(Mb)", (float)fs->n_tdata / 1e6);
        fprintf(fc, "%60s\t%.2f\n", "[Target] Target Data Rmdup(Mb)", (float)fs->n_trmdat / 1e6);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Data in all data",
                (float)fs->n_tdata / fs->n_data * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction of Target Data in mapped data",
                (float)fs->n_tdata / fs->n_mdata * 100);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Target] Len of region", a->tgt_len);
        fprintf(fc, "%60s\t%.2f\n", "[Target] Average depth", (float)fs->n_tdata / a->tgt_len);
        fprintf(fc, "%60s\t%.2f\n", "[Target] Average depth(rmdup)", (float)fs->n_trmdat / a->tgt_len);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>0.2*(Average depth)x)", tarcov->cov02x);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>0.5*(Average depth)x)", tarcov->cov05x);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>0x)", tarcov->cov);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=4x)", tarcov->cov4);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=10x)", tarcov->cov10);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=30x)", tarcov->cov30);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage (>=100x)", tarcov->cov100);
        if (f->cutoff)
        {
            char titles[40];
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                sprintf(titles, "[Target] Coverage (>=%ux)", f->cutoffs[x]);
                fprintf(fc, "%60s\t%.2f%%\n", titles, tarcov->cov_array[x]);
            }
        }
        // rmdup coverage statistics
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage(rmdup) (>0x)", tarcov->cov_rmdup);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage(rmdup) (>=4x)", tarcov->cov4_rmdup);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage(rmdup) (>=10x)", tarcov->cov10_rmdup);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage(rmdup) (>=30x)", tarcov->cov30_rmdup);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Coverage(rmdup) (>=100x)", tarcov->cov100_rmdup);
        if (f->cutoff)
        {
            char titles[60];
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                sprintf(titles, "[Target] Coverage(rmdup) (>=%ux)", f->cutoffs[x]);
                fprintf(fc, "%60s\t%.2f%%\n", titles, tarcov->cov_array_rmdup[x]);
            }
        }
        // ratio-based coverage statistics
        if (f->depth_ratio && f->num_ratios > 0)
        {
            char titles[60];
            for (int r = 0; r < f->num_ratios; r++)
            {
                sprintf(titles, "[Target] Coverage (>%.1f*Avg)", f->ratios[r]);
                fprintf(fc, "%60s\t%.2f%%\n", titles, tarcov->cov_ratio[r]);
            }
        }
        // tgt regions
        fprintf(fc, "%60s\t%u\n", "[Target] Target Region Count", a->tgt_nreg);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[Target] Region covered > 0x", regcov->cnt);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered > 0x", regcov->cov);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 4x", regcov->cov4);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 10x", regcov->cov10);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 30x", regcov->cov30);
        fprintf(fc, "%60s\t%.2f%%\n", "[Target] Fraction Region covered >= 100x", regcov->cov100);
        if (f->cutoff)
        {
            char titles[60];
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                sprintf(titles, "[Target] Fraction Region covered (>=%ux)", f->cutoffs[x]);
                fprintf(fc, "%60s\t%.2f%%\n", titles, regcov->cov_array[x]);
            }
        }
        // tgt regions

        // flk
        fprintf(fc, "%60s\t%u\n", "[flank] flank size", flank_reg);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[flank] Len of region (not include target region)", a->flk_len);
        fprintf(fc, "%60s\t%.2f\n", "[flank] Average depth", (float)fs->n_fdata / a->flk_len);
        fprintf(fc, "%60s\t%" PRIu64 "\n", "[flank] flank Reads", fs->n_flk);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Reads in all reads",
                (float)fs->n_flk / fs->n_reads * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Reads in mapped reads",
                (float)fs->n_flk / fs->n_mapped * 100);
        fprintf(fc, "%60s\t%.2f\n", "[flank] flank Data(Mb)", (float)fs->n_fdata / 1e6);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Data in all data",
                (float)fs->n_fdata / fs->n_data * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Fraction of flank Data in mapped data",
                (float)fs->n_fdata / fs->n_mdata * 100);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>0x)", flkcov->cov);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=4x)", flkcov->cov4);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=10x)", flkcov->cov10);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=30x)", flkcov->cov30);
        fprintf(fc, "%60s\t%.2f%%\n", "[flank] Coverage (>=100x)", flkcov->cov100);
        if (f->cutoff)
        {
            char titles[40];
            for (int x = 0; x < f->num_cutoffs; x++)
            {
                sprintf(titles, "[flank] Coverage (>=%ux)", f->cutoffs[x]);
                fprintf(fc, "%60s\t%.2f%%\n", titles, flkcov->cov_array[x]);
            }
        }

    } while (0);

    // 生成 JSON 格式报告
    print_report_json(f, a, fs, tarcov, flkcov, regcov, iavg, imed);

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
    DEPTHRATIO,
    REFERENCE,
    THREADS,
    HELP
};

static struct option const long_opts[] = {{"outdir", required_argument, NULL, 'o'},
                                          {"bed", required_argument, NULL, 'p'},
                                          {"flank", required_argument, NULL, 'f'},
                                          {"maxdepth", required_argument, NULL, MAXDEPTH},
                                          {"cutoffdepth", required_argument, NULL, CUTOFF},
                                          {"isize", required_argument, NULL, INSERTSIZE},
                                          {"mapthres", required_argument, NULL, 'q'},
                                          {"uncover", required_argument, NULL, UNCOVER},
                                          {"bamout", required_argument, NULL, BAMOUT},
                                          {"depthratio", required_argument, NULL, DEPTHRATIO},
                                          {"reference", required_argument, NULL, 'T'},
                                          {"threads", required_argument, NULL, THREADS},
                                          //{"rmdup", no_argument, NULL, 'd'},
                                          {"help", no_argument, NULL, 'h'},
                                          {"version", no_argument, NULL, 'v'}};

int show_version()
{
    printf("%s\n", Version);
    return 1;
}
int xamdst(int argc, char *argv[])
{
    int n, i;
    char *probe = 0;

    // struct opt_aux opt = {.inputs = NULL, .isize_lim = 2000, .mapQ_lim = 20};
    struct opt_aux opt = init_opt_aux();
    while ((n = getopt_long(argc, argv, "o:p:f:q:l:T:h1v", long_opts, NULL)) >= 0)
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
        case DEPTHRATIO:
            // 解析用户指定的深度比例
            {
                char *ratio_token = strtok(optarg, ",");
                opt.num_ratios = 0;  // 重置为用户指定的值
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
                        fprintf(stderr, "Invalid depth ratio: %s (must be between 0 and 1)\n", ratio_token);
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
        case 'T':
            // 参考基因组路径（用于 CRAM 文件）
            opt.reference = strdup(optarg);
            break;
        case THREADS:
            opt.nthreads = atoi(optarg);
            if (opt.nthreads < 0)
            {
                fprintf(stderr, "Invalid thread count: %s\n", optarg);
                exit(EXIT_FAILURE);
            }
            break;
        case 'h':
            usage(1);
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
    if (isNull(outdir) || isNull(probe))
        usage(1);
    if (export_target_bam && check_filename_isbam(export_target_bam))
    {
        fprintf(stderr, "--bamout must be a bam file: %s", export_target_bam);
        goto freeall;
    }

    n = argc - optind;
    mkdirp(outdir, 0755);
    // capable of deals with several bam/cram files
    aux_t *aux;
    aux = aux_init();
    if (isZero(n))
    {
        aux->data = (htsFile **)needmem(sizeof(htsFile *));
        aux->data[0] = hts_open("-", "r");
        if (aux->data[0] == NULL)
            errabort("Failed to open stdin for reading");
        // Enable multi-threading for stdin reading
        if (opt.nthreads > 0)
        {
            if (hts_set_threads(aux->data[0], opt.nthreads) != 0)
            {
                fprintf(stderr, "Warning: Failed to set %d threads for stdin\n", opt.nthreads);
            }
        }
        // Set reference for CRAM if provided
        if (opt.reference)
        {
            if (hts_set_fai_filename(aux->data[0], opt.reference) != 0)
            {
                errabort("Failed to set reference file: %s", opt.reference);
            }
        }
        aux->h = sam_hdr_read(aux->data[0]);
        if (aux->h == NULL)
            errabort("Failed to read header from stdin");
        aux->ndata = 1;
        opt.nfiles = 0;
    }
    else
    {
        aux->data = (htsFile **)needmem(n * sizeof(htsFile *));
        opt.nfiles = n;
        opt.inputs = (char **)needmem(n * sizeof(char *));
        for (i = 0; i < n; ++i)
        {
            sam_hdr_t *h_tmp;
            if (STREQ(argv[optind + i], "-"))
            {
                aux->data[i] = hts_open("-", "r");
                stdin_lock = 1;
            }
            else
            {
                aux->data[i] = hts_open(argv[optind + i], "r");
            }
            if (aux->data[i] == NULL)
                errabort("%s: %s", argv[optind + i], strerror(errno));

            // Enable multi-threading for BAM/CRAM reading
            if (opt.nthreads > 0)
            {
                if (hts_set_threads(aux->data[i], opt.nthreads) != 0)
                {
                    fprintf(stderr, "Warning: Failed to set %d threads for %s\n", opt.nthreads, argv[optind + i]);
                }
            }

            // Check if CRAM and reference is needed
            const htsFormat *fmt = hts_get_format(aux->data[i]);
            if (fmt->format == cram)
            {
                if (opt.reference == NULL)
                {
                    errabort("CRAM file requires a reference genome. Use -T/--reference to specify.");
                }
                if (hts_set_fai_filename(aux->data[i], opt.reference) != 0)
                {
                    errabort("Failed to set reference file: %s", opt.reference);
                }
            }
            else if (opt.reference)
            {
                // Also set reference for BAM if provided (useful for some operations)
                hts_set_fai_filename(aux->data[i], opt.reference);
            }
            
            h_tmp = sam_hdr_read(aux->data[i]);
            if (h_tmp == NULL)
                errabort("Failed to read header from %s", argv[optind + i]);
            if (i == 0)
                aux->h = h_tmp;
            else
                sam_hdr_destroy(h_tmp);
            opt.inputs[i] = strdup(argv[optind + i]);
        }
        aux->ndata = n;
    }
    // 多 BAM/CRAM 文件支持说明：
    // 当前实现已支持多个 BAM/CRAM 文件输入，它们会被顺序处理并合并统计结果
    // 注意：所有文件必须使用相同的参考基因组和排序方式
    // 第一个文件的 header 会被用于输出
    if (export_target_bam)
    {
        bamoutfp = hts_open(export_target_bam, "wb");
        if (bamoutfp == NULL)
            errabort("%s : %s", export_target_bam, strerror(errno));
        // Enable multi-threading for BAM writing
        if (opt.nthreads > 0)
        {
            if (hts_set_threads(bamoutfp, opt.nthreads) != 0)
            {
                fprintf(stderr, "Warning: Failed to set %d threads for output BAM\n", opt.nthreads);
            }
        }
        // Write BAM header using htslib compatible format
        if (sam_hdr_write(bamoutfp, aux->h) < 0)
            errabort("Failed to write header to %s", export_target_bam);
    }
    h_chrlength_init();
    header2chrhash(aux->h);
    load_bed_init(probe, aux);
    chrhash_destroy();
    freemem(probe);
    if (aux->c_isize->n < opt.isize_lim)
    {
        aux->c_isize->a = realloc(aux->c_isize->a, opt.isize_lim * sizeof(unsigned));
    }
    for (i = 0; i < opt.isize_lim; ++i)
        aux->c_isize->a[i] = 0;
    // aux->c_isize->m = opt.isize_lim;
    aux->nchr = sam_hdr_nref(aux->h);
    struct bamflag fs = {};
    load_bamfiles(&opt, aux, &fs);
    print_report(&opt, aux, &fs);
    aux_destroy(aux);
    for (i = 0; i < opt.nfiles; ++i)
        freemem(opt.inputs[i]);
    freemem(opt.inputs);
    if (export_target_bam)
        hts_close(bamoutfp);
    freemem(opt.cutoffs);
    freemem(opt.ratios);
    freemem(opt.reference);
freeall:
    freemem(export_target_bam);
    freemem(outdir);
    return 0;
}

/* main */
int main(int argc, char *argv[])
{
    return xamdst(argc, argv);
}
