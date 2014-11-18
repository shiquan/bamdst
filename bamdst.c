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

#include "commons.h"
#include "count.h"
#include "bam.h"
#include "sam_header.h"
#include "bedutil.h"
#include "khash.h"
#include "kstring.h"
#include "knetfile.h"

static char const *program_name = "bamdst";
static char const *Version = "1.0.0";

static int flank_reg = 200;

extern bedHandle_t *bedHand;

#if !defined OPEN_MAX && defined NR_OPEN
# define OPEN_MAX NR_OPEN
#endif
#if !defined OPEN_MAX
# define OPEN_MAX 20
#endif

/* only accepted one stdin pipeline */
static bool stdin_lock = FALSE;

/* The number of threads after which there are 
   diminishing performance gains. */
enum { DEFAULT_MAX_THREADS = 8 };

/* duplicate will be removed if rmdup_mark is TRUE */
static bool rmdup_mark = FALSE;

struct opt_aux
	{
	int nfiles;
	char *outdir;
	char **inputs;
	//	int flk_reg;
	int isize_lim;
	int cutoff;
	int mapQ_lim;
	int maxdepth;
	};

/* @FLANK_REGION coverage of flank region list in report.
   @INSERTSIZE_LIMIT the insert size bigger than this will not be calculated.
   @MAXDEPTH_LIMIT the depth bigger than this will not be calculated, 
                   0 means no limit.
   @CUTOFF  if you want know coverage of specity depth, set this value.
   @MAPQ_LIMIT  >=mapQ_limit list in report.
*/

static struct opt_aux *
opt_init()
	{
	struct opt_aux *aux;
	aux = (struct opt_aux*)needmem(sizeof(struct opt_aux));
	aux->nfiles = 0;
	aux->outdir = NULL;
	aux->inputs = NULL;
	//aux->flk_reg = 200;
	aux->isize_lim = 2000;
	aux->cutoff = 0; //
	aux->maxdepth = 0;
	aux->mapQ_lim = 20;
	return aux;
	}

struct depnode
	{
	unsigned len;
	unsigned start;
	unsigned stop;
	unsigned *vals;
	unsigned *cnts;
	struct depnode *next;
	};

void opt_destroy(struct opt_aux *opt)
	{
	freemem(opt->outdir);
	int i;
	for (i = 0; i < opt->nfiles; ++i) freemem(opt->inputs[i]);
	freemem(opt);
	}

static const char * init_debugmsg[] =
	{
	"Success",
	"Trying to allocated an unempty node",
	"Trying to allocated a zero memory"
	"END of list"
	};

#define INIT_DEBUG(x) do {											\
 int _a = x;														\
 if(_a)																\
	 {																\
	 warnings("%s : %d %s",__FILE__, __LINE__, init_debugmsg[_a]);	\
	 }																\
 } while(0)

/* node->vals is the depth value of each loc 
 * node->cnts is the count of covered reads  
 * init the memory before use it */
static int 
depnode_init(struct depnode *node)
	{
	if ( isNull(node) ) return 3;
	if ( node->len ) return 1;
	node->len = node->stop - node->start + 1;
	if ( isZero(node->len)) return 2;
	node->vals = (unsigned*)needmem((node->len+1) * sizeof(unsigned));
	node->cnts = (unsigned*)needmem((node->len+1) * sizeof(unsigned));
	memset(node->vals, 0, sizeof *(node->vals));
	memset(node->cnts, 0, sizeof *(node->cnts));
	return 0;
	};

/* delete node and make the node point to the next node */
#define del_node(node) do {							\
 if (node)											\
	 {												\
	 struct depnode *tmpnode = node;				\
	 node = node->next;								\
	 freemem(tmpnode->vals);						\
	 freemem(tmpnode->cnts);						\
	 freemem(tmpnode);								\
	 }												\
 } while(0)

static struct depnode *
bed_depnode_list(bedreglist_t *bed)
	{
	struct depnode *node;
	struct depnode *header;
	struct depnode *tmpnode;
	int i;
	for ( i = 0; i < bed->m; ++i)
		{
		node = (struct depnode*)needmem(sizeof(struct depnode));
		node->start = (uint32_t)(bed->a[i] >> 32);
		node->stop = (uint32_t)bed->a[i];
		node->len = 0; //(uint32_t)bed->a[i] - node->start + 1;
		node->vals = NULL;
		node->cnts = NULL;
		if ( isZero(i)) header = node;
		else tmpnode->next = node;
		tmpnode = node;
		}
	INIT_DEBUG(depnode_init(header));
	return header;
	}

struct _aux
	{
	int nchr, ndata, maxdep;
	uint64_t tgt_len, flk_len;
	unsigned tgt_nreg;
	bamFile *data;
	bam_header_t *h;
	regHash_t *h_tgt;
	regHash_t *h_flk;
	count32_t *c_dep;
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
	count32_init(a->c_isize);
	count32_init(a->c_flkdep);
	count32_init(a->c_reg);
	return a;
	}

void aux_destroy(struct _aux *a)
	{
	int i;
	for (i = 0; i < a->ndata; ++i) bam_close(a->data[i]);
	free(a->data);
	bedHand->destroy((void *)a->h_tgt, destroy_void);
	bedHand->destroy((void *)a->h_flk, destroy_void);
	bam_header_destroy(a->h);
	count_destroy(a->c_dep);
	count_destroy(a->c_flkdep);
	count_destroy(a->c_isize);
	count_destroy(a->c_reg);
	free(a);
	}

typedef struct
    {
    uint64_t n_reads, n_mapped, n_pair_map, n_pair_all, n_pair_good;
    uint64_t n_sgltn, n_read1, n_read2;
    uint64_t n_dup;
    uint64_t n_diffchr, n_pstrand, n_mstrand;
    uint64_t n_qcfail;
    uint64_t n_data, n_mdata;
    uint64_t n_qual;
    /* delete n_uniq
	 * ref:  https://www.biostars.org/p/59281/ */
    uint64_t n_tgt, n_flk, n_tdata, n_fdata;
    }
bamflag_t;

bamflag_t *bamflag_init()
	{
	bamflag_t *fs = (bamflag_t*)needmem(sizeof(bamflag_t));
	fs->n_reads = fs->n_mapped = fs->n_pair_good = fs->n_pair_good = 0;
	fs->n_sgltn = fs->n_read1 = fs->n_read2 = fs->n_pair_map = 0;
	fs->n_dup = fs->n_diffchr = fs->n_pstrand = fs->n_mstrand = 0;
	fs->n_qcfail = fs->n_data = fs->n_mdata = fs->n_qual = 0;
	fs->n_tgt = fs->n_flk = fs->n_tdata = fs->n_fdata = 0;
	return fs;
	}

#define flagstat(s, c, ret) do {										\
 ++(s)->n_reads;														\
 (s)->n_data += (c)->l_qseq;											\
 if ((c)->flag & BAM_FQCFAIL)											\
	 {																	\
	 ++(s)->n_qcfail;													\
	 ret = 0;															\
	 }																	\
 else																	\
	 {																	\
	 ret = 1;															\
	 if ((c)->flag & BAM_FPAIRED)										\
		 {																\
		 ++(s)->n_pair_all;												\
		 if ((c)->flag & BAM_FPROPER_PAIR) ++(s)->n_pair_good;			\
		 if ((c)->flag & BAM_FREAD1) ++(s)->n_read1;					\
		 if ((c)->flag & BAM_FREAD2) ++(s)->n_read2;					\
		 if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP))	\
			 ++(s)->n_sgltn;											\
		 if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP))	\
			 {															\
			 ++(s)->n_pair_map;											\
			 if ((c)->mtid != (c)->tid) ++(s)->n_diffchr;				\
			 }															\
		 }																\
	 if (!((c)->flag & BAM_FUNMAP))										\
		 {																\
		 ++(s)->n_mapped;												\
		 (s)->n_mdata += (c)->l_qseq;									\
		 }																\
	 if ((c)->flag & BAM_FDUP)											\
		 {																\
		 ++(s)->n_dup;													\
		 ret = 2;														\
		 }																\
	 if ((c)->flag & BAM_FSECONDARY) ret = 3;							\
	 if ((c)->flag & BAM_FREVERSE) ++(s)->n_mstrand;					\
	 else ++(s)->n_pstrand;												\
	 }																	\
 } while(0)

static void emit_try_help(void)
	{
	fprintf (stderr, "Try '%s --help' for more information.\n", program_name);
	}

void usage(int status)
	{
	if (status == 0)
		emit_try_help();
	else
		{
		printf ("\
USAGE : %s [OPTION] -p <probe.bed> -o <output_dir> [in1.bam [in2.bam ... ]]\n\
   or : %s [OPTION] -p <probe.bed> -o <output_dir> -\n\
",
				program_name, program_name);
		puts ("\
Option -o and -p are mandatory:\n\
  -o, --outdir         output dir\n\
  -p, --bed            probe or target regions file, the region file will \n\
                       be merged before calculate depths\n\
");
		puts ("\
Ordering options:\n\
   -f, --flank [200]   flank n bp of each region\n\
   --maxdepth[0]       set the max depth to stat the cumu distribution.\n\
   --cutoffdepth       list the coverage of above depths\n\
   -h, --help          print this help info\n\
\n");
		/*-d, --rmdup         remove dup reads when calculate depth\n	\*/
		puts ("\
* Five essential files would be created in the output file :\n\n\
 - coverage.report     a report of the coverage information and reads \n\
                       information of whole target regions\n\
 - target.detail       print depth value of each base in fasta-like format\n\
 - cumu.plot           distribution data of depth values\n\
 - insert.plot         distribution data of inferred insert size \n\
 - chromosome.report   coverage information for each chromosome\n\
 - exon.tsv            mean depth, median depth and coverage of each region\n\
 - depth.tsv.gz        raw depth, rmdup depth, coverage depth of each position\n\
");
		}
	exit(EXIT_SUCCESS);
	}

#include"ksort.h"

KSORT_INIT_GENERIC(uint32_t)

static float median_cal(const uint32_t * array, int l)
	{
	if (isNull(l)) return 0;
	float med = 0;
	uint32_t *tmp;
	tmp = (uint32_t*)needmem(l *sizeof(uint32_t));
	memcpy(tmp, array, l * sizeof(uint32_t));
	ks_introsort(uint32_t, l, tmp);
	med = l & 1 ? tmp[(l >> 1) + 1] : (float)(tmp[l >> 1] + tmp[(l >> 1) - 1]) / 2;
	free(tmp);
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
	return avg;
	}

static float coverage_cal(const uint32_t * array, int l)
	{
	if (isNull(l)) return 0;
	float cov = 0;
	int i;
	for (i = 0; i < l; ++i)
		if (array[i]) cov++;
	cov /= (float)l;
	return cov;
	}

const static char *report_total[] =
	{
	"[Total] Raw Reads (All reads)", "[Total] QC Fail reads",
	"[Total] Raw Data(Mb)",
	"[Total] Paired Reads", "[Total] Mapped Reads",
	"[Total] Fraction of Mapped Reads",
	"[Total] Mapped Data(Mb)", "[Total] Fraction of Mapped Data(Mb)",
	"[Total] Properly paired",
	"[Total] Fraction of Properly paired", "[Total] Read and mate paired",
	"[Total] Fraction of Read and mate paired", "[Total] Singletons",
	"[Total] Read and mate map to diff chr", "[Total] Read1", "[Total] Read2",
	"[Total] forward strand reads", "[Total] backward strand reads",
	"[Total] PCR duplicate reads", "[Total] Fraction of PCR duplicate reads",
	"[Total] Map quality cutoff value",
	"[Total] MapQuality above cutoff reads",
	"[Total] Fraction of MapQ reads in all reads",
	"[Total] Fraction of MapQ reads in mapped reads"
	};

const static char *report_tar[] =
	{
	"[Target] Target Reads", "[Target] Fraction of Target Reads in all reads",
	"[Target] Fraction of Target Reads in mapped reads",
	"[Target] Target Data(Mb)",
	"[Target] Fraction of Target Data in all data",
	"[Target] Fraction of Target Data in mapped data",
	"[Target] Len of region", "[Target] Average depth",
	"[Target] Coverage (>0x)",
	"[Target] Coverage (>=4x)", "[Target] Coverage (>=10x)",
	"[Target] Coverage (>=30x)", "[Target] Coverage (>=100x)",
	"[Target] Target Region Count",
	"[Target] Region covered > 0x",
	"[Target] Fraction Region covered > 0x",
	//"[Target] Region covered >= 4x",
	"[Target] Fraction Region covered >= 4x",
	//"[Target] Region covered >= 10x",
	"[Target] Fraction Region covered >= 10x",
	//"[Target] Region covered >= 30x",
	"[Target] Fraction Region covered >= 30x",
	//"[Target] Region covered >= 100x",
	"[Target] Fraction Region covered >= 100x",
	"[flank] flank size", "[flank] Len of region", "[flank] Average depth",
	"[flank] flank Reads", "[flank] Fraction of flank Reads in all reads",
	"[flank] Fraction of flank Reads in mapped reads",
	"[flank] flank Data(Mb)",
	"[flank] Fraction of flank Data in all data",
	"[flank] Fraction of flank Data in mapped data", "[flank] Coverage (>0x)",
	"[flank] Coverage (>=4x)", "[flank] Coverage (>=10x)",
	"[flank] Coverage (>=30x)", "[flank] Coverage (>=100x)"
	};

// FIXME: need broken when bed file is truncated
int load_bed_init(char const *fn, aux_t * a)
	{
	bedHand->read(fn, a->h_tgt, 0, 0);
	bedHand->merge(a->h_tgt);
	inf_t *inf1 = bedHand->stat(a->h_tgt);
	a->tgt_len = inf1->length;
	a->tgt_nreg = inf1->total;
	bedHand->read(fn, a->h_flk, flank_reg, flank_reg);
	bedHand->merge(a->h_flk);
	bedHand->diff(a->h_flk, a->h_tgt);
	inf_t *inf2 = bedHand->stat(a->h_flk);
	a->flk_len = inf2->length;
	mustfree(inf1);
	mustfree(inf2);
	return 1;
	}

int match_pos(struct depnode * header, uint32_t pos)
	{
	struct depnode *tmp = header;
	while (tmp && pos > tmp->start+ tmp->len)
		{
		tmp = tmp->next;
		if (tmp && isZero(tmp->len)) depnode_init(tmp);
		}
	if (isNull(tmp)) return 1; // this chromosome is finished, skip in next loop
	//if (pos > tmp->start + tmp->len - 1) return 1;
	if (pos >= tmp->start) tmp->vals[pos - tmp->start]++;
	return 0;
	}

int readcore(struct depnode * header, bam1_t * b)
	{
	struct depnode *tmp = header;
	int i;
	bam1_core_t *c = &b->core;
	if (c->pos < tmp->start) return 0;
	while (tmp && c->pos > tmp->stop)
		{
		tmp = tmp->next;
		}
	if (isNull(tmp)) return 0;
	uint32_t *cigar = bam1_cigar(b);
	uint32_t end = bam_calend(c, cigar);
	if (end >= tmp->start - 1)
		{
		int j, pos, l, s;
		pos = c->pos + 1;
		if (pos > tmp->start)
			tmp->cnts[pos - tmp->start]++;
		for (i = 0; i < c->n_cigar; ++i)
			{
			s = cigar[i] & 0xf;
			l = cigar[i] >> BAM_CIGAR_SHIFT;
			if (s == BAM_CMATCH)
				for (j = 0; j < l; ++j)
					{
					if (match_pos(tmp, pos))
						{
						pos += l - j;
						break;
						}
					pos++;
					}
			else if (s == BAM_CDEL)
				pos++;
			}
		return 1;
		}
	return 0;
	}

int print_tgtdep(char *name, FILE *fp, struct depnode *node, aux_t *a)
	{
	if (isNull(node->vals))
		{
		errabort("node->vals is null");
		}
	int j;
	float avg, med, cov;
	avg = avg_cal(node->vals, node->len);
	med = median_cal(node->vals, node->len);
	cov = coverage_cal(node->vals, node->len);
	count_increase(a->c_reg, (int)avg, uint32_t);
	fprintf(fp, ">%s:%d-%d:NA:NA\t[%.2f,%.1f,%.4f]\n",
			name, node->start, node->stop, avg, med, cov);
	for (j = 0; j < node->len; ++j)
		{
		fprintf(fp, "%u ", node->vals[j]);
		// count_increase will alloc memory space automatically
		count_increase(a->c_dep, node->vals[j], uint32_t);
		}
	putc('\n', fp);
	return 1;
	}

int print_zerotgtdep(char *name, FILE *fp, struct depnode *node, aux_t *a)
	{
	int j;
	fprintf(fp, ">%s:%d-%d:NA:NA\t[0,0,0]\n",
			name, node->start, node->stop);
	for (j = 0; j < node->len; ++j)
		fputs("0 ", fp);
	putc('\n', fp);
	count_increaseN(a->c_dep, 0, node->len, uint32_t);
	count_increase(a->c_reg, 0, uint32_t);
	return 1;
	}
	
int load_bamfiles(struct opt_aux *f, aux_t * a, bamflag_t * fs)
	{
	int i, j;
	int ret;
	int tid = -1;
	int lstpos = 0;
	bam_header_t *h = a->h;
	bedreglist_t *tar;
	bedreglist_t *flk;
	if (f->outdir) chdir(f->outdir);
	FILE *fp = open_wfile("target.dep");
	struct depnode *tgt_node = NULL;
	struct depnode *flk_node = NULL;
	char *name;
	int dbgret;
	for (i = 0; i < a->ndata; ++i)
		{
		bamFile dat = a->data[i];
		bool goto_next_chromosome = FALSE;
		while (1)
			{
			bam1_t *b;
			b = (bam1_t*)needmem(sizeof(bam1_t));
			ret = bam_read1(dat, b);
			if (ret == -1)
				{
				mustfree(b); break;				//normal end
				}
			if (ret == -2) errabort("%d bam file is truncated!\n", i + 1);
			bam1_core_t *c = &b->core;
			flagstat(fs, c, ret);
			if (c->qual > f->mapQ_lim) fs->n_qual++;
			if (rmdup_mark && ret > 1) goto endcore;
			if (c->tid == -1) goto endcore;
			/* stat the insertsize */
			if (c->isize > 0 && c->isize < f->isize_lim)
				{
				count_increase(a->c_isize, c->isize, uint32_t);
				}
			
			if (tid == c->tid )
				{
				if (goto_next_chromosome) goto endcore;
				if (tid == c->tid && lstpos > c->pos)
					errabort("The bam file is not sorted!");
				}
			lstpos = c->pos;
			if (tid != c->tid) // FIXME: need multi thread to improve it
				{
				goto_next_chromosome = FALSE;
				if (tgt_node && tid >= 0)
					{
					print_tgtdep(name, fp, tgt_node, a);
					del_node(tgt_node);
					INIT_DEBUG(depnode_init(tgt_node));
					while (tgt_node)
						{
						print_zerotgtdep(name, fp, tgt_node, a);
						del_node(tgt_node);
						INIT_DEBUG(depnode_init(tgt_node));
						}
					}
				if (flk_node && tid >= 0)
					{
					int j;
					for (j = 0; j < flk_node->len; ++j)
						count_increase(a->c_flkdep, flk_node->vals[j], uint32_t);
					del_node(flk_node);
					INIT_DEBUG(depnode_init(flk_node));
					
					while (flk_node)
						{
						count_increaseN(a->c_flkdep, 0, flk_node->len, uint32_t);
						//a->c_flkdep->a[0] += flk_node->len;
						del_node(flk_node);
						INIT_DEBUG(depnode_init(flk_node));						
						}
					}
				tid = c->tid;
				name = h->target_name[tid];
				debug("%s",name);
				if (tid > a->nchr)
					errabort("chromosome %s is not in bam header!"
							 "It must be use different bam headers"
							, name);
				khiter_t k;
				k = kh_get(reg, a->h_tgt, name);
				if (k == kh_end(a->h_tgt))
					{
					tgt_node = NULL;
					flk_node = NULL;
					goto_next_chromosome = TRUE;
					goto endcore;
					}
				tar = &kh_val(a->h_tgt, k);
				flk = &kh_val(a->h_flk, k);
				if (tar->flag || flk->flag)
					errabort("bam files are not properly sorted\n");
				tgt_node = bed_depnode_list(tar);
				flk_node = bed_depnode_list(flk);
				}
			while (flk_node && flk_node->stop < lstpos+1)
				{
				int j;
				for (j = 0; j < flk_node->len; ++j)
					count_increase(a->c_flkdep, flk_node->vals[j], uint32_t);
				del_node(flk_node);
				if (flk_node && isZero(flk_node->len))
					INIT_DEBUG(depnode_init(flk_node));
				}
			if (flk_node && readcore(flk_node, b)) fs->n_flk++;
			while (tgt_node && tgt_node->stop < lstpos +1)
				{
				print_tgtdep(name, fp, tgt_node, a);
				del_node(tgt_node);
				if (tgt_node && isZero(tgt_node->len))
					INIT_DEBUG(depnode_init(tgt_node));
				}
			if (tgt_node && readcore(tgt_node, b)) fs->n_tgt++;
			endcore:
			bam_destroy1(b);
			}
		}
	khiter_t k, l;
	for (k = 0; k < kh_end(a->h_tgt); ++k)
		{
		if (kh_exist(a->h_tgt, k))
			{
			char *name;
			name = (char *)kh_key(a->h_tgt, k);
			tar = &kh_val(a->h_tgt, k);
			if (tar->flag == 1)
				continue;
			flk = &kh_val(a->h_flk, k);
			tgt_node = bed_depnode_list(tar);
			flk_node = bed_depnode_list(flk);
			while (tgt_node)
				{
				print_zerotgtdep(name, fp, tgt_node, a);
				del_node(tgt_node);
				}
			while (flk_node)
				{
				count_increaseN(a->c_flkdep, 0, flk_node->len, uint32_t);
				del_node(flk_node);
				}
			}
		}
	fclose(fp);
	return 1;
	}
 
struct regcov
	{
	uint64_t cnt, cnt4, cnt10, cnt30, cnt100, cntx;
	float    cov, cov4, cov10, cov30, cov100, covx;
	};

struct regcov *
regcov_init()
	{
	struct regcov *c;
	c = (struct regcov*)needmem(sizeof(struct regcov));
	return c;
	}
 
void cntcov_cal(struct opt_aux *f,
				struct regcov * cov, count32_t * cnt, uint64_t *data)
	{
	uint64_t rawcnt = 0;
	int i;
	*data = 0;
	cov->cnt = cov->cnt4 = cov->cnt10 = 0;
	cov->cnt30 = cov->cnt100 = cov->cntx = 0;
	cov->cov = cov->cov4 = cov->cov10 = 0;
	cov->cov30 = cov->cov100 = cov->covx = 0;
	for (i = 0; i < cnt->m; ++i)
		rawcnt += cnt->a[i];
	for (i = 0; i < cnt->m; ++i)
		{
		(*data) += cnt->a[i] * i;
		if (i < 100)
			{
			cov->cnt100 += cnt->a[i];
			if (i < 30)
				{
				cov->cnt30 += cnt->a[i];
				if (i < 10)
					{
					cov->cnt10 += cnt->a[i];
					if (i < 4) cov->cnt4 += cnt->a[i];
					}
				}
			} //faster
		if (f->cutoff && i < f->cutoff)
			cov->cntx += cnt->a[i];
		}
	cov->cnt = rawcnt - (uint64_t)cnt->a[0];
	cov->cnt4 = rawcnt - cov->cnt4;
	cov->cnt10 = rawcnt - cov->cnt10;
	cov->cnt30 = rawcnt - cov->cnt30;
	cov->cnt100 = rawcnt - cov->cnt100;
	cov->cov = (float)cov->cnt / rawcnt;
	cov->cov4 = (float)cov->cnt4 / rawcnt;
	cov->cov10 = (float)cov->cnt10 / rawcnt;
	cov->cov30 = (float)cov->cnt30 / rawcnt;
	cov->cov100 = (float)cov->cnt100 / rawcnt;
	if (f->cutoff)
		{
		cov->cntx = rawcnt - cov->cntx;
		cov->covx = (float)cov->cntx / rawcnt;
		}
	}

int print_report(struct opt_aux *f, aux_t * a, bamflag_t * fs)
	{
	int i;
	if (f->outdir) chdir(f->outdir);
	FILE *finsert;
	FILE *fdep;
	finsert = open_wfile("insertsize.plot");
	fdep = open_wfile("depth_distribution.plot");
	struct regcov *tarcov = regcov_init();
	struct regcov *flkcov = regcov_init();
	struct regcov *regcov = regcov_init();
	uint64_t icnt = 0;
	uint64_t dcnt = 0;
	for (i = 0; i < a->c_isize->m; ++i)	 icnt += a->c_isize->a[i];
	for (i = 0; i < a->c_isize->m; ++i)
		fprintf(finsert, "%d\t%d\t%f\n",
				i, a->c_isize->a[i], (float)a->c_isize->a[i] / icnt);
	for (i = 0; i < a->c_dep->m; ++i)	 dcnt += a->c_dep->a[i];
	for (i = 0; i < a->c_dep->m; ++i)
		fprintf(fdep, "%d\t%d\t%f\n",
				i, a->c_dep->a[i], (float)a->c_dep->a[i] / dcnt);
	fclose(fdep);
	fclose(finsert);
	cntcov_cal(f, tarcov, a->c_dep, &fs->n_tdata);
	cntcov_cal(f, regcov, a->c_reg, &fs->n_fdata);
	//merge_cnt(a->c_flkdep, a->c_dep);
	cntcov_cal(f, flkcov, a->c_flkdep, &fs->n_fdata);

	FILE *fc = open_wfile("coverage.report");

	{	
	//total
	fprintf(fc, "%60s\t%lld\n", report_total[0], fs->n_reads);
	fprintf(fc, "%60s\t%lld\n", report_total[1], fs->n_qcfail);
	fprintf(fc, "%60s\t%.2f\n", report_total[2], (float)fs->n_data / 1e6);
	fprintf(fc, "%60s\t%lld\n", report_total[3], fs->n_pair_all);
	fprintf(fc, "%60s\t%lld\n", report_total[4], fs->n_mapped);
	fprintf(fc, "%60s\t%.2f\n", report_total[5], (float)fs->n_mapped / fs->n_reads);
	fprintf(fc, "%60s\t%.2f\n", report_total[6], fs->n_mdata / 1e6);
	fprintf(fc, "%60s\t%.2f\n", report_total[7], (float)fs->n_mdata / fs->n_data);
	fprintf(fc, "%60s\t%lld\n", report_total[8], fs->n_pair_good);
	fprintf(fc, "%60s\t%.2f\n", report_total[7], (float)fs->n_pair_good / fs->n_reads);
	fprintf(fc, "%60s\t%lld\n", report_total[9], fs->n_pair_map);
	fprintf(fc, "%60s\t%.2f\n", report_total[10], (float)fs->n_pair_map / fs->n_reads);
	fprintf(fc, "%60s\t%lld\n", report_total[11], fs->n_sgltn);
	fprintf(fc, "%60s\t%lld\n", report_total[13], fs->n_diffchr);
	fprintf(fc, "%60s\t%lld\n", report_total[14], fs->n_read1);
	fprintf(fc, "%60s\t%lld\n", report_total[15], fs->n_read2);
	fprintf(fc, "%60s\t%lld\n", report_total[16], fs->n_pstrand);
	fprintf(fc, "%60s\t%lld\n", report_total[17], fs->n_mstrand);
	fprintf(fc, "%60s\t%lld\n", report_total[18], fs->n_dup);
	fprintf(fc, "%60s\t%.2f\n", report_total[19], (float)fs->n_dup / fs->n_reads);
	fprintf(fc, "%60s\t%d\n", report_total[20], f->mapQ_lim);
	fprintf(fc, "%60s\t%lld\n", report_total[21], fs->n_qual);
	fprintf(fc, "%60s\t%.2f\n", report_total[22], (float)fs->n_qual / fs->n_reads);
	fprintf(fc, "%60s\t%.2f\n", report_total[23], (float)fs->n_qual / fs->n_mapped);
	//tgt
	fprintf(fc, "%60s\t%lld\n", report_tar[0], fs->n_tgt);
	fprintf(fc, "%60s\t%.2f\n", report_tar[1], (float)fs->n_tgt / fs->n_reads);
	fprintf(fc, "%60s\t%.2f\n", report_tar[2], (float)fs->n_tgt / fs->n_mapped);
	fprintf(fc, "%60s\t%.2f\n", report_tar[3], (float)fs->n_tdata / 1e6);
	fprintf(fc, "%60s\t%.2f\n", report_tar[4], (float)fs->n_tdata / fs->n_data);
	fprintf(fc, "%60s\t%.2f\n", report_tar[5], (float)fs->n_tdata / fs->n_mdata);
	fprintf(fc, "%60s\t%lld\n", report_tar[6], a->tgt_len);
	fprintf(fc, "%60s\t%.2f\n", report_tar[7], (float)fs->n_tdata / a->tgt_len);
	fprintf(fc, "%60s\t%.4f\n", report_tar[8], tarcov->cov);
	fprintf(fc, "%60s\t%.4f\n", report_tar[9], tarcov->cov4);
	fprintf(fc, "%60s\t%.4f\n", report_tar[10], tarcov->cov10);
	fprintf(fc, "%60s\t%.4f\n", report_tar[11], tarcov->cov30);
	fprintf(fc, "%60s\t%.4f\n", report_tar[12], tarcov->cov100);
	//tgt regions
	fprintf(fc, "%60s\t%u\n", report_tar[13], a->tgt_nreg);
	fprintf(fc, "%60s\t%lld\n", report_tar[14], regcov->cnt);
	fprintf(fc, "%60s\t%.4f\n", report_tar[15], regcov->cov);
	fprintf(fc, "%60s\t%.4f\n", report_tar[16], regcov->cov4);
	fprintf(fc, "%60s\t%.4f\n", report_tar[17], regcov->cov10);
	fprintf(fc, "%60s\t%.4f\n", report_tar[18], regcov->cov30);
	fprintf(fc, "%60s\t%.4f\n", report_tar[19], regcov->cov100);
	//flk
	fprintf(fc, "%60s\t%u\n", report_tar[20], flank_reg);
	fprintf(fc, "%60s\t%lld\n", report_tar[21], a->flk_len);
	fprintf(fc, "%60s\t%.2f\n", report_tar[22], (float)fs->n_fdata / a->flk_len);
	fprintf(fc, "%60s\t%lld\n", report_tar[23], fs->n_flk);
	fprintf(fc, "%60s\t%.2f\n", report_tar[24], (float)fs->n_flk / fs->n_reads);
	fprintf(fc, "%60s\t%.2f\n", report_tar[25], (float)fs->n_flk / fs->n_mapped);
	fprintf(fc, "%60s\t%.2f\n", report_tar[26], (float)fs->n_fdata / 1e6);
	fprintf(fc, "%60s\t%.2f\n", report_tar[27], (float)fs->n_fdata / fs->n_data);
	fprintf(fc, "%60s\t%.2f\n", report_tar[28], (float)fs->n_fdata / fs->n_mdata);
	fprintf(fc, "%60s\t%.4f\n", report_tar[29], flkcov->cov);
	fprintf(fc, "%60s\t%.4f\n", report_tar[30], flkcov->cov4);
	fprintf(fc, "%60s\t%.4f\n", report_tar[31], flkcov->cov10);
	fprintf(fc, "%60s\t%.4f\n", report_tar[32], flkcov->cov30);
	fprintf(fc, "%60s\t%.4f\n", report_tar[33], flkcov->cov100);
	}
	mustfree(tarcov);
	mustfree(regcov);
	mustfree(flkcov);
	fclose(fc);
	return 1;
	}

enum
	{
	MAXDEPTH,
	CUTOFF,
	INSERTSIZE,
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
	{"rmdup", no_argument, NULL, 'd'},
	{"help", no_argument, NULL, 'h'}
	};

int bamdst(int argc, char *argv[])
	{
	int n, i;
	char *probe = 0;
	//bool help = FALSE;
	struct opt_aux *opt = opt_init();
	while ((n = getopt_long(argc, argv, "o:p:f:q:l:dh", long_opts, NULL)) >= 0)
		{
		switch (n)
			{
			//output dir, must have right to write
			case 'o': opt->outdir = strdup(optarg); break;
				//capture region or just the region you interesting
			case 'p': probe = optarg; break;
				//flk the region for more information, default is 200 bp
			case 'f': flank_reg = atoi(optarg); break;
				//max depth to considered in the cumulative distribution of depths
			case MAXDEPTH: opt->maxdepth = atoi(optarg); break;
			case CUTOFF: opt->cutoff = atoi(optarg); break;
			case INSERTSIZE: opt->isize_lim = atoi(optarg); break;
			case 'q': opt->mapQ_lim = atoi(optarg); break;
			case 'h': usage(1); break;
			case 'd': rmdup_mark = TRUE; break;
			default: usage(0); 
				//more help
			}
		
		}
	if (isNull(opt->outdir) || isNull(probe))
		errabort("out dir and bed file are mandatory!");
	n = argc - optind;
	//capable of deals with severl bam files
	aux_t * aux;
	aux = aux_init();
	if (isZero(n))
		{
		aux->data = (bamFile*)needmem(sizeof(bamFile));
		aux->data[0] = bgzf_dopen(fileno(stdin), "r");
		aux->h = bam_header_read(aux->data[0]);
		aux->ndata = 1;
		opt->nfiles = 0;
		}
	else
		{
		aux->data = (bamFile*)needmem(n * sizeof(bamFile));
		opt->nfiles = n;
		opt->inputs = (char**)needmem(n * sizeof(char*));
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
				errabort("%s: %s\n", argv[optind + i], strerror(errno));
			h_tmp = bam_header_read(aux->data[i]);
			if (i == 0)	aux->h = h_tmp;
			else bam_header_destroy(h_tmp);
			opt->inputs[i] = strdup(argv[optind+i]);
			}
		aux->ndata = n;
		}
	load_bed_init(probe, aux);
	aux->c_isize->a = calloc(opt->isize_lim, sizeof(unsigned));
	for (i = 0; i < opt->isize_lim; ++i) aux->c_isize->a[i] = 0;
	aux->c_isize->n = opt->isize_lim;
	aux->nchr = aux->h->n_targets;
	bamflag_t *fs = bamflag_init();
	load_bamfiles(opt, aux, fs);
	print_report(opt, aux, fs);
	mustfree(fs);
	aux_destroy(aux);
	opt_destroy(opt);
	return 1;
	}

/* main */
int main(int argc, char *argv[])
	{
	bamdst(argc, argv);
	return 1;
	}
