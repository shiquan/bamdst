/*
   Copyright (c) 2016-2026 The bamdst Authors

   Adapted from PISA (https://github.com/shiquan/PISA)

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

#ifndef GTF_H
#define GTF_H

#include <stdlib.h>
#include "commons.h"
#include "dict.h"

#define EXONIC    0
#define INTRONIC  1
#define GTF_STRAND_FWD 0
#define GTF_STRAND_REV 1
#define GTF_STRAND_UNK -1

/* GTF coverage level (for gtf_to_regHash) */
#define GTF_LEVEL_GENE   1
#define GTF_LEVEL_EXON   3

enum feature_type {
    feature_unknow = -1,
    feature_gene,
    feature_transcript,
    feature_CDS,
    feature_start_codon,
    feature_stop_codon,
    feature_5UTR,
    feature_3UTR,
    feature_inter,
    feature_inter_CNS,
    feature_intron_CNS,
    feature_exon,
    feature_5UTR_alias,
    feature_3UTR_alias,
    feature_Selenocysteine,
    feature_antisense,
};

struct attr {
    int id;
    char *val;
    struct attr *next;
};

struct gtf {
    int seqname;
    int source;
    enum feature_type type;
    int start;
    int end;
    int strand;     /* 0 forward, 1 reverse */
    int gene_id;
    int gene_name;
    int transcript_id;
    int coding;     /* set to 1 if CDS exists */
    struct attr *attr;
    struct gtf *ext; /* next record (for same gene_name lookup) */
    int n_gtf, m_gtf;
    struct gtf **gtf;
};

struct gtf_ctg {
    struct dict *gene_idx;
    int n_gtf, m_gtf;
    struct gtf **gtf;
};

struct gtf_spec {
    struct dict *name;     /* contig/chromosome names */
    struct dict *gene_name;
    struct dict *gene_id;
    struct dict *transcript_id;
    struct dict *sources;
    struct dict *attrs;
    struct dict *features;
};

const char *get_feature_name(enum feature_type type);
char *GTF_seqname(struct gtf_spec *G, int id);
char *GTF_genename(struct gtf_spec *G, int id);
char *GTF_transid(struct gtf_spec *G, int id);

struct gtf_spec *gtf_read(const char *fname, int filter);
struct gtf_spec *gtf_read_lite(const char *fname);
void gtf_destroy(struct gtf_spec *G);
void gtf_dump(struct gtf_spec *G, const char *fname, struct dict *);
struct gtf *gtf_query_gene(struct gtf_spec *G, const char *name);
struct gtf *gtf_query_tx(struct gtf_spec *G, const char *name);

/* Convert parsed GTF to bamdst's region hash (bedutil regHash_t).
 *  level: GTF_LEVEL_GENE or GTF_LEVEL_EXON
 *  Returns a regHash_t mapping chromosome name → merged intervals.
 *  Caller must destroy with bedHand->destroy(). */
#include "bedutil.h"
regHash_t *gtf_to_regHash(struct gtf_spec *G, int level);

/* Strand-aware conversion: separate forward- and reverse-strand intervals.
 *  h_fwd: forward-strand genes (GTF_STRAND_FWD / '+')
 *  h_rev: reverse-strand genes (GTF_STRAND_REV / '-')
 *  Intervals with unknown strand are omitted.
 *  Caller must destroy both hashes with bedHand->destroy(). */
void gtf_to_regHash_stranded(struct gtf_spec *G, int level,
                              regHash_t **h_fwd, regHash_t **h_rev);

#endif
