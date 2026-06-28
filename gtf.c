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

#include "commons.h"
#include "khash.h"
#include "kseq.h"
#include "kstring.h"
#include "ksort.h"
#include "dict.h"
#include "gtf.h"
#include "bedutil.h"
#include "number.h"
#include <zlib.h>
#include <sys/time.h>
#include <time.h>

KSTREAM_INIT(gzFile, gzread, 8193)

KHASH_MAP_INIT_INT(attr, char*)

/* ---- local helpers (adapted from PISA utils.h) ---- */

static inline double realtime(void)
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

#define LOG_print(line, ...) do {                                       \
    time_t _t;                                                          \
    time(&_t);                                                          \
    struct tm *_tm = localtime(&_t);                                    \
    char _ts[32];                                                       \
    strftime(_ts, 32, "%Y-%m-%d %H:%M:%S", _tm);                       \
    fprintf(stderr, "[%s] " line "\n", _ts, ##__VA_ARGS__);            \
} while(0)

/* ---- feature type names ---- */

static const char *feature_type_names[] = {
    "gene",
    "transcript",
    "CDS",
    "start_codon",
    "stop_codon",
    "5UTR",
    "3UTR",
    "inter",
    "inter_CNS",
    "intron_CNS",
    "exon",
    "five_prime_utr",
    "three_prime_utr",
    "Selenocysteine"
};

const char *get_feature_name(enum feature_type type)
{
    assert(type > -1);
    return feature_type_names[type];
}

/* ---- gtf struct helpers ---- */

void gtf_reset(struct gtf *gtf)
{
    memset(gtf, 0, sizeof(struct gtf));
    gtf->seqname = gtf->source = gtf->start = gtf->end =
        gtf->gene_id = gtf->gene_name = gtf->transcript_id = -1;
}

struct gtf *gtf_create(void)
{
    struct gtf *g = malloc(sizeof(*g));
    gtf_reset(g);
    return g;
}

void gtf_attr_clear(struct attr *attr)
{
    if (attr) {
        if (attr->next) gtf_attr_clear(attr->next);
        if (attr->val) free(attr->val);
        free(attr);
    }
}

void gtf_clear(struct gtf *gtf)
{
    int i;
    for (i = 0; i < gtf->n_gtf; ++i) {
        gtf_clear(gtf->gtf[i]);
        free(gtf->gtf[i]);
    }
    if (gtf->n_gtf) free(gtf->gtf);
    if (gtf->attr != NULL) gtf_attr_clear(gtf->attr);
}

void gtf_copy(struct gtf *dest, struct gtf *src)
{
    gtf_reset(dest);
    dest->seqname = src->seqname;
    dest->source = src->source;
    dest->type = src->type;
    dest->start = src->start;
    dest->end = src->end;
    dest->strand = src->strand;
    dest->gene_id = src->gene_id;
    dest->gene_name = src->gene_name;
    dest->transcript_id = src->transcript_id;
    dest->attr = src->attr;
    src->attr = NULL;
}

/* ---- sorting ---- */

static int cmpfunc1(const void *_a, const void *_b)
{
    const struct gtf *a = *(const struct gtf**)_a;
    const struct gtf *b = *(const struct gtf**)_b;
    if (a->seqname != b->seqname) return (a->seqname > b->seqname) - (a->seqname < b->seqname);
    if (a->start != b->start) return (a->start > b->start) - (a->start < b->start);
    return (a->end > b->end) - (a->end < b->end);
}

/* ---- GFF attribute parser ---- */

struct attr_pair {
    char *key;
    char *val;
};

static struct attr_pair *split_gff(kstring_t *str, int *_n)
{
    int i = 0;
    int n = 0, m = 0;
    struct attr_pair *pair = NULL;

    int j = str->l - 1;
    if (j < 0) { *_n = 0; return NULL; }
    while (j >= 0 && (isspace(str->s[j]) || str->s[j] == ';')) j--;
    if (j < 0) { *_n = 0; return NULL; }
    str->l = j + 1;
    str->s[str->l] = '\0';

    for (;;) {
        if (i >= str->l) break;
        if (n == m) {
            m += 4;
            pair = realloc(pair, sizeof(struct attr_pair)*m);
        }

        kstring_t name = {0,0,0};
        kstring_t val = {0,0,0};
        while (i < str->l && !isspace(str->s[i]) && str->s[i] != ';') {
            kputc(str->s[i], &name);
            ++i;
        }

        while (isspace(str->s[i]) || str->s[i] == ';') ++i;

        if (str->s[i] == '"') {
            ++i;
            for (; i < str->l;) {
                if (str->s[i] == '"' || i+1 == str->l) {
                    i++;
                    i++;
                    break;
                }
                kputc(str->s[i], &val);
                i++;
            }
        }

        while (i < str->l && (isspace(str->s[i]) || str->s[i] == ';')) ++i;
        if (name.l == 0) {
            warnings("Empty key. %s", str->s);
            continue;
        }

        pair[n].key = name.s;
        pair[n].val = val.s;
        n++;
    }
    *_n = n;
    return pair;
}

static struct attr_pair *bend_pair(char *s, int *n)
{
    if (s == NULL) return NULL;
    kstring_t str = {0,0,0};
    kputs(s, &str);

    struct attr_pair *p = split_gff(&str, n);
    free(str.s);
    return p;
}

/* ---- push a GTF record into the gene→transcript→exon tree ---- */

static int gtf_push(struct gtf_spec *G, struct gtf_ctg *ctg,
                     struct gtf *gtf, int feature)
{
    char *gene_id = dict_name(G->gene_id, gtf->gene_id);
    int gene_idx = dict_push(ctg->gene_idx, gene_id);

    struct gtf *gene_gtf = dict_query_value(ctg->gene_idx, gene_idx);

    if (feature == feature_gene && gene_gtf != NULL) {
        warnings("Duplicated gene record? %s",
                 dict_name(G->gene_name, gtf->gene_name));
        if (gene_gtf->start < 0 || gene_gtf->start > gtf->start)
            gene_gtf->start = gtf->start;
        if (gene_gtf->end < 0 || gene_gtf->end < gtf->end)
            gene_gtf->end = gtf->end;
        gene_gtf->strand = gtf->strand;
        gene_gtf->gene_id = gtf->gene_id;
        gene_gtf->gene_name = gtf->gene_name;
        gene_gtf->seqname = gtf->seqname;
        return 1;
    }

    if (gene_gtf == NULL) {
        if (ctg->n_gtf == ctg->m_gtf) {
            ctg->m_gtf = ctg->m_gtf == 0 ? 4: ctg->m_gtf*2;
            ctg->gtf = realloc(ctg->gtf, sizeof(struct gtf*)*ctg->m_gtf);
        }
        ctg->gtf[ctg->n_gtf] = gtf_create();
        gene_gtf = ctg->gtf[ctg->n_gtf++];

        dict_assign_value(ctg->gene_idx, gene_idx, gene_gtf);

        gtf_reset(gene_gtf);
        if (feature == feature_gene) {
            gtf_copy(gene_gtf, gtf);
            return 0;
        }
        gene_gtf->type = feature_gene;
        gene_gtf->strand = gtf->strand;
        gene_gtf->seqname = gtf->seqname;
    }

    if (gene_gtf->start < 0 || gene_gtf->start > gtf->start)
        gene_gtf->start = gtf->start;
    if (gene_gtf->end < 0 || gene_gtf->end < gtf->end)
        gene_gtf->end = gtf->end;

    if (gtf->transcript_id == -1)
        errabort("No transcript found. %s:%s:%d:%d",
                 feature_type_names[feature],
                 dict_name(G->name, gtf->seqname), gtf->start, gtf->end);

    if (gene_gtf->gene_id == -1) gene_gtf->gene_id = gtf->gene_id;
    if (gene_gtf->gene_name == -1) gene_gtf->gene_name = gtf->gene_name;

    struct gtf *tx_gtf = NULL;
    int i;
    for (i = 0; i < gene_gtf->n_gtf; ++i) {
        struct gtf *tx_gtf0 = gene_gtf->gtf[i];
        if (gtf->transcript_id == tx_gtf0->transcript_id) {
            tx_gtf = tx_gtf0;
            break;
        }
    }
    if (feature == feature_transcript && tx_gtf != NULL) {
        warnings("Duplicated transcript record? %s",
                 dict_name(G->transcript_id, gtf->transcript_id));
        if (tx_gtf->start < 0 || tx_gtf->start > gtf->start)
            tx_gtf->start = gtf->start;
        if (tx_gtf->end < 0 || tx_gtf->end < gtf->end)
            tx_gtf->end = gtf->end;
        tx_gtf->strand = gene_gtf->strand;
        tx_gtf->gene_id = gene_gtf->gene_id;
        tx_gtf->gene_name = gene_gtf->gene_name;
        tx_gtf->seqname = gene_gtf->seqname;
        return 1;
    }

    if (tx_gtf == NULL) {
        if (gene_gtf->n_gtf == gene_gtf->m_gtf) {
            gene_gtf->m_gtf = gene_gtf->m_gtf == 0? 4 : gene_gtf->m_gtf*2;
            gene_gtf->gtf = realloc(gene_gtf->gtf,
                                     gene_gtf->m_gtf *sizeof(struct gtf*));
        }
        gene_gtf->gtf[gene_gtf->n_gtf] = gtf_create();
        tx_gtf = gene_gtf->gtf[gene_gtf->n_gtf++];

        gtf_reset(tx_gtf);

        if (feature == feature_transcript) {
            gtf_copy(tx_gtf, gtf);
            if (tx_gtf->strand == -1) tx_gtf->strand = gene_gtf->strand;
            if (tx_gtf->gene_name == -1) tx_gtf->gene_name = gene_gtf->gene_name;
            if (tx_gtf->gene_id == -1) tx_gtf->gene_id = gene_gtf->gene_id;
            return 0;
        }

        tx_gtf->type = feature_transcript;
        tx_gtf->strand = gene_gtf->strand;
        tx_gtf->seqname = gene_gtf->seqname;
    }

    if (tx_gtf->gene_id == -1) tx_gtf->gene_id = gene_gtf->gene_id;
    if (tx_gtf->gene_name == -1) tx_gtf->gene_name = gene_gtf->gene_name;
    if (tx_gtf->transcript_id == -1) tx_gtf->transcript_id = gtf->transcript_id;

    if (tx_gtf->start < 0 || tx_gtf->start > gtf->start)
        tx_gtf->start = gtf->start;
    if (tx_gtf->end < 0 || tx_gtf->end < gtf->end)
        tx_gtf->end = gtf->end;

    if (gtf->gene_id >= 0 && tx_gtf->gene_id != gtf->gene_id) {
        errabort("Inconsitance gene id between transcript %s and gene %s.",
                 dict_name(G->transcript_id, tx_gtf->transcript_id),
                 dict_name(G->gene_id, gtf->gene_id));
    }
    if (gtf->gene_name >= 0 && tx_gtf->gene_name != gtf->gene_name) {
        errabort("Inconsitance gene name between transcript %s and gene %s.",
                 dict_name(G->transcript_id, tx_gtf->transcript_id),
                 dict_name(G->gene_name, gtf->gene_name));
    }
    if (gtf->transcript_id >= 0 &&
        tx_gtf->transcript_id != gtf->transcript_id) {
        errabort("Inconsitance transcript id between transcript %s and record %s.",
                 dict_name(G->transcript_id, tx_gtf->transcript_id),
                 dict_name(G->transcript_id, gtf->transcript_id));
    }

    if (tx_gtf->n_gtf == tx_gtf->m_gtf) {
        tx_gtf->m_gtf = tx_gtf->m_gtf == 0 ? 4 : tx_gtf->m_gtf*2;
        tx_gtf->gtf = realloc(tx_gtf->gtf, tx_gtf->m_gtf*sizeof(struct gtf*));
    }
    tx_gtf->gtf[tx_gtf->n_gtf] = gtf_create();
    struct gtf *exon_gtf = tx_gtf->gtf[tx_gtf->n_gtf++];
    gtf_copy(exon_gtf, gtf);

    if (exon_gtf->type == feature_CDS) {
        tx_gtf->coding = 1;
        gene_gtf->coding = 1;
    }

    if (exon_gtf->gene_id == -1) exon_gtf->gene_id = gene_gtf->gene_id;
    if (exon_gtf->gene_name == -1) exon_gtf->gene_name = gene_gtf->gene_name;
    if (exon_gtf->transcript_id == -1) exon_gtf->transcript_id = tx_gtf->transcript_id;

    return 0;
}

/* ---- GTF line parser ---- */

#define FILTER_ATTRS  2
#define FILTER_TRANS  1

static int parse_str(struct gtf_spec *G, kstring_t *str, int filter)
{
    int n;
    int *s = ksplit(str, '\t', &n);
    if (n != 9) errabort("Unknown format. %s", str->s);

    char *feature = str->s + s[2];

    int qry = dict_query(G->features, feature);
    if (qry == -1) {
        free(s);
        return 1;
    }

    if ((filter & 0x3) & FILTER_TRANS) {
        if (qry != feature_gene && qry != feature_exon &&
            qry != feature_transcript) {
            free(s);
            return 0;
        }
    }

    struct gtf gtf;
    gtf_reset(&gtf);
    gtf.seqname = dict_push(G->name, str->s + s[0]);
    if (s[1] != 1 || str->s[s[1]] != '.')
        gtf.source = dict_push(G->sources, str->s + s[1]);
    gtf.type = qry;
    gtf.start = str2int(str->s + s[3]);
    gtf.end = str2int(str->s + s[4]);
    char *strand = str->s + s[6];
    gtf.strand = strand[0] == '-' ? 1 : 0;
    char *attr = str->s + s[8];

    struct gtf_ctg *ctg = dict_query_value(G->name, gtf.seqname);
    if (ctg == NULL) {
        ctg = malloc(sizeof(struct gtf_ctg));
        memset(ctg, 0, sizeof(struct gtf_ctg));
        ctg->gene_idx = dict_init();
        dict_set_value(ctg->gene_idx);
        dict_assign_value(G->name, gtf.seqname, ctg);
    }

    int i;
    int n0 = 0;
    struct attr_pair *pair = bend_pair(attr, &n0);
    for (i = 0; i < n0; ++i) {
        struct attr_pair *pp = &pair[i];
        if (strcmp(pp->key, "gene_id") == 0)
            gtf.gene_id = dict_push(G->gene_id, pp->val);
        else if (strcmp(pp->key, "gene_name") == 0)
            gtf.gene_name = dict_push(G->gene_name, pp->val);
        else if (strcmp(pp->key, "gene") == 0)
            gtf.gene_name = dict_push(G->gene_name, pp->val);
        else if (strcmp(pp->key, "transcript_id") == 0)
            gtf.transcript_id = dict_push(G->transcript_id, pp->val);
        else {
            if ((filter & 0x3) & FILTER_ATTRS) {
                struct attr *attr_new = malloc(sizeof(struct attr));
                attr_new->id = dict_push(G->attrs, pp->key);
                attr_new->val = NULL;
                attr_new->next = NULL;
                if (gtf.attr == NULL) gtf.attr = attr_new;
                else {
                    struct attr *tmp;
                    for (tmp = gtf.attr; tmp->next; tmp = tmp->next);
                    tmp->next = attr_new;
                }
                if (pp->val != NULL) attr_new->val = strdup(pp->val);
            }
        }
        free(pp->key);
        if (pp->val) free(pp->val);
    }

    free(pair);
    free(s);

    if (gtf.gene_id == -1 && gtf.gene_name == -1) {
        warnings("Record %s:%s:%d-%d has no gene_name and gene_id. Skip.",
                 dict_name(G->name, gtf.seqname),
                 feature_type_names[qry], gtf.start, gtf.end);
        gtf_clear(&gtf);
        return 1;
    }
    if (gtf_push(G, ctg, &gtf, qry)) {
        warnings("Failed to push record, %s:%s:%d-%d",
                 dict_name(G->name, gtf.seqname),
                 feature_type_names[qry], gtf.start, gtf.end);
    }
    gtf_clear(&gtf);
    return 0;
}

/* ---- recursive sort of the gene→transcript→exon tree ---- */

static void gtf_sort(struct gtf *gtf)
{
    int i;
    for (i = 0; i < gtf->n_gtf; ++i)
        gtf_sort(gtf->gtf[i]);

    if (gtf->n_gtf) {
        qsort((const struct gtf**)gtf->gtf, gtf->n_gtf, sizeof(struct gtf*),
              cmpfunc1);
        assert(gtf->start < gtf->end);
    }
}

/* ---- build index (gene lookup dicts — no region_index) ---- */

static int gtf_build_index(struct gtf_spec *G)
{
    int i;
    int total_gene = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        assert(ctg);
        qsort((const struct gtf**)ctg->gtf, ctg->n_gtf, sizeof(struct gtf*),
              cmpfunc1);
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            gtf_sort(ctg->gtf[j]);
            struct gtf *g = ctg->gtf[j];
            struct gtf *g0 = dict_query_value(G->gene_name, g->gene_name);
            if (g0 == NULL) {
                dict_assign_value(G->gene_name, g->gene_name, g);
            } else {
                for (;;) {
                    if (g0->ext == NULL) {
                        g0->ext = g;
                        break;
                    }
                    g0 = g0->ext;
                }
            }
            int k;
            for (k = 0; k < g->n_gtf; ++k) {
                struct gtf *tx = g->gtf[k];
                struct gtf *tx0 = dict_query_value(G->transcript_id,
                                                    tx->transcript_id);
                if (tx0 == NULL) {
                    dict_assign_value(G->transcript_id, tx->transcript_id, tx);
                } else {
                    for (;;) {
                        if (tx0->ext == NULL) {
                            tx0->ext = tx;
                            break;
                        }
                        tx0 = tx0->ext;
                    }
                }
            }
        }
        total_gene += ctg->n_gtf;
    }
    return total_gene;
}

/* ---- init / read / destroy ---- */

struct gtf_spec *gtf_spec_init(void)
{
    struct gtf_spec *G = malloc(sizeof(*G));
    memset(G, 0, sizeof(*G));
    G->name          = dict_init();
    G->gene_name     = dict_init();
    G->gene_id       = dict_init();
    G->transcript_id = dict_init();
    G->sources       = dict_init();
    G->attrs         = dict_init();
    G->features      = dict_init();

    dict_set_value(G->name);
    dict_set_value(G->gene_name);
    dict_set_value(G->transcript_id);
    int i;
    int l = sizeof(feature_type_names)/sizeof(feature_type_names[0]);
    for (i = 0; i < l; ++i)
        dict_push(G->features, (char*)feature_type_names[i]);

    return G;
}

struct gtf_spec *gtf_read(const char *fname, int f)
{
    LOG_print("GTF loading..");
    double t_real;
    t_real = realtime();

    gzFile fp;
    fp = gzopen(fname, "r");
    if (fp == NULL) errabort("%s : %s.", fname, strerror(errno));

    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    int line = 0;
    struct gtf_spec *G = gtf_spec_init();

    while (ks_getuntil(ks, 2, &str, &ret) >= 0) {
        line++;
        if (str.l == 0) {
            warnings("Line %d is empty. Skip.", line);
            continue;
        }
        if (str.s[0] == '#') continue;
        parse_str(G, &str, f);
    }
    free(str.s);
    gzclose(fp);
    ks_destroy(ks);

    if (dict_size(G->name) == 0) {
        gtf_destroy(G);
        return NULL;
    }

    int n_gene = gtf_build_index(G);
    LOG_print("Load %d genes.", n_gene);
    LOG_print("Load time : %.3f sec", realtime() - t_real);
    return G;
}

struct gtf_spec *gtf_read_lite(const char *fname)
{
    return gtf_read(fname, 1);
}

void gtf_destroy(struct gtf_spec *G)
{
    int i;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        if (ctg->gene_idx) dict_destroy(ctg->gene_idx);

        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            gtf_clear(ctg->gtf[j]);
            free(ctg->gtf[j]);
        }

        free(ctg->gtf);
        free(ctg);
    }
    dict_destroy(G->name);
    dict_destroy(G->gene_name);
    dict_destroy(G->gene_id);
    dict_destroy(G->transcript_id);
    dict_destroy(G->sources);
    dict_destroy(G->attrs);
    dict_destroy(G->features);
    free(G);
}

char *GTF_seqname(struct gtf_spec *G, int id)
{
    return dict_name(G->name, id);
}

char *GTF_genename(struct gtf_spec *G, int id)
{
    return dict_name(G->gene_name, id);
}

char *GTF_transid(struct gtf_spec *G, int id)
{
    return dict_name(G->transcript_id, id);
}

/* ---- dump GTF to file ---- */

void write_gtf_fp(struct gtf_spec *G, struct gtf *gtf, FILE *fp,
                   struct dict *keys)
{
    fputs(dict_name(G->name, gtf->seqname), fp);
    if (gtf->source == -1) fputs("\t.\t", fp);
    else fprintf(fp, "\t%s\t", dict_name(G->sources, gtf->source));
    fputs(feature_type_names[gtf->type], fp);
    fprintf(fp, "\t%d\t%d\t.\t%c\t.\t", gtf->start, gtf->end,
            "+-"[gtf->strand]);
    fprintf(fp, "gene_name \"%s\"; gene_id \"%s\";",
            dict_name(G->gene_name, gtf->gene_name),
            dict_name(G->gene_id, gtf->gene_id));

    if (gtf->transcript_id >= 0) {
        fprintf(fp, " transcript_id \"%s\";",
                dict_name(G->transcript_id, gtf->transcript_id));
    }

    if (gtf->attr) {
        struct attr *tmp = gtf->attr;
        for (tmp = gtf->attr; tmp; tmp = tmp->next) {
            char *name = dict_name(G->attrs, tmp->id);
            fprintf(fp, " %s", name);
            if (tmp->val) fprintf(fp, " \"%s\";", tmp->val);
            else fputc(';', fp);
        }
    }

    fputc('\n', fp);

    int i;
    for (i = 0; i < gtf->n_gtf; ++i)
        write_gtf_fp(G, gtf->gtf[i], fp, keys);
}

void gtf_dump(struct gtf_spec *G, const char *fname, struct dict *keys)
{
    FILE *fp = fname == NULL ? stdout : fopen(fname, "w");
    if (fp == NULL) errabort("%s : %s.", fname, strerror(errno));
    int i, j;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gtf = ctg->gtf[j];
            write_gtf_fp(G, gtf, fp, keys);
        }
    }
    if (fp != stdout) fclose(fp);
}

/* ---- gene / transcript lookup ---- */

struct gtf *gtf_query_gene(struct gtf_spec *G, const char *name)
{
    return (struct gtf*)dict_query_value2(G->gene_name, name);
}

struct gtf *gtf_query_tx(struct gtf_spec *G, const char *name)
{
    return (struct gtf*)dict_query_value2(G->transcript_id, name);
}

/* ================================================================
 * gtf_to_regHash — convert parsed GTF tree to bamdst regHash_t
 *
 * Walks the gene→transcript→exon tree and pushes intervals into
 * a bedutil regHash_t.  This is the bridge between GTF parsing
 * and bamdst's existing coverage pipeline.
 * ================================================================ */

/* external bedHandle from bedutil.c */
extern bedHandle_t *bedHand;

static void gtf_push_interval(regHash_t *h, const char *chr,
                               uint32_t beg, uint32_t end)
{
    khiter_t k;
    int ret;
    k = kh_put(reg, h, strdup(chr), &ret);
    bedreglist_t *b;
    if (ret) {
        /* new chromosome entry */
        b = (bedreglist_t*)needmem(sizeof(bedreglist_t));
        memset(b, 0, sizeof(bedreglist_t));
        kh_val(h, k) = *b;
    }
    b = &kh_val(h, k);

    if (b->n == 0) {
        b->n = 2;
        b->a = (uint64_t *)needmem(b->n * sizeof(uint64_t));
    } else if (b->m == b->n) {
        b->n = b->n + 1024;
        b->a = (uint64_t *)enlarge_empty_mem((void *)b->a,
                     b->m * sizeof(uint64_t), b->n * sizeof(uint64_t));
    }
    /* pack: high 32 bits = begin, low 32 bits = end (0-based) */
    b->a[b->m++] = (uint64_t)beg << 32 | (uint32_t)end;
}

regHash_t *gtf_to_regHash(struct gtf_spec *G, int level)
{
    regHash_t *h = kh_init(reg);
    int i;

    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        char *seqname = dict_name(G->name, i);
        int j;

        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gene = ctg->gtf[j];

            if (level == GTF_LEVEL_GENE) {
                /* push gene span (0-based start, 1-based end) */
                gtf_push_interval(h, seqname,
                                  (uint32_t)(gene->start - 1),
                                  (uint32_t)gene->end);
            } else if (level == GTF_LEVEL_EXON) {
                int k;
                for (k = 0; k < gene->n_gtf; ++k) {
                    struct gtf *tx = gene->gtf[k];
                    int e;
                    for (e = 0; e < tx->n_gtf; ++e) {
                        struct gtf *ex = tx->gtf[e];
                        if (ex->type == feature_exon) {
                            gtf_push_interval(h, seqname,
                                              (uint32_t)(ex->start - 1),
                                              (uint32_t)ex->end);
                        }
                    }
                }
            }
        }
    }

    return h;
}

/*
 * Strand-aware GTF → regHash_t conversion.
 * Splits intervals by strand: forward ('+' / unstranded) and reverse ('-').
 */
void gtf_to_regHash_stranded(struct gtf_spec *G, int level,
                              regHash_t **h_fwd, regHash_t **h_rev)
{
    *h_fwd = kh_init(reg);
    *h_rev = kh_init(reg);
    int i;

    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        char *seqname = dict_name(G->name, i);
        int j;

        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gene = ctg->gtf[j];
            /* determine target hash by strand:
             * GTF_STRAND_FWD (0) → h_fwd, GTF_STRAND_REV (1) → h_rev */
            regHash_t *target = (gene->strand == GTF_STRAND_REV)
                              ? *h_rev : *h_fwd;

            if (level == GTF_LEVEL_GENE) {
                gtf_push_interval(target, seqname,
                                  (uint32_t)(gene->start - 1),
                                  (uint32_t)gene->end);
            } else if (level == GTF_LEVEL_EXON) {
                int k;
                for (k = 0; k < gene->n_gtf; ++k) {
                    struct gtf *tx = gene->gtf[k];
                    int e;
                    for (e = 0; e < tx->n_gtf; ++e) {
                        struct gtf *ex = tx->gtf[e];
                        if (ex->type == feature_exon) {
                            /* exon strand follows the gene strand */
                            regHash_t *et = (gene->strand == GTF_STRAND_REV)
                                          ? *h_rev : *h_fwd;
                            gtf_push_interval(et, seqname,
                                              (uint32_t)(ex->start - 1),
                                              (uint32_t)ex->end);
                        }
                    }
                }
            }
        }
    }
}
