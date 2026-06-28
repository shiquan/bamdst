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

#include "dict.h"
#include "khash.h"
#include "kseq.h"
#include "kstring.h"
#include <zlib.h>

typedef struct dict_val {
    int idx;
    int ref;
} dict_val_t;

KHASH_MAP_INIT_STR(name, dict_val_t)

KSTREAM_INIT(gzFile, gzread, 8193)

struct dict {
    int n, m;
    char **name;
    kh_name_t *dict;
    uint32_t *count;
    int assign_value_flag;
    void **value;
};

struct dict *dict_init()
{
    struct dict *D = malloc(sizeof(*D));
    memset(D, 0, sizeof(*D));
    D->dict = kh_init(name);
    return D;
}
void dict_set_value(struct dict *D)
{
    if (D->assign_value_flag == 1) errabort("Double assign value to a dict.");
    D->assign_value_flag = 1;
    if (D->m > 0) {
        D->value = malloc(D->m*sizeof(void*));
        int i;
        for (i = 0; i < D->m; ++i) D->value[i] = NULL;
    }
}
void *dict_query_value(struct dict *D, int idx)
{
    if (idx < 0 || idx >= D->n) return NULL;
    if (D->value == NULL) return NULL;

    return D->value[idx];
}

void *dict_query_value2(struct dict *D, const char *key)
{
    int idx = dict_query(D, key);
    if (idx == -1) return NULL; // failed to query
    return D->value[idx];
}


int dict_assign_value(struct dict *D, int idx, void *val)
{
    if (idx < 0 || idx >= D->n) return 1;
    D->value[idx] = val;
    return 0;
}

int dict_delete_value(struct dict *D, int idx)
{
    if (idx < 0 || idx >= D->n) return 1;
    void *p = D->value[idx];
    D->value[idx] = NULL;
    free(p);
    return 0;
}

char *dict_name(const struct dict *D, int idx)
{
    assert(idx >= 0 && idx < D->n);
    return D->name[idx];
}

int dict_size(const struct dict *D)
{
    return D->n;
}
uint32_t dict_count(const struct dict *D, int idx)
{
    return D->count[idx];
}
uint32_t dict_count_sum(const struct dict *D)
{
    uint32_t sum = 0;
    int i;
    for (i = 0; i < D->n; ++i) sum+=D->count[i];
    return sum;
}
void dict_destroy(struct dict *D)
{
    if (D == NULL) return;
    int i;
    for (i = 0; i < D->n; ++i) free(D->name[i]);
    // values are actually points, need free pointed values manually
    if (D->assign_value_flag && D->value) free(D->value);
    if (D->n > 0) {
        free(D->name);
        free(D->count);
    }
    kh_destroy(name,D->dict);
    free(D);
}


int dict_query(const struct dict *D, char const *key)
{
    if (key == NULL) errabort("Trying to query an empty key.");
    khint_t k;
    k = kh_get(name, D->dict, key);
    if (k == kh_end(D->dict)) return -1;
    return kh_val(D->dict, k).idx;
}

struct dict_val *dict_query0(const struct dict *D, char const *key)
{
    if (key == NULL) errabort("Trying to query an empty key.");
    khint_t k;
    k = kh_get(name, D->dict, key);
    if (k == kh_end(D->dict)) return NULL;
    return &kh_val(D->dict, k);
}

int dict_query2(const struct dict *D, char const *key)
{
    if (key == NULL) errabort("Trying to query an empty key.");
    khint_t k;
    k = kh_get(name, D->dict, key);
    if (k == kh_end(D->dict)) return -1;
    return kh_val(D->dict, k).ref;
}

static int dict_push0(struct dict *D, char const *key, int idx)
{
    if (key == NULL) errabort("Trying to push an empty key.");
    int ret;
    ret = dict_query(D, key);
    if (ret != -1) {
// D->count[ret]++;
        return ret;
    }
    if (D->n == D->m) {
        int new_m = D->m == 0 ? 1024 : D->m<<1;
        uint32_t *tmp_count = realloc(D->count, sizeof(uint32_t)*new_m);
        if (tmp_count == NULL) errabort("Failed to allocate memory.");
        D->count = tmp_count;
        int i;
        for (i = D->n; i < new_m; ++i) D->count[i] = 0;
        char **tmp_name = realloc(D->name, sizeof(char*)*new_m);
        if (tmp_name == NULL) errabort("Failed to allocate memory.");
        D->name = tmp_name;
        if (D->assign_value_flag) {
            void **tmp_value = realloc(D->value, sizeof(void *)*new_m);
            if (tmp_value == NULL) errabort("Failed to allocate memory.");
            D->value = tmp_value;
            int j;
            for (j = D->n; j < new_m; ++j) D->value[j] = NULL; // reset
        }
        D->m = new_m;
    }

    ret = idx >= 0 ? idx : D->n;

    D->name[D->n] = strdup(key);
    khint_t k;
    int ret0;
    k = kh_put(name, D->dict, D->name[D->n], &ret0);

    kh_val(D->dict, k).idx = D->n;
    kh_val(D->dict, k).ref = ret;

    D->n++;
    return ret;
}

int dict_push(struct dict *D, char const *key)
{
    if (key == NULL) {
        warnings("Try to push empty key! Skip ..");
        return -1;
    }
    int ret;
    ret = dict_push0(D, key, -1);
    D->count[ret]++;
    return ret;
}
// push new key without increase count
int dict_push1(struct dict *D, char const *key)
{
    int ret;
    ret = dict_push0(D, key, -1);
    return ret;
}

// point error tolerated key to original key
int dict_push2(struct dict *D, char const *key, int idx)
{
    if (key == NULL) {
        warnings("Try to push empty key! Skip ..");
        return -1;
    }
    int ret;
    ret = dict_push0(D, key, idx);
    if (ret != -1) return ret; // already present
    D->count[ret]++;
    return ret;
}
// allow space in the key??
int dict_read(struct dict *D, const char *fname, int allow_space)
{
    gzFile fp;
    fp = gzopen(fname, "r");
    if (fp == NULL) errabort("%s : %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    while (ks_getuntil(ks, 2, &str, &ret)>=0){
        if (str.l == 0) continue;
        if (str.s[0] == '#') continue;
        if (strcmp(str.s, "Barcode") == 0) {
            warnings("\"Barcode\" in %s looks like a title, skip it. ", fname);
            continue; // emit header
        }
        if (allow_space == 0) {
            char *p = str.s;
            char *e = str.s + str.l;
            while (p < e && !isspace(*p)) p++;
            *p = '\0';
        }
        dict_push1(D, str.s); // v0.4, init whitelist but not increase count
    }
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);

    if (dict_size(D) == 0) return 1;
    return 0;
}
// read tab with value, col 1 is key, col 2 is val
int dict_read2(struct dict *D, const char *fname, int *val)
{
    *val = 0;
    gzFile fp;
    fp = gzopen(fname, "r");
    if (fp == NULL) errabort("%s : %s.", fname, strerror(errno));
    kstream_t *ks = ks_init(fp);
    kstring_t str = {0,0,0};
    int ret;
    while (ks_getuntil(ks, 2, &str, &ret)>=0){
        if (str.l == 0) continue;
        if (str.s[0] == '#') continue;
        if (strcmp(str.s, "Barcode") == 0) {
            warnings("\"Barcode\" in %s looks like a title, skip it. ", fname);
            continue; // emit header
        }

        int n;
        int *s = ksplit(&str, '\t', &n);
        int idx = dict_push1(D, str.s + s[0]); // init whitelist but not increase count

        if (n > 1) {
            if (*val == 0) {
                dict_set_value(D);
                *val = 1;
            }

            char *v = strdup(str.s + s[1]);
            dict_assign_value(D, idx, v);
        }
        free(s);
    }
    if (str.m) free(str.s);
    ks_destroy(ks);
    gzclose(fp);

    if (dict_size(D) == 0) return 1;
    return 0;
}

char **dict_names(struct dict *D)
{
    return D->name;
}
// hamming distance
static int check_similar(char *a, char *b, int mis)
{
    int l0 = strlen(a);
    assert(strlen(b) == l0);
    int i;
    int m = 0;
    for (i = 0; i < l0; ++i) {
        if (a[i] != b[i]) m++;
        if (m >mis) return 1;
    }
    return 0;
}
// Consider errors during PCR and sequencing, sometime barcodes may contain one or more mismatches.
// this function try to compare each value, allow 1 mismatch, and return the most likely key.
char *dict_most_likely_key(struct dict *D)
{
    uint32_t count = 0;
    char *key = NULL;
    int i;
    for (i = 0; i < D->n; ++i) {
        if (key == NULL) {
            key = D->name[i];
            count = D->count[i];
        }
        else {
            if (check_similar(key, D->name[i], 1) == 0) {
                if (count < D->count[i]) {
                    count = D->count[i];
                    key = D->name[i];
                }
            }
            else {
                return NULL;
            }
        }
    }

    return key;
}
int dict_del(struct dict *D, const char *key)
{
    struct dict_val *v = dict_query0(D, key);
    if (v == NULL) return -1;
    v->ref = -2;
    return 0;
}

int dict_exist(struct dict *D, const char *key)
{
    struct dict_val *v = dict_query0(D, key);
    if (v == NULL) return -1; // not exist
    return v->ref; // deleted or ref to
}

struct dict *dict_dup(struct dict *D)
{
    struct dict *new = dict_init();
    int i;
    for (i = 0; i < D->n; ++i) {
        dict_push(new, D->name[i]);
    }

    return new;
}
