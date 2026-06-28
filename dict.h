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

#ifndef DICT_H
#define DICT_H

#include "commons.h"

struct dict;

struct dict *dict_init();
struct dict *dict_dup(struct dict *D);

void dict_destroy(struct dict *D);

int dict_query(const struct dict *D, char const *key);
int dict_query2(const struct dict *D, char const *key);

int dict_push(struct dict *D, char const *key);
int dict_push1(struct dict *D, char const *key);
int dict_push2(struct dict *D, char const *key, int idx);
int dict_read(struct dict *D, const char *fname, int allow_space);
int dict_read2(struct dict *D, const char *fname, int *val);

char *dict_name(const struct dict *D, int idx);

int dict_size(const struct dict *D);

uint32_t dict_count_sum(const struct dict *D);

uint32_t dict_count(const struct dict *D, int idx);

char **dict_names(struct dict *D);

char *dict_most_likely_key(struct dict *D);

void dict_set_value(struct dict *D);
void *dict_query_value(struct dict *D, int idx);
void *dict_query_value2(struct dict *D, const char *key);
int dict_assign_value(struct dict *D, int idx, void *val);
int dict_delete_value(struct dict *D, int idx);

int dict_del(struct dict *D, const char *key);
int dict_exist(struct dict *D, const char *key);

#endif
