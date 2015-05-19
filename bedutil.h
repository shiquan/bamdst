//
//  bedutil.h
//
// Created by shiquan on 5/18/14.
// Copyright (c) 2014 ___BGIRESEARCH___. 
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// 
//
// This program is designed to handle bed file, inspired by liheng's bedidx.c
//
#ifndef BEDUTILS_H
#define BEDUTILS_H
#include <stdint.h>
#include "khash.h"

typedef struct
  {
  int m; // used-length
  int n; // alloced-length
  int flag;
  uint64_t *a;
  int *idx;
  void *data; // self defined data
  }
bedreglist_t;

typedef struct
  {
  uint32_t length;
  void *data;
  }
chromosome_t;

KHASH_MAP_INIT_STR(reg, bedreglist_t)
KHASH_MAP_INIT_STR(chr, chromosome_t)

typedef kh_reg_t regHash_t;
typedef kh_chr_t chrHash_t;

typedef struct
  {
  int is_empty;
  regHash_t *reg;
  chrHash_t *chr;
  char **seq_names; // chromosome names
  int n_seq; // number of names
  uint32_t region; // region count
  uint32_t length; // region length
  uint32_t total;
  }
bedaux_t;

typedef struct
  {
  uint32_t region;
  uint32_t length;
  uint32_t total;
  }
inf_t;

typedef void (*bedvoid_destroy)(void *data);
void destroy_void(void *data);

typedef void (*handle_func)(regHash_t *bed);

struct _bedHandle
  {
  /* handle inf_t */
  inf_t * (*init)();
  /* handle whole struct */ 
  void (*read)(const char *fn, regHash_t *reg, int a, int b, int *ret);
  void (*merge)(regHash_t *bed);
  void (*uniq)(regHash_t *bed1, regHash_t *bed2);
  void (*diff)(regHash_t *bed1, regHash_t *bed2);
  void (*trim)(regHash_t *bed, int trim1, int trim2);
  /* handle each chromosome */
  bedreglist_t * (*reg_add)(bedreglist_t *reg1, bedreglist_t *reg2);
  void (*reg_merge)(bedreglist_t *reg1);
  bedreglist_t * (*reg_uniq)(bedreglist_t *reg1, bedreglist_t *reg2);
  bedreglist_t * (*reg_diff)(bedreglist_t *reg1, bedreglist_t *reg2);
  bedreglist_t * (*reg_trim)(bedreglist_t *reg, int trim1, int trim2);
  void (*sort_reg) (bedreglist_t *reg);
  /* save file */
  void (*save)(const char *fn, regHash_t *bed);
  /* destroy */
  void (*destroy)(regHash_t *rghsh, bedvoid_destroy func);
  /* stat */
  inf_t * (*stat)(regHash_t *bed);
  void (*base1to0) (regHash_t *bed);
  void (*base0to1) (regHash_t *bed);
  void (*pipeout)(regHash_t *bed);

  /* if the region is out of chromosome return khiter, else 
     return -1 : normal 
     return -2 : the chomosome is not contaioned in chrHash
  */
  int (*check_length)(regHash_t * bed, chrHash_t *chr);
  };

typedef struct _bedHandle bedHandle_t;

#endif
