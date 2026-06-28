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

#ifndef NUMBER_H
#define NUMBER_H

int get_numbase(const char *s);
int get_numbase_l(const char *s, int l);

int is_ieee_magic_val(const char *val);

double nondec2num(char *str, int length);

int check_num_likely(const char *str);
int check_num_likely_l(const char *str, int l);
int check_char_num(const char x);

double force2num(char *str);
double force2num_l(char *str, int l);

int str2int(const char *str);
int str2int_l(const char *str, int l);

int human2int(const char *str);
int genome2int(const char *str);

#endif
