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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "number.h"

// Return the base to use for the number in 's', this function adapt from gawk/node.c
int get_numbase(const char *s)
{
    assert(s);
    if (strlen(s) < 2)
        return 10;

    const char *str = s;
    int dec_point = '.';
    if ( str[0] != '0' )
        return 10;

    // Leading with 0x or 0X
    if ( str[1] == 'x' || str[1] == 'X')
        return 16;

    // Numbers with dec_number, 'e' or 'E' are decimal.
    for (; *str != '\0'; str++) {
        if ( *str == 'e' || *str == 'E' || *str == dec_point)
            return 10;
        else if ( !isdigit((unsigned char)*str))
            break;
    }

    if ( !isdigit((unsigned char)s[1]) || s[1] == '8' || s[1] == '9')
        return 10;
    return 8;
}

int get_numbase_l(const char *s, int l)
{
    assert(s);
    if (strlen(s) < 2)
        return 10;
    const char *str = s;
    int i = 0;
    int dec_point = '.';
    if ( str[0] != '0' )
        return 10;

    // Leading with 0x or 0X
    if ( str[1] == 'x' || str[1] == 'X')
        return 16;

    // Numbers with dec_number, 'e' or 'E' are decimal.
    for (; *str != '\0' && i < l; str++, i++) {
        if ( *str == 'e' || *str == 'E' || *str == dec_point)
            return 10;
        else if ( !isdigit((unsigned char)*str))
            break;
    }

    if ( !isdigit((unsigned char)s[1]) || s[1] == '8' || s[1] == '9')
        return 10;
    return 8;
}
// is_ieee_magic_val - return 1 for +inf, -inf, +nan, -nan
int is_ieee_magic_val(const char *val)
{
    if ( strlen(val) != 4)
        return 0;
    return ( (val[0] == '+' || val[0] == '-')  &&
             ( ( (val[1] == 'i' || val[1] == 'I')
                 && (val[2] == 'n' || val[2] == 'N')
                 && (val[3] == 'f' || val[3] == 'F') ) ||
               ( (val[1] == 'n' || val[1] == 'N')
                 && (val[2] == 'a' || val[2] == 'A')
                 && (val[3] == 'n' || val[3] == 'N') ) ) );
}

double nondec2num(char *str, int length)
{
    double retval = 0.0;
    char save;
    short val;
    char *start = str;
    if ( *str == '0' && (str[1] == 'x' || str[1] == 'X')) {
        if (length <= 2)
            return 0.0;
        for ( str+=2, length -=2; length > 0; length--, str++) {
            switch(*str) {
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                    val = *str - '0';
                    break;

                case 'a':
                case 'b':
                case 'c':
                case 'd':
                case 'e':
                case 'f':
                    val = *str - 'a' + 10;
                    break;

                case 'A':
                case 'B':
                case 'C':
                case 'D':
                case 'E':
                case 'F':
                    val = *str - 'A' + 10;
                    break;
                default:
                    goto nondec2num_done;
            }
            retval = ( retval *16) + val;
        }
    } else if (*str == '0') {
        for (; length > 0; length -- ) {
            if ( !isdigit((unsigned char)*str))
                goto nondec2num_done;
            else if ( *str == '8' || *str == '9') {
                str = start;
                goto nondec2num_decimal;
            }
            retval = (retval*8) + (*str-'0');
            str++;
        }

    } else {
      nondec2num_decimal:
        save = str[length];
        str[length] = '\0';
        retval = strtod(str, NULL);
        str[length] = save;
    }
  nondec2num_done:
    return retval;
}

int check_num_likely(const char *str)
{
    if ( str == NULL || strlen(str) == 0)
        return 0;
    for (; str != NULL &&  *str != '\0'; ++str) {
        if ( isspace(*str) )
            continue;
        if ( strlen(str) > 2) {
            if ( str[0] == '0' && (str[1] == 'X' || str[1] == 'x') )
                str += 2;
        }

        switch(*str) {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            case 'a':
            case 'b':
            case 'c':
            case 'd':
            case 'e':
            case 'f':
            case 'A':
            case 'B':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'x':
            case 'X':
            case '.':
            case '+':
            case '-':
            case ' ':
                continue;
                break;

            default:
                return 0;
        }
    }
    return 1;
}
int check_num_likely_l(const char *str, int l)
{
    if ( str == NULL || strlen(str) == 0)
        return 0;
    int i;
    for ( i = 0; str != NULL &&  *str != '\0' && i < l; ++str, ++i) {
        if ( isspace(*str) )
            continue;
        if ( strlen(str) > 2) {
            if ( str[0] == '0' && (str[1] == 'X' || str[1] == 'x') )
                str += 2;
        }

        switch(*str) {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            case 'a':
            case 'b':
            case 'c':
            case 'd':
            case 'e':
            case 'f':
            case 'A':
            case 'B':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'x':
            case 'X':
            case '.':
            case '+':
            case '-':
            case ' ':
                continue;
                break;

            default:
                return 0;
        }
    }
    return 1;
}
int check_char_num(const char x)
{
    switch(x) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
            return 1;
        default:
            return 0;
    }
    return 0;
}

double force2num(char *str)
{
    int base = 0;
    base = get_numbase(str);
    if ( base == 10 )
        return strtod(str, NULL);
    else
        return nondec2num(str, strlen(str));
}

double force2num_l(char *str, int l)
{
    int base = 0;
    base = get_numbase_l(str, l);
    char *temp = strndup(str, l);
    // char *pend;
    if ( base == 10 )  {
        double b = strtod(temp, NULL );
        free(temp);
        return b;
    } else {
        free(temp);
        return nondec2num(str, l);
    }
}
int str2int(const char *str)
{
    return (int)force2num((char*)str);
}
int str2int_l(const char *str, int l)
{
    return (int)force2num_l((char*)str, l);
}

int human2int(const char *str)
{
    char *q;
    int m = strtol(str, &q, 0);
    if (*q == 'k'||*q=='K') m<<=10;
    else if (*q == 'm'||*q=='M') m<<=20;
    else if (*q == 'g'||*q=='G') m<<=30;
    return m;
}

int genome2int(const char *str)
{
    char *q;
    int m = strtol(str, &q, 0);
    if (*q == 'k'||*q=='K') m = m *1000;
    else if (*q == 'm'||*q=='M') m = m*1000000;
    else if (*q == 'g'||*q=='G') m = m*1000000000;
    return m;
}
