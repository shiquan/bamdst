/* retrieveNregion.c -- retrieve the N's regions in the fasta file
 * Author: shiquan@genomics.cn
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdint.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static char * out = NULL;


int main(int argc, char **argv)
{
    int n, i;
    BGZF *bgzf;
    FILE *fp = stdout;
    n = argc - optind;
    if ( n > 0 ) {
	bgzf = bgzf_open(argv[1], "r");
	if (!bgzf) {
	    fprintf(stderr, "load fasta file failed %s\n", argv[0]);
	    return -1;
	}
	if ( n == 2 ) {
	    fp = fopen(argv[2], "w");
	    if (!fp) {
		fprintf(stderr, "write bed file failed %s\n", argv[1]);
		return -1;
	    }
	}
    } else {
	fprintf(stderr, "Usage: retrieveNreg <in.fasta> [out.bed]\n");
	return -1;
    }
    int c;
    uint64_t l1, l2, l3;
    int m_name, l_name;
    char *name;
    m_name = l_name = 0;
    name = 0;
    l1 = l2 = l3 = 0;
    while ( (c = bgzf_getc(bgzf)) >= 0)
    {
	if (c == '\n') continue;
	if (c == '>') { // fasta header
	    if (l2 != l3) fprintf(fp,"%s\t%llu\t%llu\n", name, l1, l3);
	    l1 = l2 = l3 = 0;
	    l_name = 0;
	    while ( (c = bgzf_getc(bgzf)) >= 0 && c!= '\n')
	    {
		if ( m_name < l_name +2) {
		    m_name = l_name+2;
		    kroundup32(m_name);
		    name = (char*)realloc(name, m_name);
		}
		name[l_name++] = c;
	    }
	    name[l_name] = '\0';
	    if (c != '\n') while ( (c=bgzf_getc(bgzf))>0 && c != '\n');
	} else {
	 
	    if ( c == 'N' ) {
		if (l1 != l2) l1 = l2;
		l3++;
	    } else {
		if (l2 != l3) {
		    fprintf(fp, "%s\t%llu\t%llu\n", name, l1, l3);
		    l2 = ++l3;
		} else {
		    l3++; l2++;
		}
	    }
	}
    }
    free(name);
    bgzf_close(bgzf);
    fclose(fp);
    return 0;
}
