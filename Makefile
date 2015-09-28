CC=		gcc
CFLAGS=		-g -Wall -O2 
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DBGZF_CACHE
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o bam_index.o sam_header.o bedutil.o commons.o
PROG=		bamdst
INCLUDES=	-Isamlib/ -I.
SUBDIRS=	. samlib
LIBPATH=        -L. 

.SUFFIXES:.c .o
.PHONY: all

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:clean $(PROG) 

.PHONY:all  clean

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bamdst:lib $(AOBJS) samlib/bam.h
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LDFLAGS) bamdst.c $(LIBPATH) $(INCLUDES) -lm -lbam -lz

bgzf.o:bgzf.c bgzf.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) bgzf.c -o $@

kstring.o:kstring.c kstring.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) kstring.c -o $@

bam.o:samlib/bam.h samlib/bam_endian.h kstring.h samlib/sam_header.h samlib/bam.c
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) samlib/bam.c -o $@

bam_import.o:samlib/bam.h kseq.h khash.h samlib/bam_import.c
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) samlib/bam_import.c -o $@

bam_index.o:samlib/bam.h khash.h ksort.h  samlib/bam_endian.h samlib/bam_index.c  
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) samlib/bam_index.c -o $@

sam_header.o:samlib/sam_header.h khash.h samlib/sam_header.c
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) samlib/sam_header.c -o $@

bam_aux.o:samlib/bam.h samlib/bam_aux.c khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) samlib/bam_aux.c -o $@

commons.o:commons.c commons.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) commons.c -o $@	

bedutil.o:bedutil.c bedutil.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) bedutil.c -o $@	

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a target.dep *.plot *.report *.tsv.gz uncover.bed
