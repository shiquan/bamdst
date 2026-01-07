CC=		gcc
CFLAGS=		-g -Wall -O2 -std=gnu99 
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

# htslib support - use pkg-config if available, otherwise use default paths
HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null || echo "")
HTSLIB_LIBS := $(shell pkg-config --libs htslib 2>/dev/null || echo "-lhts")

# Source files
SOURCES=	xamdst.c kstring.c bedutil.c commons.c
PROG=		xamdst
INCLUDES=	-I. $(HTSLIB_CFLAGS)

.SUFFIXES:.c .o
.PHONY: all

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:clean $(PROG)

.PHONY:all clean

xamdst: $(SOURCES:.c=.o)
		$(CC) $(CFLAGS) -o $@ $(SOURCES:.c=.o) $(INCLUDES) -lm $(HTSLIB_LIBS) -lz -lpthread

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a target.dep *.plot *.report *.tsv.gz uncover.bed

