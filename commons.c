#include "commons.h"

static void stderr_print(char *format, ...)
  {
  va_list args;
  va_start(args, format);
  if (format != NULL)
    {
    fflush(stdout);
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    }
  va_end(args);
  }

void writeout(char *format, ...)
  {
  va_list args;
  va_start(args, format);
  fflush(stdout);
  if (format != NULL)
    {
    vfprintf(stdout, format, args);
    }
  va_end(args);
  }

void warnings(char *format, ...)
  {
  va_list args;
  va_start(args, format);
  if (format != NULL)
    {
    fflush(stdout);
    fprintf(stderr, "[Warnings] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    }
  va_end(args);
  }

void errabort(char *format, ...)
  {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "@@ ERROR ABORT @@ \n");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  va_end(args);
  exit(-1);
  }

void debug(char const *format, ...)
  {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "[DEBUG] ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  va_end(args);
  }

void LOG(char * format, ...)
  {
  va_list args;
  va_start(args, format);
  if (format != NULL)
    {
    fflush(stdout);
    fprintf(stderr, "[LOG] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    fflush(stderr);
    }
  va_end(args);
  }

FILE *
stream_open(char const *file, char const *mode)
  {
  FILE *fp = NULL;
  if (*mode == 'r')
    {
    if (STREQ(file, "-"))
      fp = stdin;
    else
      fp = fopen(file, mode);
    //fadvise(fp, FADVISE_SEQUENTIAL);
    }
  else if (*mode == 'w')
    {
    if (file && ftruncate(STDOUT_FILENO, 0) != 0)
      errabort("%s : error truncating", file);
    fp = stdout;
    }
  else
    errabort("Unexpected mode passed to stream_open: %s", mode);
  return fp;
  }

FILE *
open_wfile(char const *file)
  {
  FILE *fp;
  if (STREQ(file, "-"))
    fp = stderr;
  else
    {
    fp = fopen(file, "w");
    if ( isNull(fp) )
      errabort("%s : %s", file, strerror(errno));
    }
  return fp;
  }
  
FILE *
open_rfile(char const *file)
  {
  return stream_open(file, "r");
  }

gzFile 
safe_gzopen(char const *file)
  {
  gzFile fp;
  if (STREQ(file, "-"))
    fp = gzdopen(fileno(stdin), "r");
  else
    fp = gzopen(file, "r");
  if ( isNull(fp) )
    errabort("%s : %s", file, strerror(errno));
  return fp;
  }

void write_file(FILE *fp, char *format, ...)
  {
  va_list args;
  va_start(args, format);
  vfprintf(fp, format, args);
  va_end(args);
  }

#define NEEDMEM_LIMIT 500000000
#define CHECK_BIT (sizeof(size_t) / 4)
/* 128*8*1024*1024 == 1073741824 == 2^30 on 32 bit machines,size_t == 4 bytes*/
/* on 64 bit machines, size_t = 8 bytes, 2^30 * 2 * 2 * 2 * 2 = 2^34 == 16 Gb */
static size_t maxalloc =
  (size_t)1073741824 * CHECK_BIT *CHECK_BIT *CHECK_BIT *CHECK_BIT;
			
void *needmem(size_t size)
  {
  void *pt;
  if (isZero(size) || size > NEEDMEM_LIMIT)
    errabort("[NeedMEM] trying to allocate %llu bytes (limits: %llu)",
	     (unsigned long long)size, NEEDMEM_LIMIT);
  if (isNull(pt = malloc(size)))
    errabort("[needMem] Out of memory - request size %llu bytes : %s",
	     (unsigned long long)size, strerror(errno));
  memset(pt, 0, size);
  return pt;
  }

void *neednmem(size_t n, size_t m)
  {
  return needmem(n * m);
  }

void freemem(void *pt)
  {
  if (isNull(pt)) return;
  free(pt);
  }

void *needlargemem(size_t size)
  {
  void *pt;
  if (isZero(size) || size >= maxalloc)
    errabort("[needLargeMem] trying to allocate %llu bytes (limit: %llu)",
	     (unsigned long long)size, (unsigned long long)maxalloc);
  if (isNull(pt = malloc(size)))
    errabort("[needLargeMem] Out of memory - request size %llu bytes: %s",
	     (unsigned long long)size, strerror(errno));
  return pt;
  }

void *resizemem(void *vp, size_t size)
  {
  void *pt;
  if ( isNull(vp))
    errabort("%s : %d [resizeMEM] Trying to allocate an empty point."
	     , __FILE__, __LINE__);
  if (isZero(size) || size > NEEDMEM_LIMIT)
    errabort("[resizeMEM] trying to allocate %llu bytes (limits: %llu)",
	     (unsigned long long)size, NEEDMEM_LIMIT);
  if (isNull(pt = realloc(vp, size)))
    errabort("[resizeMEM] Out of memory - request size %llu bytes : %s",
	     (unsigned long long)size, strerror(errno));
  return pt;
  }

void *resizelargemem(void *vp, size_t size)
  {
  void *pt;
  if ( isNull(vp))
    errabort("%s : %d [resizeMEM] Trying to allocate an empty point."
	     , __FILE__, __LINE__); // need use macro improve this
	
  if (isZero(size) || size > maxalloc)
    errabort("[resizeLargeMem] trying to allocate %llu bytes (limit: %llu)",
	     (unsigned long long)size, (unsigned long long)maxalloc);
  if (isNull(pt = realloc(vp, size)))
    errabort("[resizeLargeMem] Out of memory - request size %llu bytes: %s",
	     (unsigned long long)size, strerror(errno));
  return pt;
  }

void *enlarge_empty_mem(void *vp, size_t old_size, size_t new_size)
  {
  void *pt;
  pt = resizemem(vp, new_size);
  memset(((char*)pt) + old_size, 0, new_size - old_size);
  return pt;
  }

void *enlarge_empty_largemem(void *vp, size_t old_size, size_t new_size)
  {
  void *pt;
  pt = resizelargemem(vp, new_size);
  memset(((char*)pt)+old_size, 0, new_size-old_size);
  return pt;
  }
