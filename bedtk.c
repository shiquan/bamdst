#include "bedutil.h"
#include "commons.h"

extern bedHandle_t *bedHand;


int mergeHelp(void)
  {
  fprintf(stderr, "Usage : program   merge   <in1.bed> [in2.bed] ...\n"
	  "==================================================\n"
	  "  -r     add n bp before each regions\n"
	  "  -l     add n bp after each regions\n"
	  "  -o     output file\n"
	  "  -h     show this message\n"
	  );
  return 1;
  }

int mergeBed(int argc, char * argv[]) 
  {
  int help = 0;
  int n, i;
  char *out = 0;
  int add1 = 0, add2 = 0;
  while ((n = getopt(argc, argv, "o:hr:l:")) >= 0)
    {
    switch(n)
      {
      case 'o': out = optarg; break;
      case 'h': help = 1; break;
      case 'r': add1 = atoi(optarg); break;
      case 'l': add2 = atoi(optarg); break;
      }
    }
  if(help) return mergeHelp();
  assert(add1 >= 0 && add2 >= 0);
  n = argc - optind;
  regHash_t * reghash;
  reghash = kh_init(reg);
  for (i = 0; i < n; ++i)
    {
    bedHand->read(argv[optind+i], reghash, add1, add2);
    }
  inf_t *itmp = bedHand->stat(reghash);
  bedHand->merge(reghash);
  inf_t *inf = bedHand->stat(reghash);
  inf->region = itmp->total - inf->total;
  if (out) bedHand->save(out, reghash);
  writeout("Merged %u regions.\n"
	   "Total regions is %u.\n"
	   "The length of the regions is %u bp.\n",
	   inf->region, inf->total, inf->length);
  bedHand->destroy(reghash, destroy_void);
  freemem(itmp);
  freemem(inf);
  return 1;
  }

void uniqHelp()
  {
  fprintf(stderr,"Usage: uniq <in1.bed> <in2.bed>\n"
	  "   -o  Set output.\n"
	  "   -h  See this information."
	  );
  }

int uniqBed(int argc, char * argv[])
  {
  int help = 0;
  int n, i;
  char * out = 0;
  while ((n = getopt(argc, argv, "o:h")) >= 0)
    {
    switch(n)
      {
      case 'o': out = optarg; break;
      case 'h': help = 1; break;
      default: uniqHelp();
      }
    if (help) uniqHelp();
    }
  n = argc - optind;
  if (n < 2) errabort("At least set 2 files!");
  regHash_t * reghash = kh_init(reg);
  bedHand->read(argv[optind], reghash, 0, 0);
  bedHand->merge((void*)reghash);
  for (i = 1; i < n; ++i)
    {
    regHash_t * reghash1 = kh_init(reg);
    bedHand->read(argv[optind+i],reghash1, 0, 0);
    bedHand->merge(reghash1);	
    bedHand->uniq(reghash, reghash1); 
    bedHand->destroy(reghash1, destroy_void);
    }
  inf_t *inf = bedHand->stat(reghash);
  if (out) bedHand->save(out, reghash);
  writeout("The length of the uniq regions is %d bp.\n", inf->length);
  freemem(inf);
  bedHand->destroy(reghash, destroy_void);
  return 1;
  }

void diffHelp()
  {
  fprintf(stderr,"Usage: diff <in1.bed> <in2.bed>\n"
	  "   -o  Set output file.\n"
	  "   -h  See this message.\n"
	  );
  }

int diffBed(int argc, char * argv[]) 
  {
  int help = 0;
  char *out = 0;
  int n, i;
  while ((n = getopt(argc, argv, "o:h")) >= 0)
    {
    switch(n)
      {
      case 'o': out = optarg; break;
      case 'h': help = 1; break;
      default: diffHelp();
      }
    if (help) diffHelp();
    }
  n = argc - optind;
  if(n < 2) errabort("At least 2 files!");
  regHash_t * reghash, * reghash1;
  reghash = kh_init(reg);
  reghash1 = kh_init(reg);
  bedHand->read(argv[optind], (void *)reghash, 0, 0);
  bedHand->merge(reghash);
  for (i = 1; i < n; ++i)
    {
    bedHand->read(argv[optind+i], (void *)reghash1, 0, 0);
    }
  bedHand->merge(reghash1);
  bedHand->diff(reghash, reghash1);
  inf_t * inf = bedHand->stat(reghash);
  writeout("There are %d bp in the first file but not in others.", inf->length);
  freemem(inf);
  bedHand->destroy(reghash1, destroy_void);
  if (out) bedHand->save(out, reghash);
  bedHand->destroy(reghash, destroy_void);
  return 1;
  }

void trimHelp()
  {
  fprintf(stderr,"Usage: trim -r <trim1> -l <trim2> <in1.bed> [in2.bed]\n"
	  "    -o <FILE>      Set output file.\n"
	  "    -r <n>         Trim n bp in the front of each region.\n"
	  "    -l <n>         Trim n bp in the end of each region.\n"
	  "    -h             See this message.\n"
	  " Notices: the regions will be merged firstly."
	  );
  }

int trimBed(int argc, char *argv[])
  {
  int help = 0;
  int n;
  int i;
  char *out = 0;
  int trim1 = 0;
  int trim2 = 0;
  while ((n = getopt(argc, argv, "o:r:l:h")) >= 0)
    {
    switch(n)
      {
      case 'o' : out = optarg; break;
      case 'l' : trim1 = atoi(optarg); break;
      case 'r' : trim2 = atoi(optarg); break;
      case 'h' : help = 1; break;
      default :
	errabort("%c : unknown option!", (char)n);
	break;
      }
    if (help) trimHelp();
    }
  if (trim1 == 0 && trim2 == 0) errabort("You must set a trim size!");
  if (trim1 < 0 || trim2 < 0) errabort("Trim size must be a postive int!");
  n = argc - optind;
  if (n < 1) errabort("Sucker! Set the bed file(s)!");
  regHash_t *rghsh = kh_init(reg);
  for (i = 0; i < n; ++i)
    {
    bedHand->read(argv[optind+i], rghsh, 0, 0);	 
    }
  bedHand->merge(rghsh);
  inf_t *itmp = bedHand->stat(rghsh);
  bedHand->trim(rghsh, trim1, trim2);
  inf_t *inf = bedHand->stat(rghsh);
  if (out) bedHand->save(out, rghsh);
  unsigned trim_length = itmp->length - inf->length;
  writeout("Trimmed %d bp\n"
	   "The length of whole regions is %d bp\n"
	   ,trim_length, inf->length);
  freemem(itmp);
  freemem(inf);
  bedHand->destroy(rghsh, destroy_void);
  return 1;
  }

static int usage(void)
  {
  fprintf(stderr, "bedutils   merge -o <FILE> <in1.bed [in2.bed ...]>\n"
	  "           diff  -o <FILE> <in1.bed> <in2.bed> [in3.bed ...]\n"
	  "           uniq  -o <FILE> <in1.bed> <in2.bed> [in3.bed ...]\n"
	  "           trim  -l xx -n xx <in1.bed> [in3.bed ...]\n"
	  );
  return 1;
  }

int main(int argc, char *argv[]) 
  {
  if (argc < 2) return usage();
  else if (STREQ(argv[1], "merge")) return mergeBed(argc-1, argv+1);
  else if (STREQ(argv[1], "diff")) return diffBed(argc-1, argv+1);
  else if (STREQ(argv[1], "uniq")) return uniqBed(argc-1, argv+1);
  else if (STREQ(argv[1], "trim")) return trimBed(argc-1, argv+1);
  else return usage();
  return 1;
  }

