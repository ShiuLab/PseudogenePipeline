/* copyright (c) 1996, 1997, 1998, 1999, 2002  William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: comp_lib.c,v 1.100 2007/04/26 18:36:36 wrp Exp $ */

/*
 * Concurrent read version
 *
 *	Feb 20, 1998 modifications for prss3
 *
 *	December, 1998 - DNA searches are now down with forward and reverse
 *			 strands
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include <limits.h>
#include <float.h>
#include <math.h>

#ifdef UNIX
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#endif

#include "defs.h"
#include "mm_file.h"

#include "mw.h"			/* defines beststr */
#include "structs.h"		/* mngmsg, libstruct */
#include "param.h"		/* pstruct, thr_str, buf_head, rstruct */

#define XTERNAL
#include "uascii.h"

char *mp_verstr="34.26";

/********************************/
/* global variable declarations */
/********************************/
char gstring2[MAX_STR];                  /* string for label */
char gstring3[MAX_STR];
char hstring1[MAX_STR];

extern int max_workers;

#ifdef SUPERFAMNUM
int nsfnum;
int sfnum[10];
extern int sfn_cmp(int *q, int *s);
int nsfnum_n;
int sfnum_n[10];
#endif

/********************************/
/* extern variable declarations */
/********************************/
extern char *prog_func;		/* function label */
extern char *verstr, *iprompt0, *iprompt1, *iprompt2, *refstr;

/********************************/
/*extern function declarations  */
/********************************/
/* open sequence file (getseq.c) */
extern int getseq(char *filen, int *sascii,
		  unsigned char *seq, int maxs,
		  char *libstr, int n_libstr,
		  long *sq0ff);

struct lmf_str *openlib(char *, int, int *, int, struct lmf_str *);

void set_shuffle(struct mngmsg m_msg);
void closelib(struct lmf_str *m_fptr);

void irand(int);
int nrand(int);

extern int ann_scan(unsigned char *, int, struct mngmsg *, int );
extern int scanseq(unsigned char *seq, int n, char *str);
extern void re_ascii(int *qascii, int *sascii);
extern int recode(unsigned char *seq, int n, int *qascii, int nsq);
extern void revcomp(unsigned char *seq, int n, int *c_nt);

extern void init_ascii(int is_ext, int *sascii, int is_dna);
extern void qshuffle(unsigned char *aa0, int n0, int nm0);
extern void free_pam2p(int **);

/* initialize environment (doinit.c) */
extern void initenv (int argc, char **argv, struct mngmsg *m_msg,
		     struct pstruct *ppst, unsigned char **aa0);

/* print timing information */
extern void ptime (FILE *, time_t);

#ifdef COMP_MLIB 
#define QGETLIB (q_file_p->getlib)
#endif

#define GETLIB (m_file_p->getlib)

/* calculation functions */
extern void
init_work(unsigned char *aa0, int n0,
	  struct pstruct *ppst, void **f_arg );
#ifndef COMP_THR
extern void
do_work(unsigned char *aa0, int n0, unsigned char *aa1, int n1, int frame, 
	struct pstruct *ppst, void *f_str, int qr_flg, struct rstruct *rst);
#endif

extern void
close_work(unsigned char *aa0, int n0, struct pstruct *ppst, void **f_arg);
extern void
get_param (struct pstruct *pstr, char *pstring1, char *pstring2);

#ifdef COMP_THR
#ifndef PRSS
void
save_best(struct buf_head *cur_buf, struct mngmsg, struct pstruct pst, 
	  FILE *fdata, int *, struct hist_str *, void **);
#else
void
save_best(struct buf_head *cur_buf, struct mngmsg, struct pstruct pst, 
	  FILE *fdata, int *, struct hist_str *, void **, int *, int *);
#endif
#endif

/* statistics functions */
extern int
process_hist(struct stat_str *sptr, int nstat, 
	     struct mngmsg m_msg,
	     struct pstruct pst,
	     struct hist_str *hist, void **, int);
extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );
extern double (*find_zp)(int score, double escore, int length, double comp,void *);

void last_stats(const unsigned char *, int, 
		struct stat_str *sptr, int nstats,
		struct beststr **bestp_arr, int nbest,
		struct mngmsg m_msg, struct pstruct pst, 
		struct hist_str *histp, void *);

int last_calc( unsigned char **a0, unsigned char *a1, int maxn,
	       struct beststr **bestp_arr, int nbest,
	       struct mngmsg m_msg, struct pstruct *ppst, 
	       void **f_str, void *rs_str);

void scale_scores(struct beststr **bestp_arr, int nbest,
		  struct db_str,struct pstruct pst, void *);

#ifndef COMP_THR
extern int shuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int, int *);
#endif

extern void set_db_size(int, struct db_str *, struct hist_str *);

/* display functions */
extern void
showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1,
	  int maxn, struct beststr **bestp_arr, int nbest,
	  int qlib, struct mngmsg *m_msg,struct pstruct pst,
	  struct db_str db, char *gstring2, void **f_str);

extern void
showalign (FILE *fp, unsigned char **aa0, unsigned char *aa1,
	   int maxn, struct beststr **bestp_arr, int nbest,
	   int qlib, struct mngmsg m_msg,struct pstruct pst,
	   char *gstring2, void **f_str);

/* misc functions */
void h_init(struct pstruct *, struct mngmsg *, char *);		/* doinit.c */
void last_init(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void last_params(unsigned char *, int, struct mngmsg *, struct pstruct *);

void s_abort(char *, char *);		/* compacc.c */

/* initfa/sw.c */
void resetp(struct mngmsg *, struct pstruct *); 

void gettitle(char *, char *, int);	/* nxgetaa.c */
void libchoice(char *lname, int, struct mngmsg *); /* lib_sel.c */
void libselect(char *lname, struct mngmsg *);	/* lib_sel.c */
void query_parm(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void selectbestz(struct beststr **, int, int);

/* compacc.c */
void prhist(FILE *, struct mngmsg, struct pstruct, 
	    struct hist_str hist, int nstats, struct db_str, char *);
void printsum(FILE *, struct db_str db);
int reset_maxn(struct mngmsg *, int);	/* set m_msg.maxt, maxn from maxl */

FILE *outfd;			/* Output file */

/* this information is global for fsigint() */
extern time_t s_time();			/* fetches time */
time_t tstart, tscan, tprev, tdone;	/* Timing */
#ifdef COMP_MLIB
time_t ttscan, ttdisp;
#endif
time_t tdstart, tddone;

static struct db_str qtt = {0l, 0l, 0};

#ifdef COMP_THR
/***************************************/
/* thread global variable declarations */
/***************************************/

/* functions for getting/sending buffers to threads (thr_sub.c) */
extern void init_thr(int , struct thr_str *, struct mngmsg, struct pstruct *,
		     unsigned char *, int);
extern void start_thr(void);
extern void get_rbuf(struct buf_head **cur_buf, int max_wor_buf);
extern void put_rbuf(struct buf_head *cur_buf, int max_work_buf);
extern void put_rbuf_done(int nthreads, struct buf_head *cur_buf, 
			  int max_work_buf);
#undef XTERNAL
#include "thr.h"
struct buf_head buf_list[NUM_WORK_BUF];
#endif

/* these variables must be global for comp_thr.c so that savebest()
   can use them */

static struct beststr 
    *best,		/* array of best scores */
    *bestp,
    **bestp_arr;	/* array of pointers */
static int nbest;	/* number of best scores */

static struct stat_str *stats, *qstats;	/* array of scores for statistics */

/* these variables are global so they can be set both by the main()
   program and savebest() in threaded mode.
*/
static int nstats, nqstats, kstats;
static double zbestcut;		/* cut off for best z-score */
static int bestfull;		/* index for selectbest() */
static int stats_done=0;	/* flag for z-value processing */
void fsigint();

int
main (int argc, char *argv[]) 
{
  unsigned char *aa0[6], *aa0s, *aa1, *aa1ptr, *aa1s;
  int n1, n1s;  /* n1s needed for PRSS so that when getlib() returns -1 (because no more
		   library sequences, we have a valid n1 for shuffling */

  int *n1tot_ptr=NULL, *n1tot_cur;
  int n1tot_cnt=0;
  int n1tot_v, aa1_loff;

  long qoffset;		/* qoffset is the equivalent of loffset */
  			/* m_msg.sq0off is the l_off equivalent */

  long loffset, l_off;	/* loffset is the coordinate of first residue
			   when lcont > 0; l_off is not used in the
			   main loop, only in showbest and showalign */
  char lib_label[MAX_FN];
  char pgm_abbr[MAX_SSTR];
  char qlabel[MAX_FN];
#ifdef COMP_MLIB
  char q_bline[MAX_STR];
  fseek_t qseek;
  int qlib;
  struct lmf_str *q_file_p;
  int sstart, sstop, is;
#endif
  int id;
  struct lmf_str *m_file_p;

  int t_best, t_rbest, t_qrbest;	/* best score of two/six frames */
  double t_escore, t_rescore, t_qrescore; /* best evalues of two/six frames */
  int i_score;
#ifdef PRSS
  int s_score[3];
  int s_n1;
#endif

  struct pstruct pst;
  void *f_str[6], *qf_str;	/* different f_str[]'s for different
  				   translation frames, or forward,reverse */
  int have_f_str=0;

#ifdef COMP_THR
  long ntbuff;
  int max_buf_cnt, ave_seq_len, buf_siz;
  int max_work_buf;
  struct buf_head *cur_buf;
  struct buf_str *cur_buf_p;
  int nseq;
  struct thr_str *work_info;
#endif

  struct mngmsg m_msg;		/* Message from host to manager */
  int iln, itt;			/* index into library names */
  char rline[MAX_FN];
  char argv_line[MAX_STR];
  int t_quiet;

  struct rstruct  rst;		/* results structure */
  struct rstruct  rrst;		/* results structure for shuffle*/
  int i;

  FILE *fdata=NULL;		/* file for full results */
  char libstr[MAX_UID];		/* string for labeling full results */
  char *libstr_p;		/* choose between libstr and ltitle */
  int n_libstr;			/* length of libstr */
  int jstats;
  int leng;			/* leng is length of the descriptive line */
  int maxn;			/* size of the library sequence examined */
  int maxl;			/* size of library buffer */
  fseek_t lmark;		/* seek into library of current sequence */
  int qlcont;			/* continued query sequence */
  int lcont, ocont, maxt;	/* continued sequence */
  int igncnt=0;			/* count for ignoring sequences warning */
  int ieven=0;			/* tmp for wshuffle */
  double zscore;			/* tmp value */
  char *bp;			/* general purpose string ptr */
  
  /* Initialization */

#if defined(UNIX)
  m_msg.quiet= !isatty(1);
#else
  m_msg.quiet = 0;
#endif

#ifdef PGM_DOC
  argv_line[0]='#'; argv_line[1]='\0';
  for (i=0; i<argc; i++) {
    strncat(argv_line," ",sizeof(argv_line)-strlen(argv_line)-1);
    if (strchr(argv[i],' ')) {
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
    }
    else {
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
    }
  }
  argv_line[sizeof(argv_line)-1]='\0';
#endif

  /* first initialization routine - nothing is known */
  h_init(&pst, &m_msg, pgm_abbr);
  
  m_msg.db.length = qtt.length = 0l;
  m_msg.db.entries = m_msg.db.carry = qtt.entries = qtt.carry = 0;
  m_msg.pstat_void = NULL;
  m_msg.hist.entries = 0;

  for (iln=0; iln<MAX_LF; iln++) m_msg.lb_mfd[iln]=NULL;

  f_str[0] = f_str[1] = NULL;

  aa0[0] = NULL;
  /* second initialiation - get commmand line arguments */
  initenv (argc, argv, &m_msg, &pst,&aa0[0]);

#ifdef COMP_THR
  /* now have max_workers - allocate work_info[] */
  if (max_workers >= MAX_WORKERS) max_workers = MAX_WORKERS;
  if ((work_info=
       (struct thr_str *)calloc(max_workers,sizeof(struct thr_str)))==NULL) {
    fprintf(stderr, " cannot allocate work_info[%d]\n",max_workers);
    exit(1);
  }
#else
  max_workers = 1;
#endif

#ifndef PRSS
  /* label library size limits */
  if (m_msg.n1_low > 0 && m_msg.n1_high < BIGNUM) 
    sprintf(lib_label,"library (range: %d-%d)",m_msg.n1_low,m_msg.n1_high);
  else if (m_msg.n1_low > 0) 
    sprintf(lib_label,"library (range: >%d)",m_msg.n1_low);
  else if (m_msg.n1_high < BIGNUM)
    sprintf(lib_label,"library (range: <%d)",m_msg.n1_high);
  else
    strncpy(lib_label,"library",sizeof(lib_label));
#else
  sprintf(lib_label,"shuffled sequence");
#endif
  lib_label[sizeof(lib_label)-1]='\0';

  tstart = tscan = s_time();
  tdstart = time(NULL);

  /* Allocate space for the query and library sequences */
  /* pad aa0[] with an extra 32 chars for ALTIVEC padding */
  if (aa0[0]==NULL) {
    if ((aa0[0] = (unsigned char *)malloc((m_msg.max_tot+1+32)*sizeof(unsigned char)))
	== NULL)
      s_abort ("Unable to allocate query sequence", "");
    *aa0[0]=0;
    aa0[0]++;
  }
  aa0[5]=aa0[4]=aa0[3]=aa0[2]=aa0[1]=aa0[0];

  /* make room for random sequence -
     also used as storage for COMP_THR library overlaps 
  */
  if ((aa1s = (unsigned char *)malloc((m_msg.max_tot+1+32)*sizeof (char))) == NULL) {
    s_abort ("Unable to allocate shuffled library sequence", "");
  }
  *aa1s=0;
  aa1s++;

  irand(0);

  if (m_msg.markx & MX_HTML) {
#ifdef HTML_HEAD    
    fprintf(stdout,"<html>\n<head>\n<title>%s Results</title>\n</head>\n<body>\n",prog_func);
#endif
    fprintf(stdout,"<pre>\n");
  }

#ifdef PGM_DOC
    fputs(argv_line,stdout);
    fputc('\n',stdout);
#endif  

  fprintf(stdout,"%s\n",iprompt0);
  fprintf(stdout," %s%s\n",verstr,refstr);
  if (m_msg.markx & MX_HTML) fputs("</pre>\n",stdout);

  /* Query library */
  if (m_msg.tname[0] == '\0') {
      if (m_msg.quiet == 1)
	s_abort("Query sequence undefined","");
    l1:	fputs (iprompt1, stdout);
      fflush  (stdout);
      if (fgets (m_msg.tname, MAX_FN, stdin) == NULL)
	s_abort ("Unable to read query library name","");
      m_msg.tname[MAX_FN-1]='\0';
      if ((bp=strchr(m_msg.tname,'\n'))!=NULL) *bp='\0';
      if (m_msg.tname[0] == '\0') goto l1;
  }
  
  /* Fetch first sequence */
  qoffset = 0l;
  qlcont = 0;
#ifdef COMP_MLIB
  /* Open query library */
  if ((q_file_p= openlib(m_msg.tname, m_msg.qdnaseq,qascii,!m_msg.quiet,NULL))==NULL) {
    s_abort(" cannot open library ",m_msg.tname);
  }
  qlib = 0;
  m_msg.n0 = 
    QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	     &qseek, &qlcont,q_file_p,&m_msg.sq0off);
  if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
  strncpy(qlabel,m_msg.qtitle,sizeof(qlabel));
  if (bp != NULL) *bp = ' ';
  qlabel[sizeof(qlabel)-1]='\0';

  /* if annotations are included in sequence, remove them */
  if (m_msg.ann_flg) {
    m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg,m_msg.qdnaseq);
  }

  if (m_msg.term_code && !(m_msg.qdnaseq==SEQT_DNA || m_msg.qdnaseq==SEQT_RNA) &&
      aa0[0][m_msg.n0-1]!='*') {
    aa0[0][m_msg.n0++]='*';
    aa0[0][m_msg.n0]=0;
  }

  /* check for subset */
  if (q_file_p->opt_text[0]!='\0') {
    if (q_file_p->opt_text[0]=='-') {
      sstart=0; sscanf(&q_file_p->opt_text[1],"%d",&sstop);
    }
    else {
      sscanf(&q_file_p->opt_text[0],"%d-%d",&sstart,&sstop);
      sstart--;
      if (sstop <= 0 ) sstop = BIGNUM;
    }
    for (id=0,is=sstart; is<min(m_msg.n0,sstop); ) aa0[0][id++]=aa0[0][is++];
    aa0[0][id]=0;
    m_msg.n0 = min(m_msg.n0,sstop)-sstart;
    if (m_msg.sq0off==1) m_msg.sq0off = sstart+1;
  }

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
  /* for ALTIVEC, must pad with 15 NULL's */
  for (id=0; id<SEQ_PAD; id++) {aa0[0][m_msg.n0+id]=0;}
#endif

  if (qlcont) {
    qoffset += m_msg.n0 - m_msg.sq0off;
  }
  else {
    qoffset = 0l;
  }

#else
  m_msg.n0 = getseq (m_msg.tname, qascii, aa0[0], m_msg.max_tot,
		     m_msg.qtitle, sizeof(m_msg.qtitle),
		     &m_msg.sq0off);
  strncpy(qlabel,m_msg.tname,sizeof(qlabel));
  qlabel[sizeof(qlabel)-1]='\0';

  /* if annotations are included in sequence, remove them */
  if (m_msg.ann_flg) {
    m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg,m_msg.qdnaseq);
  }
#endif

  if (m_msg.n0 > MAXTST) {
    fprintf(stderr," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    fprintf(stdout," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    aa0[0][MAXTST]='\0';
    m_msg.n0=MAXTST;
  }

  if (m_msg.qdnaseq == SEQT_UNK) {

  /* do automatic sequence recognition,but only for sequences > 20 residues */
    if (m_msg.n0 > 20 &&
	(float)scanseq(aa0[0],m_msg.n0,"ACGTUNacgtun")/(float)m_msg.n0 >0.85) {
      pascii = nascii;
      m_msg.qdnaseq = SEQT_DNA;
    }
    else {	/* its protein */
      pascii = aascii;
      m_msg.qdnaseq = SEQT_PROT;
    }
    /* modify qascii to use encoded version 
       cannot use memcpy() because it loses annotations 
    */
    re_ascii(qascii,pascii);
    init_ascii(pst.ext_sq_set,qascii,m_msg.qdnaseq);
    m_msg.n0 = recode(aa0[0],m_msg.n0,qascii, pst.nsqx);
  }

  if (m_msg.n0 <= 0)
    s_abort ("Query sequence length <= 0: ", m_msg.tname);

#ifdef SUPERFAMNUM
  m_msg.nqsfnum = nsfnum;
  for (i=0; i <= nsfnum & i<10; i++) m_msg.qsfnum[i] = sfnum[i];
  m_msg.nqsfnum_n = nsfnum_n;
  for (i=0; i <= nsfnum_n & i<10; i++) m_msg.qsfnum_n[i] = sfnum_n[i];
#endif

  resetp (&m_msg, &pst);

#ifndef COMP_MLIB
  gettitle(m_msg.tname,m_msg.qtitle,sizeof(m_msg.qtitle));
  if (m_msg.tname[0]=='-' || m_msg.tname[0]=='@') {
    strncmp(m_msg.tname,m_msg.qtitle,sizeof(m_msg.tname));
    if ((bp=strchr(m_msg.tname,' '))!=NULL) *bp='\0';
  }
#endif

  /* get library file names */

#ifndef PRSS
  if (strlen (m_msg.lname) == 0) {
    if (m_msg.quiet == 1) s_abort("Library name undefined","");
    libchoice(m_msg.lname,sizeof(m_msg.lname),&m_msg);
  }
  
  libselect(m_msg.lname, &m_msg);
#else
  if (strlen (m_msg.lname) == 0) {
    if (m_msg.quiet == 1) s_abort("Shuffle sequence undefined","");
l2:    fputs(iprompt2,stdout);
    fflush(stdout);
    if (fgets (m_msg.lname, MAX_FN, stdin) == NULL)
      s_abort ("Unable to read shuffle file name","");
    m_msg.lname[MAX_FN-1]='\0';
    if ((bp=strchr(m_msg.lname,'\n'))!=NULL) *bp='\0';
    if (m_msg.lname[0] == '\0') goto l2;
  }
  m_msg.lbnames[0]= m_msg.lname;
  m_msg.nln = 1;
  m_msg.nshow = 0;
#endif

  /* Get additional parameters here */
  if (!m_msg.quiet) query_parm (&m_msg, &pst);
  
  last_init(&m_msg, &pst);

  /* Allocate space for saved scores */
  if ((best = 
       (struct beststr *)calloc((MAXBEST+1),sizeof(struct beststr)))==NULL)
    s_abort ("Cannot allocate best struct","");
  if ((bestp_arr = 
       (struct beststr **)malloc((MAXBEST+1)*sizeof(struct beststr *)))==NULL)
    s_abort ("Cannot allocate bestp_arr","");
  
  /* Initialize bestp_arr */
  for (nbest = 0; nbest < MAXBEST+1; nbest++)
    bestp_arr[nbest] = &best[nbest];
  best++; bestp_arr++;
  best[-1].score[0]=best[-1].score[1]=best[-1].score[2]= INT_MAX;
  best[-1].zscore=FLT_MAX;	/* for Z-scores, bigger is best */
  best[-1].escore=FLT_MIN;	/* for E()-values, lower is best */

  if ((stats =
       (struct stat_str *)calloc(MAXSTATS,sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate stats struct","");

#ifdef UNIX
  /* set up signals now that input is done */
  signal(SIGHUP,SIG_IGN);
#endif

#ifdef COMP_THR
  /* Set up buffers for reading the library:

     We will start by using a 2 Mbyte buffer for each worker.  For
     proteins, that means 5,000 sequences of length 400 (average).
     For DNA, that means 2,000 sequences of length 1000.  At the
     moment, those are good averages.
  */

  if (m_msg.ldnaseq== SEQT_DNA) {
    max_buf_cnt = MAX_NT_BUF;
    ave_seq_len = AVE_NT_LEN;
  }
  else {
    max_buf_cnt = MAX_AA_BUF;
    ave_seq_len = AVE_AA_LEN;
  }

  /* however - buffer sizes should be a function of the number of
     workers so that all the workers are kept busy.  Assuming a 10,000
     entry library is the smallest we want to schedule, then
  */

  if (max_buf_cnt > 10000/max_workers) 
    max_buf_cnt = 10000/(2*max_workers);

  max_buf_cnt /= m_msg.thr_fact;

  /* finally, max_work_buf should be mod 6 for tfasta */
  max_buf_cnt -= (max_buf_cnt % 6);

  max_work_buf = 2*max_workers;

  /* allocate space for library buffers and results */

  buf_siz=max_buf_cnt*ave_seq_len;
  if (buf_siz < m_msg.max_tot) buf_siz = m_msg.max_tot;
  for (i=0; i<max_work_buf; i++) {
    if ((buf_list[i].buf =(struct buf_str *)calloc((size_t)(max_buf_cnt+1),
						   sizeof(struct buf_str)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",i,max_buf_cnt+1);
      exit(1);
    }
    buf_list[i].buf_cnt=0;
    buf_list[i].have_results=0;
    if ((buf_list[i].start =
         (unsigned char *)calloc((size_t)(buf_siz),sizeof(unsigned char)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer %d\n",i);
      exit(1);
    }

    /* make certain there is a '\0' at the beginning */
    buf_list[i].start++;

    reader_buf[i] = &buf_list[i];
  }

  /* initialization of global variables for threads/buffers */

  num_worker_bufs = 0;
  num_reader_bufs = max_work_buf;
  reader_done = 0;
  worker_buf_workp = 0;
  worker_buf_readp = 0;
  reader_buf_workp = 0;
  reader_buf_readp = 0;

  start_thread = 1;	/* keeps threads from starting */
#endif

  /* Label the output */
  if ((bp = (char *) strchr (m_msg.lname, ' ')) != NULL) *bp = '\0';
  if (m_msg.ltitle[0] == '\0') {
    strncpy(m_msg.ltitle,m_msg.lname,sizeof(m_msg.ltitle));
    m_msg.ltitle[sizeof(m_msg.ltitle)-1]='\0';
  }

#ifdef COMP_MLIB
  printf("Query library %s vs %s library\n", m_msg.tname,m_msg.lname);
  if (m_msg.nln > 0) printf("searching %s library\n\n",m_msg.lbnames[0]);
#endif

#ifdef COMP_MLIB
  while(1) {
    m_msg.db.length = 0l;
    m_msg.db.entries = m_msg.db.carry = 0;
    qlib++;
    stats_done = 0;
#endif

  maxl = m_msg.max_tot - m_msg.n0 -2;	/* maxn = max library sequence space */

  maxn = reset_maxn(&m_msg,maxl);
  pst.maxlen = maxn;

  outfd = stdout;  
  nbest = 0;
  zbestcut = -FLT_MAX;
  nstats = 0;

  /* get the last parameters */
  last_params(aa0[0],m_msg.n0, &m_msg, &pst);

  /*
     if our function returns approximate E()-scores, we do not need to
     work with raw scores and later calculate z-scores.  When
     approx. E()-scores are calculated, we still need various
     statistics structures, but we can get them immediately.  In this
     case, find_zp() must produce a z_score (large positive is good)
     from an e_score.
  */

  if (m_msg.escore_flg) {
    pst.zsflag_f = process_hist(stats,nstats,m_msg,pst,
				&m_msg.hist,&m_msg.pstat_void,0);
    stats_done=1;
  }

#ifndef COMP_THR
  if (m_msg.qshuffle) {
    if ((aa0s=(unsigned char *)calloc(m_msg.n0+2,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate aa0s[%d]\n",m_msg.n0+2);
      exit(1);
    }
    *aa0s='\0';
    aa0s++;
    memcpy(aa0s,aa0[0],m_msg.n0);
    qshuffle(aa0s,m_msg.n0,m_msg.nm0);
  }

  /* previous versions of FASTA have stored the reverse complement in
     the same array as the forward query sequence.  This version
     changes that, by allocating separate space for the reverse complement,
     and thus reducing the demand for a large MAXLIB/MAXTRN for long queries
  */
  if (m_msg.qframe == 2) {
    if ((aa0[1]=(unsigned char *)calloc(m_msg.n0+2,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate aa0[1][%d]\n",m_msg.n0+2);
      exit(1);
    }
    *aa0[1] = '\0';
    aa0[1]++;
    memcpy(aa0[1],aa0[0],m_msg.n0+1);
    revcomp(aa0[1],m_msg.n0,&pst.c_nt[0]);
  }
  /* set aa1 for serial - threaded points aa1 to buffer */

  aa1 = aa0[0] + m_msg.n0+1;	/* modified now that aa0[1] is done separately */
  *aa1++='\0';
#else
  init_thr(max_workers, work_info, m_msg, &pst, aa0[0], max_work_buf);
#endif

  if (m_msg.qshuffle && qstats==NULL) {
    if ((qstats =
	 (struct stat_str *)calloc(m_msg.shuff_max+1,sizeof(struct stat_str)))==NULL)
      s_abort ("Cannot allocate qstats struct","");
  }
  nqstats = 0;

  if (m_msg.markx & MX_HTML) fputs("<pre>\n",stdout);
#ifndef PRSS
  /* rline[] is a tmp string */
  if (m_msg.qdnaseq == SEQT_DNA || m_msg.qdnaseq == SEQT_RNA) {
    strncpy(rline,(m_msg.qframe==1)? " (forward-only)" : "\0",sizeof(rline));
    rline[sizeof(rline)-1]='\0';
  }
  else rline[0]='\0';

  leng = (int)strlen(m_msg.qtitle);
  if (leng > 50) leng -= 10;

  sprintf (&m_msg.qtitle[leng], " %d %s", m_msg.n0, m_msg.sqnam);
  m_msg.seqnm = 0;


#ifdef COMP_MLIB
  printf("%3d>>>%s - %d %s%s\n vs  %.60s %s\n", qlib,
	 m_msg.qtitle, m_msg.n0, m_msg.sqnam,
	 (m_msg.revcomp ? " (reverse complement)" : rline),
	 m_msg.ltitle,lib_label);
#else
  printf("%.50s: %d %s%s\n %s\n vs  %.60s %s\n",
	 qlabel, m_msg.n0, m_msg.sqnam,
	 (m_msg.revcomp ? " (reverse complement)" : rline),
	 m_msg.qtitle,m_msg.ltitle,lib_label);
#endif
  libstr_p = &libstr[0];
  n_libstr=sizeof(libstr);
#else		/* PRSS */
  libstr_p = &m_msg.ltitle[0];
  n_libstr= sizeof(m_msg.ltitle);
  set_shuffle(m_msg);	/* set count/width parameters in llgetaa.c */ 
#endif

  fflush (outfd);

  tprev = s_time();
  
  if (m_msg.dfile[0] && (fdata=fopen(m_msg.dfile,"w"))!=NULL)
    fprintf(fdata,"%3d\t%-50s\n",m_msg.n0,m_msg.qtitle);

  qtt.length += m_msg.n0;
  qtt.entries++;

#ifdef COMP_THR
  start_thr();

  /* now open the library and start reading */
  /* get a buffer and fill it up */
  get_rbuf(&cur_buf,max_work_buf);

  cur_buf->buf_cnt = 0;
  cur_buf->have_results = 0;
  cur_buf->buf[0].aa1b = cur_buf->start;
  ntbuff = 0;
  nseq = 0;
#else /* ! COMP_THR */
  /* initialize the comparison function, returning f_str */
  init_work (aa0[0], m_msg.n0, &pst, &f_str[0]);
  have_f_str=1;

  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0];
  if (m_msg.qframe == 2) {
    init_work ( aa0[1], m_msg.n0, &pst, &f_str[1]);
  }
  if (m_msg.qshuffle) {
    init_work ( aa0s, m_msg.n0, &pst, &qf_str);
  }
#endif	/* COMP_THR */

  /* open the library - start the search */

  for (iln = 0; iln < m_msg.nln; iln++) {
    if ((m_msg.lb_mfd[iln] = m_file_p=
	 openlib(m_msg.lbnames[iln], m_msg.ldnaseq, lascii, !m_msg.quiet, m_msg.lb_mfd[iln]))
	==NULL) {
      fprintf(stderr," cannot open library %s\n",m_msg.lbnames[iln]);
      continue;
    }
#if !defined(PRSS) && !defined(COMP_MLIB)
    else
      printf ("searching %s %s\n",m_msg.lbnames[iln],lib_label);
#endif

    loffset = 0l;
    lcont = 0;
    ocont = 0;
    n1tot_v = n1tot_cnt = 0;
    n1tot_cur = n1tot_ptr = NULL;

    /* get next buffer to read into */
    maxt = maxn;

#ifndef COMP_THR
    aa1ptr = aa1;
#else
    /* read sequence directly into buffer */
    aa1ptr = aa1 = cur_buf->buf[nseq].aa1b;
#endif

    while ((n1=GETLIB(aa1ptr,maxt,libstr_p,n_libstr,&lmark,&lcont,m_file_p,&l_off))>=0) {

      if (n_libstr <= MAX_UID) {
	if ((bp=strchr(libstr_p,' '))!=NULL) *bp='\0';
      }

      if (m_msg.term_code && !lcont &&
	  m_msg.ldnaseq==SEQT_PROT && aa1ptr[n1-1]!=m_msg.term_code) {
	aa1ptr[n1++]=m_msg.term_code;
	aa1ptr[n1]=0;
      }

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
      /* for ALTIVEC, must pad with 15 NULL's */
      for (id=0; id<SEQ_PAD; id++) {aa1ptr[n1+id]=0;}
#endif

#ifdef DEBUG
      if (aa1[-1]!='\0' || aa1ptr[n1]!='\0') {
	fprintf(stderr,"%s: aa1[%d] missing NULL boundaries: %d %d\n",libstr_p,n1,aa1[-1],aa1ptr[n1]);
      }
#endif

      /* check for a continued sequence and provide a pointer to 
	 the n1_tot array if lcont || ocont */
      n1tot_v += n1;
      if (lcont && !ocont) {	/* get a new pointer */
	if (n1tot_cnt <= 0) {
	  if ((n1tot_ptr=calloc(1000,sizeof(int)))==NULL) {
	    fprintf(stderr," cannot allocate n1tot_ptr\n");
	    exit(1);
	  }
	  else {n1tot_cnt=1000;}
	}
	n1tot_cnt--;
	n1tot_cur = n1tot_ptr++;
      }

      if (n1tot_v < m_msg.n1_low || n1tot_v > m_msg.n1_high) {
	goto loop2;
      }

      m_msg.db.entries++;
      m_msg.db.length += n1;
      if (m_msg.db.length > LONG_MAX) {
	m_msg.db.length -= LONG_MAX; m_msg.db.carry++;
      }

#ifdef DEBUG
      /* This finds most reasons for core dumps */
      if (pst.debug_lib)
	for (i=0; i<n1; i++)
	  if (aa1[i]>=pst.nsqx) 
	      {fprintf(stderr,
		       "%s residue[%d/%d] %d range (%d) lcont/ocont: %d/%d\n%s\n",
		       libstr,i,n1,aa1[i],pst.nsqx,lcont,ocont,aa1ptr+i);
	      aa1[i]=0;
	      n1=i-1;
	      break;
	      }
#endif

      /* don't count long sequences more than once */
      if (aa1!=aa1ptr) {n1 += m_msg.loff; m_msg.db.entries--;}

#ifdef PROGRESS
      if (!m_msg.quiet) 
	if (m_msg.db.entries % 200 == 199) {
	  fputc('.',stderr);
	  if (m_msg.db.entries % 10000 == 9999) fputc('\n',stderr);
	  else if (m_msg.db.entries % 1000 == 999) fputc(' ',stderr);

	}
#endif

      if (n1<=1) {
	/*	if (igncnt++ <10)
		fprintf(stderr,"Ignoring: %s\n",libstr);
	*/
	goto loop2;
      }

#ifdef PRSS
      if (lmark==0) {
	n1s = n1;
	memcpy(aa1s,aa1,n1s);
	m_msg.db.entries=0;
	m_msg.db.length=0;
      }
#endif

      /* if COMP_THR - fill and empty buffers */
#ifdef COMP_THR
      ntbuff += n1+1;

      for (itt=m_msg.revcomp; itt<=m_msg.nitt1; itt++) {

	cur_buf->buf_cnt++;
	cur_buf_p = &(cur_buf->buf[nseq++]);
	cur_buf_p->n1  = n1;
	cur_buf_p->n1tot_p = n1tot_cur;
	cur_buf_p->lseek = lmark;
	cur_buf_p->cont = ocont+1;
	cur_buf_p->m_file_p = (void *)m_file_p;
	cur_buf_p->frame = itt;
	memcpy(cur_buf_p->libstr,libstr,MAX_UID);
#ifdef SUPERFAMNUM
	cur_buf_p->nsfnum = nsfnum;
	if ((cur_buf_p->sfnum[0]=sfnum[0])>0 &&
	    (cur_buf_p->sfnum[1]=sfnum[1])>0 &&
	    (cur_buf_p->sfnum[2]=sfnum[2])>0 &&
	    (cur_buf_p->sfnum[3]=sfnum[3])>0 &&
	    (cur_buf_p->sfnum[4]=sfnum[4])>0 &&
	    (cur_buf_p->sfnum[5]=sfnum[5])>0 &&
	    (cur_buf_p->sfnum[6]=sfnum[6])>0 &&
	    (cur_buf_p->sfnum[7]=sfnum[7])>0 &&
	    (cur_buf_p->sfnum[8]=sfnum[8])>0 &&
	    (cur_buf_p->sfnum[9]=sfnum[9])>0) ;
#endif

	/* this assumes that max_buf_cnt is guaranteed %6=0 so that
	   additional pointers to the same buffer can be used 
	   nseq now points to next buffer
	*/

	cur_buf->buf[nseq].aa1b = cur_buf->buf[nseq-1].aa1b;
      } /* for (itt .. */

      /* make a copy of the overlap (threaded only) */
      if (lcont) {
	memcpy(aa1s,&aa1[n1-m_msg.loff],m_msg.loff);
      }

      /* if the buffer is filled */
      if (nseq >= max_buf_cnt || ntbuff >= buf_siz - maxn) {

	/* provide filled buffer to workers */
	put_rbuf(cur_buf,max_work_buf);

	/* get an empty buffer to fill */
	get_rbuf(&cur_buf,max_work_buf);

	/* "empty" buffers have results that must be processed */
	if (cur_buf->buf_cnt && cur_buf->have_results) {
	  save_best(cur_buf,m_msg,pst,fdata,m_msg.qsfnum,&m_msg.hist,
		    &m_msg.pstat_void
#ifdef PRSS
		    ,s_score,&s_n1
#endif
		    );

	}

	/* now the buffer is truly empty, fill it up */
	cur_buf->buf_cnt = 0;
	cur_buf->have_results = 0;
	/* point the first aa1 ptr to the buffer start */
	aa1=cur_buf->buf[0].aa1b = cur_buf->start;
	ntbuff = 0;
	nseq=0;
      }
      else {	/* room left in current buffer, increment ptrs */
	aa1=cur_buf->buf[nseq].aa1b = cur_buf->buf[nseq-1].aa1b+n1+1;
      }
#else /* if !COMP_THR - do a bunch of searches */

      /* t_best and t_rbest are used to save the best score or shuffled
	 score from all the frames */

      t_best = t_rbest = t_qrbest = -1;
      t_escore = t_rescore = t_qrescore = FLT_MAX;
      for (itt=m_msg.revcomp; itt<=m_msg.nitt1; itt++) {

	rst.score[0] = rst.score[1] = rst.score[2] = 0;
	do_work (aa0[itt], m_msg.n0,aa1,n1,itt,&pst,f_str[itt],0,&rst);

	if (rst.score[pst.score_ix] > t_best) {
	  t_best = rst.score[pst.score_ix];
	}

	if (fdata) {
	  fprintf(fdata,
		  "%-12s %5d %6d %d %.5f %.5f %4d %4d %4d %g %d %d %8lld\n",
		  libstr,
#ifdef SUPERFAMNUM
		  sfn_cmp(m_msg.qsfnum,sfnum),
#else
		  0,
#endif
		  n1,itt,
		  rst.comp,rst.H,
		  rst.score[0],rst.score[1],rst.score[2],
		  rst.escore, rst.segnum, rst.seglen, lmark);
	  fflush(fdata);
	}

#ifdef PRSS
	if (lmark==0) {
	  s_score[0] = rst.score[0];
	  s_score[1] = rst.score[1];
	  s_score[2] = rst.score[2];

	  s_n1 = n1;
	  aa1_loff = l_off;
	}
	t_best = t_rbest = rst.score[pst.score_ix];
	t_escore = t_rescore = rst.escore;
#else
	if (m_msg.qshuffle) {
	  do_work (aa0s, m_msg.n0,aa1,n1,itt,&pst,qf_str,1,&rrst);

	  if (rrst.score[pst.score_ix] > t_qrbest) 
	    t_qrbest = rrst.score[pst.score_ix];
	  if (rrst.escore < t_qrescore) 
	    t_qrescore = rrst.escore;

	  if (itt==m_msg.nitt1 && nqstats < m_msg.shuff_max) {
	    qstats[nqstats].n1 = n1;	/* save the best score */
	    qstats[nqstats].comp =  rst.comp;
	    qstats[nqstats].H = rst.H;
	    qstats[nqstats].escore = t_qrescore;
	    qstats[nqstats++].score = t_qrbest;
	    t_qrbest = -1;	/* reset t_qrbest, t_qrescore */
	    t_qrescore = FLT_MAX;
	  }
	}

	if (pst.zsflag >= 10) {
	  if (pst.zs_win > 0) wshuffle(aa1,aa1s,n1,pst.zs_win,&ieven);
	  else shuffle(aa1,aa1s,n1);
	  do_work (aa0[itt], m_msg.n0, aa1s, n1,itt,&pst,f_str[itt],0,&rrst);
	  if (rrst.score[pst.score_ix] > t_rbest) {
	    t_rbest = rrst.score[pst.score_ix];
	    t_rescore = rrst.escore;
	  }
	}
#endif
	i_score = rst.score[pst.score_ix];

/* this section saves scores for statistics calculations.  For
   comparisons that can be from one of 2 or 6 frames, it should only
   be run once, for the best of the 2 or 6 scores.  t_rbest,t_rescore
   have the best of the 2 or 6 scores from the frames.  For proteins,
   this is run for every score.

*/   
#ifdef PRSS	/* don't save the first score (unshuffled) with PRSS */
	if (lmark > 0) {
#endif   	

	if (itt == m_msg.nitt1) {
	  if (nstats < MAXSTATS) {
	    stats[nstats].n1 = n1;	/* save the best score */
	    stats[nstats].comp =  rst.comp;
	    stats[nstats].H = rst.H;
	    if (pst.zsflag >=10) {
	      t_best = t_rbest;
	      t_escore = t_rescore;
	    }
	    stats[nstats].escore = t_escore;
	    stats[nstats++].score = t_best;
	    t_best = t_rbest = -1;	/* reset t_rbest, t_best */
	    t_escore = t_rescore = FLT_MAX;
	  }
	  else if (pst.zsflag >= 0) {
	    if (!stats_done) {
	      pst.zsflag_f = process_hist(stats,nstats,m_msg,pst,
					  &m_msg.hist,&m_msg.pstat_void,0);
	      stats_done = 1;
	      kstats = nstats;
	      for (i=0; i<MAXBEST; i++) {
		bestp_arr[i]->zscore = 
		  (*find_zp)(bestp_arr[i]->score[pst.score_ix],
			     bestp_arr[i]->escore, bestp_arr[i]->n1,
			     bestp_arr[i]->comp, m_msg.pstat_void);
	      }
	      zbestcut = bestp_arr[nbest-1]->zscore;
	    }

#ifdef SAMP_STATS
/* older versions saved the first MAXSTATS scores, and ignored the
   rest in the statistics.  With SAMP_STATS, scores after MAX_STATS
   are sampled at random, and included in the sample set and the
   statistics parameters are re-derived at the end of the run using
   the sampled scores.

   It would be faster not to do the nrand(); if(jstats < MAXSTATS)
   less often.
*/
	    if (!m_msg.escore_flg) {	/* only for zscores */
	      jstats = nrand(++kstats); /* no mod % 0 */
	      if (jstats < MAXSTATS) {
		stats[jstats].n1 = n1;	/* save the best score */
		stats[jstats].comp =  rst.comp;
		stats[jstats].H = rst.H;
		if (pst.zsflag >=10) t_best = t_rbest;
		stats[jstats].score = t_best;
	      }
	    }
#endif
	  }	/* ( nstats >= MAXSTATS) && zsflag >= 0 */
	}	/* itt1 == nitt1 */
#ifdef PRSS
	}
#endif

	/* this section completes work on the current score */
	if (stats_done) { /* stats_done > 0 => nstats >= MAXSTATS */
	  zscore=(*find_zp)(i_score, rst.escore, n1, rst.comp,
			    m_msg.pstat_void);

	  if (itt == m_msg.nitt1) {
	    if (pst.zsflag >= 10) t_best = t_rbest;
	    
	    addhistz((*find_zp)(t_best, t_escore, n1, rst.comp, 
				m_msg.pstat_void),
		     &m_msg.hist);
	    t_best = t_rbest = -1;
	  }
	}
	else zscore = (double) i_score;

#ifndef PRSS
	if (zscore > zbestcut ) {
	  if (nbest >= MAXBEST) {
	    bestfull = nbest-MAXBEST/4;
	    selectbestz(bestp_arr,bestfull-1,nbest);
	    zbestcut = bestp_arr[bestfull-1]->zscore;
	    nbest = bestfull;
	  }

	  bestp = bestp_arr[nbest++];
	  bestp->score[0] = rst.score[0];
	  bestp->score[1] = rst.score[1];
	  bestp->score[2] = rst.score[2];
	  bestp->comp =  rst.comp;
	  bestp->H = rst.H;
	  bestp->zscore = zscore;
	  bestp->escore = rst.escore;
	  bestp->segnum = rst.segnum;
	  bestp->seglen = rst.seglen;
	  bestp->lseek = lmark;
	  bestp->cont = ocont+1;
	  bestp->m_file_p = m_file_p;
	  bestp->n1 = n1;
	  bestp->n1tot_p=n1tot_cur;
	  bestp->frame = itt;
	  memcpy(bestp->libstr,libstr,MAX_UID);
#ifdef SUPERFAMNUM
	  bestp->nsfnum = nsfnum;
	  if ((bestp->sfnum[0]=sfnum[0])>0 &&
	      (bestp->sfnum[1]=sfnum[1])>0 &&
	      (bestp->sfnum[2]=sfnum[2])>0 &&
	      (bestp->sfnum[3]=sfnum[3])>0 &&
	      (bestp->sfnum[4]=sfnum[4])>0 &&
	      (bestp->sfnum[5]=sfnum[5])>0 &&
	      (bestp->sfnum[6]=sfnum[6])>0 &&
	      (bestp->sfnum[7]=sfnum[7])>0 &&
	      (bestp->sfnum[8]=sfnum[8])>0 &&
	      (bestp->sfnum[9]=sfnum[9])>0) ;
#endif
	}
#else 	/* PRSS */
	if (lmark == 0) {
	  bestp = bestp_arr[nbest++];
	  bestp->score[0] = rst.score[0];
	  bestp->score[1] = rst.score[1];
	  bestp->score[2] = rst.score[2];
	  bestp->comp =  rst.comp;
	  bestp->H = rst.H;
	  bestp->zscore = zscore;
	  bestp->escore = rst.escore;
	  bestp->segnum = rst.segnum;
	  bestp->seglen = rst.seglen;
	  bestp->lseek = lmark;
	  bestp->cont = 0;
	  bestp->m_file_p = m_file_p;
	  bestp->n1 = n1;
	  bestp->n1tot_p=n1tot_cur;
	  bestp->frame = itt;
	  memcpy(bestp->libstr,libstr,MAX_UID);
	  bestp->nsfnum = 0;
	}
#endif
      }
#endif

    loop2: 
      if (lcont) {
	maxt = m_msg.maxt3;
#ifndef COMP_THR
	memcpy(aa1,&aa1[n1-m_msg.loff],m_msg.loff);
#else
	memcpy(aa1,aa1s,m_msg.loff);
#endif
	aa1ptr= &aa1[m_msg.loff];
	loffset += n1 - m_msg.loff;
	ocont = lcont;
      }
      else {
	maxt = maxn;
	aa1ptr=aa1;
	if (ocont) *n1tot_cur = n1tot_v;
	ocont = 0;
	loffset = 0l;
	n1tot_v = 0;
	n1tot_cur = NULL;
      }
    } /* end while((n1=getlib())) */
  } /* end iln=1..nln */

  /* all done */

#ifdef COMP_THR
  /* check last buffers for any results */
  put_rbuf_done(max_workers,cur_buf,max_work_buf);

  for (i=0; i < num_reader_bufs; i++) {
    reader_buf_readp = (reader_buf_readp+1)%(max_work_buf);
    if (reader_buf[reader_buf_readp]->buf_cnt > 0 && 
	reader_buf[reader_buf_readp]->have_results) {
	  save_best(reader_buf[reader_buf_readp],m_msg,pst,fdata,m_msg.qsfnum,
		    &m_msg.hist, &m_msg.pstat_void
#ifdef PRSS
		    ,s_score,&s_n1
#endif
		    );
    }
  }
#endif

#ifdef PROGRESS
  if (!m_msg.quiet)
    if (m_msg.db.entries >= 200) {fprintf(stderr," Done!\n");}
#endif

  m_msg.nbr_seq = m_msg.db.entries;
  get_param(&pst, gstring2,gstring3);

/* *************************** */
/* analyze the last results    */
/* *************************** */
    
#ifndef PRSS
#ifndef SAMP_STATS
  if (!stats_done && nstats > 0) {
#endif
    pst.zsflag_f = process_hist(stats,nstats,m_msg,pst,&m_msg.hist,
				&m_msg.pstat_void,stats_done);
    if (m_msg.pstat_void != NULL) {
      stats_done = 1;
      for (i = 0; i < nbest; i++) {
	bestp_arr[i]->zscore =
	  (*find_zp)(bestp_arr[i]->score[pst.score_ix],
		     bestp_arr[i]->escore, bestp_arr[i]->n1, 
		     bestp_arr[i]->comp, m_msg.pstat_void);
      }
#ifndef SAMP_STATS
    }
    else pst.zsflag = -1;
#endif
  }
#else	/* PRSS */
  if (pst.zsflag < 10) pst.zsflag += 10;
  pst.zsflag_f = process_hist(stats,nstats,m_msg,pst,
			      &m_msg.hist, &m_msg.pstat_void,0);
  stats_done = 1;
  for (i = 0; i < nbest; i++) {
    bestp_arr[i]->zscore = (*find_zp)(bestp_arr[i]->score[pst.score_ix],
				      bestp_arr[i]->escore, bestp_arr[i]->n1,
				      bestp_arr[i]->comp, m_msg.pstat_void);
  }
#endif

  if (pst.zdb_size <= 1) pst.zdb_size = m_msg.db.entries;

#ifdef COMP_THR
  /* before I call last_calc/showbest/showalign, I need init_work() to
     get an f_str. This duplicates some code above, which is used in
     the non-threaded version
  */

  if (!have_f_str) {
    init_work(aa0[0],m_msg.n0,&pst,&f_str[0]);
    have_f_str = 1;
    f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] =  f_str[0];

    if (m_msg.qframe == 2) {
      if ((aa0[1]=(unsigned char *)calloc((size_t)m_msg.n0+2,
					  sizeof(unsigned char)))==NULL) {
	fprintf(stderr," cannot allocate aa0[1][%d] for alignments\n",
		m_msg.n0+2);
      }
      *aa0[1]='\0';
      aa0[1]++;
      memcpy(aa0[1],aa0[0],m_msg.n0+1);
      revcomp(aa0[1],m_msg.n0,&pst.c_nt[0]);
      init_work(aa0[1],m_msg.n0,&pst,&f_str[1]);
    }

    /* I also need a "real" aa1 */
    aa1 = buf_list[0].start;
#ifdef PRSS
    /* for PRSS - I need the original second (non-shuffled) sequence */
    memcpy(aa1,aa1s,n1s+1);
#endif
    }
#endif

/* now we have one set of scaled scores for in bestp_arr  -
   for FASTS/F, we need to do some additional processing */

  if (!m_msg.qshuffle) {
    last_stats(aa0[0], m_msg.n0, stats,nstats, bestp_arr,nbest,
	       m_msg, pst, &m_msg.hist, &m_msg.pstat_void);
  }
  else {
    last_stats(aa0[0], m_msg.n0,
	       qstats,nqstats, bestp_arr,nbest, m_msg, pst, 
	       &m_msg.hist, &m_msg.pstat_void);
  }

  /* here is a contradiction: if pst.zsflag < 0, then m_msg.pstat_void
     should be NULL; if it is not, then process_hist() has been called */
  if (pst.zsflag < 0 && m_msg.pstat_void != NULL) pst.zsflag = 1;

  if (m_msg.last_calc_flg) {
    /* last_calc may need coefficients from last_stats() */
    nbest = last_calc(aa0, aa1, maxn, bestp_arr, nbest, m_msg, &pst,
		      f_str, m_msg.pstat_void);
  }

  scale_scores(bestp_arr,nbest,m_msg.db,pst,m_msg.pstat_void);

  get_param(&pst, gstring2,gstring3);

#ifdef PRSS
  /*   gettitle(m_msg.lname,m_msg.ltitle,sizeof(m_msg.ltitle)); */
  printf("%.50s - %s %d %s%s\n vs %.60s - %s shuffled sequence\n",
	 m_msg.tname, m_msg.qtitle,m_msg.n0, m_msg.sqnam,
	 (m_msg.revcomp ? " (reverse complement)" : "\0"),
	 m_msg.lname,m_msg.ltitle);
#endif

  prhist (stdout, m_msg, pst, m_msg.hist, nstats, m_msg.db, gstring2);

  tscan = s_time();
  printf (" Scan time: ");
  ptime(stdout,tscan-tprev);
  printf ("\n");
#ifdef COMP_MLIB
  ttscan += tscan-tprev;
#endif

 l3:
  if (!m_msg.quiet) {
    printf("Enter filename for results [%s]: ", m_msg.outfile);
    fflush(stdout);
  }

  rline[0]='\0';
  if (!m_msg.quiet && fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
  if ((bp=strchr(rline,'\n'))!=NULL) *bp = '\0';
  if (rline[0]!='\0') strncpy(m_msg.outfile,rline,sizeof(m_msg.outfile));
  if (m_msg.outfile[0]!='\0') {
    if ((outfd=fopen(m_msg.outfile,"w"))==NULL) {
      fprintf(stderr," could not open %s\n",m_msg.outfile);
      if (!m_msg.quiet) goto l3;
      else goto l4;
    }

#ifdef PGM_DOC
    fputs(argv_line,outfd);
    fputc('\n',outfd);
#endif  
    fputs(iprompt0,outfd);
    fprintf(outfd," %s%s\n",verstr,refstr);

    fprintf(outfd," %s%s, %d %s\n vs %s %s\n",
	    qlabel, (m_msg.revcomp ? "-" : "\0"), m_msg.n0,
	    m_msg.sqnam, m_msg.ltitle, lib_label);

    prhist(outfd,m_msg,pst,m_msg.hist, nstats, m_msg.db, gstring2);
  }

 l4:   
  if (m_msg.markx & MX_HTML) {
      fputs("</pre>\n<p>\n<hr>\n<p>\n",outfd);
  }

  /* code from p2_complib.c to pre-calculate -m 9 alignment info -
     requires -q with -m 9 */

  if (m_msg.quiet || m_msg.markx & MX_M9SUMM) {

    /* to determine how many sequences to re-align (either for
       do_opt() or calc_id() we need to modify m_msg.mshow to get
       the correct number of alignments */

    if (m_msg.mshow_flg != 1 && pst.zsflag >= 0) {
      for (i=0; i<nbest && bestp_arr[i]->escore< m_msg.e_cut; i++) {}
      m_msg.mshow = i;
    }

#ifndef PRSS
    if (m_msg.mshow <= 0) { /* no results to display */
      fprintf(outfd,"!! No sequences with E() < %f\n",m_msg.e_cut);
      m_msg.nshow = 0;
      goto end_l;
    }
#endif
  }

#ifdef PRSS
  memcpy(aa1,aa1s,n1s);
  maxn = n1s;
  nbest = 1;
#endif

  showbest (stdout, aa0, aa1, maxn, bestp_arr, nbest, qtt.entries, &m_msg, pst,
	    m_msg.db, gstring2, f_str);

  if (outfd != stdout) {
    t_quiet = m_msg.quiet;
    m_msg.quiet = -1;	/* should guarantee 1..nbest shown */
    showbest (outfd, aa0, aa1, maxn, bestp_arr, nbest, qtt.entries, &m_msg, pst,
	      m_msg.db, gstring2, f_str);
    m_msg.quiet = t_quiet;
  }

  if (m_msg.nshow > 0) {
    rline[0]='N';
    if (!m_msg.quiet){
      printf(" Display alignments also? (y/n) [n] "); fflush(stdout);
      if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
    }
    else rline[0]='Y';

    if (toupper((int)rline[0])=='Y') {
      if (!m_msg.quiet) {
	printf(" number of alignments [%d]? ",m_msg.nshow);
	fflush(stdout);
	if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
	if (rline[0]!=0) sscanf(rline,"%d",&m_msg.nshow);
	m_msg.ashow=m_msg.nshow;
      }

      if (m_msg.markx & (MX_AMAP+ MX_HTML + MX_M9SUMM)) {
	fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		qlabel,(m_msg.revcomp ? "_rev":"\0"), m_msg.n0,
		m_msg.sqnam,m_msg.lname);
      }

      if (m_msg.markx & MX_M10FORM) {
	fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		qlabel,(m_msg.revcomp ? "-":"\0"), m_msg.n0, m_msg.sqnam,
		m_msg.lname);
	fprintf(outfd,"; pg_name: %s\n",argv[0]);
	fprintf(outfd,"; pg_ver: %s\n",mp_verstr);
	fprintf(outfd,"; pg_argv:");
	for (i=0; i<argc; i++)
	  fprintf(outfd," %s",argv[i]);
	fputc('\n',outfd);
	fputs(gstring3,outfd);
	fputs(hstring1,outfd);
      }

#ifndef PRSS
      showalign (outfd, aa0, aa1, maxn, bestp_arr, nbest, qtt.entries,
		 m_msg, pst, gstring2, f_str);
#else
      if (pst.sw_flag > 0 || (!m_msg.quiet && m_msg.nshow>0)) {
	showalign (outfd, aa0, aa1, maxn, bestp_arr, nbest, qtt.entries,
		 m_msg, pst, gstring2, f_str);
      }
#endif

      fflush(outfd);
    }
  }

  end_l:
#if defined(COMP_THR) && defined(COMP_MLIB)
    for (i=0; i<max_work_buf; i++) {
      buf_list[i].buf_cnt=0;
      buf_list[i].have_results=0;
    }

    num_worker_bufs = 0;
    num_reader_bufs = max_work_buf;
    reader_done = 0;
    worker_buf_workp = 0;
    worker_buf_readp = 0;
    reader_buf_workp = 0;
    reader_buf_readp = 0;

    start_thread = 1;	/* stop thread from starting again */
#endif

    /* clean up alignment encodings */
    for (i=0; i < m_msg.nshow; i++) {
      if (bestp_arr[i]->have_ares) {
	free(bestp_arr[i]->a_res.res);
	bestp_arr[i]->a_res.res = NULL;
	bestp_arr[i]->have_ares = 0;
      }
    }

    if (m_msg.qframe == 2) free(aa0[1]-1);

    if (have_f_str) {
      if (f_str[1]!=f_str[0]) {
	close_work (aa0[1], m_msg.n0, &pst, &f_str[1]);
      }
      close_work (aa0[0], m_msg.n0, &pst, &f_str[0]);
      have_f_str = 0;
#ifndef COMP_THR
      if (m_msg.qshuffle) close_work (aa0s, m_msg.n0, &pst, &qf_str);
#endif
      if (pst.pam_pssm) {
	free_pam2p(pst.pam2p[0]);
	free_pam2p(pst.pam2p[1]);
      }
    }

    for (iln=0; iln < m_msg.nln; iln++) {
      if (m_msg.lb_mfd[iln]!=NULL) closelib(m_msg.lb_mfd[iln]);
    }

    tddone = time(NULL);
    tdone = s_time();
    fflush(outfd);

    if (fdata) {
      fprintf(fdata,"/** %s **/\n",gstring2);
      fprintf(fdata,"%3ld%-50s\n",qtt.entries-1,m_msg.qtitle);
      fflush(fdata);
    }
    
#ifdef COMP_MLIB
    ttdisp += tdone-tscan;

    maxn = m_msg.max_tot;
    m_msg.n0 = 
      QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	       &qseek, &qlcont,q_file_p,&m_msg.sq0off);
    if (m_msg.n0 <= 0) break;
    if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
    strncpy(qlabel, m_msg.qtitle,sizeof(qlabel));
    if (bp != NULL) *bp=' ';
    qlabel[sizeof(qlabel)-1]='\0';

    if (m_msg.ann_flg) {
      m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg,m_msg.qdnaseq);
    }

    if (m_msg.term_code && m_msg.qdnaseq==SEQT_PROT &&
	aa0[0][m_msg.n0-1]!=m_msg.term_code) {
      aa0[0][m_msg.n0++]=m_msg.term_code;
      aa0[0][m_msg.n0]=0;
    }

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
    /* for ALTIVEC, must pad with 15 NULL's */
    for (id=0; id<SEQ_PAD; id++) {aa0[0][m_msg.n0+id]=0;}
#endif

#ifdef SUPERFAMNUM
    m_msg.nqsfnum = nsfnum;
    for (i=0; i <= nsfnum & i<10; i++) m_msg.qsfnum[i] = sfnum[i];
    m_msg.nqsfnum_n = nsfnum_n;
    for (i=0; i <= nsfnum_n & i<10; i++) m_msg.qsfnum_n[i] = sfnum_n[i];
#endif
  }
#endif
  if (m_msg.markx & MX_M10FORM)
      fprintf(outfd,">>><<<\n");

    tdone = s_time();
    if ( m_msg.markx & MX_HTML) fputs("<p><pre>\n",outfd); 
    printsum(outfd, m_msg.db);
    if ( m_msg.markx & MX_HTML) fputs("</pre>\n",outfd);
#ifdef HTML_HEAD
    if (m_msg.markx & MX_HTML) fprintf(outfd,"</body>\n</html>\n");
#endif
    if (outfd!=stdout) printsum(stdout,m_msg.db);

    exit(0);
}   /* End of main program */

void
printsum(FILE *fd, struct db_str ntt)
{
  double db_tt;
  char tstr1[26], tstr2[26];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */
  fprintf(fd, "\n\n%ld residues in %ld query   sequences\n", qtt.length, qtt.entries);
  if (ntt.carry == 0) 
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt.length, ntt.entries);
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.0f residues in %ld library sequences\n", db_tt, ntt.entries);
  }

#ifndef COMP_THR
  fprintf(fd," Scomplib [%s]\n start: %s done: %s\n",mp_verstr,tstr1,tstr2);
#else
  fprintf(fd," Tcomplib [%s] (%d proc)\n start: %s done: %s\n", mp_verstr,
    max_workers,tstr1,tstr2);
#endif
#ifndef COMP_MLIB
  fprintf(fd," Scan time: ");
  ptime(fd, tscan - tprev);
  fprintf (fd," Display time: ");
  ptime (fd, tdone - tscan);
#else
  fprintf(fd," Total Scan time: ");
  ptime(fd, ttscan);
  fprintf (fd," Total Display time: ");
  ptime (fd, ttdisp);
#endif
  fprintf (fd,"\n");
  fprintf (fd, "\nFunction used was %s [%s]\n", prog_func,verstr);
}

void fsigint()
{
  struct db_str db;

  db.entries = db.length = db.carry = 0;
  tdone = s_time();
  tddone = time(NULL);

  printf(" /*** interrupted ***/\n");
  if (outfd!=stdout) fprintf(outfd,"/*** interrupted ***/\n");
  fprintf(stderr,"/*** interrupted ***/\n");

  printsum(stdout,db);
  if (outfd!=stdout) printsum(outfd,db);

  exit(1);
}

#ifdef COMP_THR
void save_best(struct buf_head *cur_buf, struct mngmsg m_msg, struct pstruct pst, 
	       FILE *fdata, int *qsfnum, struct hist_str *histp,
	       void **pstat_voidp
#ifdef PRSS
	       , int *s_score, int *s_n1

#endif
	       )
{
  double zscore;
  int i_score;
  struct buf_str *p_rbuf, *cur_buf_p;
  int i, t_best, t_rbest, t_qrbest, tm_best, t_n1, sc_ix;
  double e_score, tm_escore, t_rescore, t_qrescore;
  int jstats;

  sc_ix = pst.score_ix;

  cur_buf_p = cur_buf->buf;
  
  t_best = t_rbest = t_qrbest = -1;
  tm_escore = t_rescore = t_qrescore = FLT_MAX;

  while (cur_buf->buf_cnt--) { /* count down the number of results */
    p_rbuf = cur_buf_p++;	/* step through the results buffer */

    i_score = p_rbuf->rst.score[sc_ix];
    e_score = p_rbuf->rst.escore;

    /* need to look for frame 0 if TFASTA, then save stats at frame 6 */
    if (fdata) {
      fprintf(fdata,
	      "%-12s %5d %6d %d %.5f %.5f %4d %4d %4d %g %d %d %8ld\n",
	      p_rbuf->libstr,
#ifdef SUPERFAMNUM
	      sfn_cmp(qsfnum,p_rbuf->sfnum),
#else
	      0,
#endif
	      p_rbuf->n1,p_rbuf->frame,p_rbuf->rst.comp,p_rbuf->rst.H,
	      p_rbuf->rst.score[0],p_rbuf->rst.score[1],p_rbuf->rst.score[2],
	      p_rbuf->rst.escore, p_rbuf->rst.segnum, p_rbuf->rst.seglen, p_rbuf->lseek);
    }

#ifdef PRSS
    if (p_rbuf->lseek==0) {
      s_score[0] = p_rbuf->rst.score[0];
      s_score[1] = p_rbuf->rst.score[1];
      s_score[2] = p_rbuf->rst.score[2];
      *s_n1 = p_rbuf->n1;

      bestp = bestp_arr[nbest++];
      bestp->score[0] = s_score[0];
      bestp->score[1] = s_score[1];
      bestp->score[2] = s_score[2];
      bestp->n1 = *s_n1;
      bestp->escore = p_rbuf->rst.escore;
      bestp->segnum = p_rbuf->rst.segnum;
      bestp->seglen = p_rbuf->rst.seglen;
      bestp->zscore = zscore;
      bestp->lseek = p_rbuf->lseek;
      bestp->m_file_p = p_rbuf->m_file_p;
      memcpy(bestp->libstr,p_rbuf->libstr,MAX_UID);
      bestp->n1tot_p = p_rbuf->n1tot_p;
      bestp->frame = p_rbuf->frame;

      continue;
    }
#endif

    t_n1 = p_rbuf->n1;
    if (i_score > t_best) tm_best = t_best = i_score;
    if (e_score < tm_escore) tm_escore = e_score;

    if (m_msg.qshuffle) {
      if (p_rbuf->qr_score > t_qrbest)
	t_qrbest = p_rbuf->qr_score;
      if (p_rbuf->qr_escore < t_qrescore)
	t_qrescore = p_rbuf->qr_escore;
      
      if (p_rbuf->frame == m_msg.nitt1 && nqstats < m_msg.shuff_max) {
	qstats[nqstats].n1 = p_rbuf->n1;	/* save the best score */
	qstats[nqstats].comp =  p_rbuf->rst.comp;
	qstats[nqstats].H = p_rbuf->rst.H;
	qstats[nqstats].escore = t_qrescore;
	qstats[nqstats++].score = t_qrbest;
	t_qrbest = -1;	/* reset t_qrbest, t_qrescore */
	t_qrescore = FLT_MAX;
      }
    }

    if (pst.zsflag >= 10 && p_rbuf->r_score > t_rbest) {
      t_rbest = p_rbuf->r_score;
      t_rescore = p_rbuf->r_escore;
    }

    /* statistics done for best score of set */


    if (p_rbuf->frame == m_msg.nitt1) {
      if (nstats < MAXSTATS ) {
	stats[nstats].n1 = t_n1;
	stats[nstats].comp = p_rbuf->rst.comp;
	stats[nstats].H = p_rbuf->rst.H;
	if (pst.zsflag >= 10) {
	  tm_best = t_rbest;
	  tm_escore = t_rescore;
	  t_rbest = -1;
	  t_rescore = FLT_MAX;
	}
	stats[nstats].escore  = tm_escore;
	stats[nstats++].score = tm_best;
	t_best = -1;
	tm_escore = FLT_MAX;
      }
      else if (pst.zsflag > 0) {
	if (!stats_done) {
	  pst.zsflag_f = process_hist(stats,nstats,m_msg,pst,
				      histp, pstat_voidp,0);
	  kstats = nstats;
	  stats_done = 1;
	  for (i=0; i<MAXBEST; i++) {
	    bestp_arr[i]->zscore = 
	      (*find_zp)(bestp_arr[i]->score[pst.score_ix],
			 bestp_arr[i]->escore, bestp_arr[i]->n1,
			 bestp_arr[i]->comp, *pstat_voidp);
	  }
	}
#ifdef SAMP_STATS
	else {
	  if (!m_msg.escore_flg) {
	    jstats = nrand(++kstats);
	    if (jstats < MAXSTATS) {
	      stats[jstats].n1 = t_n1;
	      stats[jstats].comp = p_rbuf->rst.comp;
	      stats[jstats].H = p_rbuf->rst.H;
	      if (pst.zsflag >= 10) {
		tm_best = t_rbest;
	      }
	      stats[jstats].score = tm_best;
	    }
	  }
	}
#endif
      }
    }

    /* best saved for every score */
    if (stats_done) {

      zscore=(*find_zp)(i_score, e_score, p_rbuf->n1,(double)p_rbuf->rst.comp,
			*pstat_voidp);

      if (p_rbuf->frame == m_msg.nitt1) {
	addhistz((*find_zp)(t_best, tm_escore, p_rbuf->n1, (double) p_rbuf->rst.comp,
			    *pstat_voidp), histp);
	t_best = t_rbest = -1;
	tm_escore = t_rescore = FLT_MAX;
      }
    }
    else zscore = (double) i_score;

#ifndef PRSS
    if (zscore > zbestcut) {
      if (nbest >= MAXBEST) {
	bestfull = nbest-MAXBEST/4;
	selectbestz(bestp_arr,bestfull-1,nbest);
	zbestcut = bestp_arr[bestfull-1]->zscore;
	nbest = bestfull;
      }
      bestp = bestp_arr[nbest++];
      bestp->score[0] = p_rbuf->rst.score[0];
      bestp->score[1] = p_rbuf->rst.score[1];
      bestp->score[2] = p_rbuf->rst.score[2];
      bestp->comp = (double) p_rbuf->rst.comp;
      bestp->H = (double) p_rbuf->rst.H;
      bestp->escore = p_rbuf->rst.escore;
      bestp->segnum = p_rbuf->rst.segnum;
      bestp->seglen = p_rbuf->rst.seglen;
      bestp->zscore = zscore;
      bestp->lseek = p_rbuf->lseek;
      memcpy(bestp->libstr,p_rbuf->libstr,MAX_UID);
      bestp->cont = p_rbuf->cont; /* not cont+1 because incremented already */
      bestp->m_file_p = p_rbuf->m_file_p;
      bestp->n1 = p_rbuf->n1;
      bestp->n1tot_p = p_rbuf->n1tot_p;
      bestp->frame = p_rbuf->frame;
      bestp->nsfnum = p_rbuf->nsfnum;
#ifdef SUPERFAMNUM
      if ((bestp->sfnum[0] = p_rbuf->sfnum[0])>0 &&
	  (bestp->sfnum[1] = p_rbuf->sfnum[1])>0 &&
	  (bestp->sfnum[2] = p_rbuf->sfnum[2])>0 &&
	  (bestp->sfnum[3] = p_rbuf->sfnum[3])>0 &&
	  (bestp->sfnum[4] = p_rbuf->sfnum[4])>0 &&
	  (bestp->sfnum[5] = p_rbuf->sfnum[5])>0 &&
	  (bestp->sfnum[6] = p_rbuf->sfnum[6])>0 &&
	  (bestp->sfnum[7] = p_rbuf->sfnum[7])>0 &&
	  (bestp->sfnum[8] = p_rbuf->sfnum[8])>0 &&
	  (bestp->sfnum[9] = p_rbuf->sfnum[9])>0) ;
#endif
    }
#endif
  }
}
#endif
