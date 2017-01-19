/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: p2_complib.c,v 1.96 2007/01/12 20:15:16 wrp Exp $ */

/*
 * pcomplib.c : Parallel library search
 * 
 *	#define FIRSTNODE 0/1 (in msg.h) can be used to reserve one node
 *	for collecting results
 *
 * Parallel specific options (from doinit.c):
 *	-J # jump to query #
 *	-I   self-comparison, do (N choose 2) comparisons
 *	-T # number of workers
 */

/* This version is modifed to read all files, query and database,
   through the manager process. Workers will now receive their
   database from the manager, rather than reading it themselves.  This
   cuts down considerably on NFS traffic, simplifies searches of
   multiple files, and allows use of clusters of slave nodes that do
   not have NFS access
*/

/* modified 5-November-2004 to ensure 15 byte (SEQ_PAD) NULL
   padding

   modified 12-December-2006 to ensure n0>0 before SEQ_PAD padding.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include <limits.h>
#include <float.h>
#include <math.h>

#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <sys/stat.h>

#ifdef PVM_SRC
#include "pvm3.h"
char *mp_verstr="34.26, January 12, 2007 PVM";
#endif

#ifdef MPI_SRC
#include "mpi.h"
char *mp_verstr="34.26, January 12, 2007 MPI";
#endif

#include "msg.h"
#include "defs.h"
#include "mm_file.h"

#include "structs.h"
#include "param.h"
#include "p_mw.h"

#define XTERNAL
#include "uascii.h"

char pgmdir[MAX_FN];
char workerpgm[MAX_FN];
char managepgm[MAX_FN];

#define XTERNAL
#include "upam.h"
#undef XTERNAL

/********************************/
/* global variable declarations */
/********************************/
char gstring2[MAX_STR];                  /* string for label */
char gstring3[MAX_STR];                  /* string for label */
char hstring1[MAX_STR];

int nsfnum;	/* number of superfamily numbers */
int sfnum[10];	/* superfamily number from types 0 and 5 */
int nsfnum_n;
int sfnum_n[10];

/********************************/
/* extern variable declarations */
/********************************/
extern char *prog_func;		/* function label */
extern char *verstr, *iprompt0, *iprompt1, *iprompt2, *refstr;

/********************************/
/*extern function declarations  */
/********************************/

void libchoice(char *lname, int, struct mngmsg *); /* lib_sel.c */
void libselect(char *lname, struct mngmsg *);	/* lib_sel.c */

extern void closelib();
/* check for DNA sequence (nxgetaa.c) */
extern int scanseq(unsigned char *seq, int n, char *str);
extern void re_ascii(int *qascii, int *sascii);
extern int recode(unsigned char *seq, int n, int *qascii, int nsq);

/* 1d to 2d pam (initxx.c) */
extern void initpam2 (struct pstruct *ppst);
/* initialize environment (doinit.c) */
extern void h_init (struct pstruct *ppst, struct mngmsg *, char *);
extern void s_abort (char *p,  char *p1);
extern void query_parm (struct mngmsg *m_msp, struct pstruct *ppst);
extern void last_init (struct mngmsg *, struct pstruct *, int);

extern void initenv (int argc, char **argv, struct mngmsg *m_msg,
		     struct pstruct *ppst, unsigned char **aa0);

/* print hist, summaries, timing information */
void prhist(FILE *, struct mngmsg, struct pstruct, struct hist_str, int nstats, struct db_str, char *);
void printsum(FILE *);
extern void ptime (FILE *, time_t);

/* reset parameters if DNA sequence (initxx.c) */
extern void resetp (struct mngmsg *, struct pstruct *);

/* read a sequence (nmgetlib.c) */
struct lmf_str *openlib(char *, int, int *, int, struct lmf_str *);

#define QGETLIB (q_file_p->getlib)
#define LGETLIB (l_file_p->getlib)

/* these functions are in scaleswn.c */
extern int process_hist(struct stat_str *sptr, int nstat,
			struct mngmsg m_msg, struct pstruct pst,
			struct hist_str *hist, void **pstat_void, int);
extern double zs_to_E(double zs, int n1, int isdna, long, struct db_str ntt);
extern double (*find_zp)(int score, double escore, int length, double comp, void *);
void addhistz(double zscore, struct hist_str *);	/* scaleswn.c */
void last_stats(const unsigned char *aa0, int n0,  
		struct stat_str *sptr, int nstats,
		struct beststr **bestp_arr, int nbest,
		struct mngmsg m_msg, struct pstruct pst,
		struct hist_str *histp, void *rs);

void selectbestz(struct beststr **, int, int);
void sortbest(struct beststr **, int, int);

void showbest (FILE *fp, struct beststr **bptr, int nbest,
	       int qlib, struct mngmsg *m_msg, struct pstruct pst,
	       struct db_str ntt, char *gstring2);

void showalign (FILE *fp, 
		struct beststr **bptr, int nbest,int qlib, struct mngmsg m_msg,
		struct pstruct pst, char *gstring2);

#ifdef PVM_SRC
char worknode[120];
int pinums[MAXNOD],hosttid;
int narch;
struct pvmhostinfo *hostp;
#endif

FILE *outfd;			/* Output file */

extern time_t s_time ();                 /* fetches time for timing */

/* this information is global for fsigint() */
time_t tstart, tscan, tprev, tdone;	/* Timing */
time_t tdstart, tddone, time();
int max_nodes, nnodes;			/* number of nodes */
int node_map[MAXWRKR], node_id[MAXWRKR];
int tot_speed,h_speed;
int  qlib = 0;	/* number of sequences scanned */
struct db_str ntt, qtt;

extern int max_workers, worker_1, worker_n;
int  wlsn [MAXWRKR + 1];	/* number of library sequences in worker */
int  clsn [MAXWRKR + 1];	/* number of 1st library sequence in worker */

int max_buf_cnt;

#ifdef PVM_SRC
#ifndef WORKERPGM
#define WORKERPGM "c34.work"
#endif
#endif

main (int argc, char *argv[])
{
  unsigned char *aa00, *aa01, *aa0p0, *aa0p1;
  unsigned char *aa1, *aa1ptr, *aa1prev;
  int aa1i, *aa1i_arr;	/* integer offset of sequence in buffer */

  int n1;
  int *n1tot_ptr=NULL, *n1tot_cur;
  int n1tot_cnt=0;
  int n1tot_v;

  long l_off;
  char nodefile[240];
  struct pstruct pst;
  int i_score;
  struct lmf_str *q_file_p;
  struct lmf_str *l_file_p;

  /* from manage code */
  struct mngmsg m_msg0, m_msg1;	/* Message from host to manager */
  struct mngmsg *m_msp0, *m_msp1;	/* alternating pointers */
  struct qmng_str qm_msg0, qm_msg1;	/* stuff updated for each query */
  char q_sqnam[4]; 
  int sstart, sstop;
    
  struct qmng_str *qm_msp0, *qm_msp1;	/* pointer to stuff updated */
  int last_msg_b[10];	/* last set of numbers */
  long curtype = ONETYPE;	/* current message type */
  int nclib;
  struct beststr *best,		/* array of best scores */
                 **bptr;	/* array of pointers */
  struct comstr bestr[BFR+1];	/* temporary structure array */
  struct comstr2 bestr2[BFR2+1];	/* temporary structure array */
  struct a_struct *aln_d_base=NULL;	/* alignment info for -m 9 */
  int qres_bufsize;		/* buffer size for results */
  struct stat_str *stats=NULL, *qstats=NULL;
  int best_flag = 1;		/* bptr[] must be re-initialized */
  int fast_flag = 0;		/* send new sequences before old displayed */
  int nstats, nqstats, kstats, jstats;
  int nbest, nres;		/* number of best scores */
  double zbestcut = -BIGNUM;	/* z-value cutoff */
  int lcnt;			/* counters */
  int nopt;
  int i, j, k, is, id, iw, ires, naa0 = 0;

  FILE *fdata=NULL;		/* file for full results */
  struct sql *desptr;
  struct sql *ldes;		/* descriptive lines for all lib sequences */
  char *bline_buf, *bline_bufp;
  char *bline_buf_mx;	/* buffer for blines */
  char q_bline[256];
  char t_bline[256];
  int max_bline_b, bline_inc;
  int *n1_arr, *m_seqnm_arr;
  unsigned char *aa1_buf;

  char tlibstr[11];		/* used only for fdata *.res files */
  
  int node, snode, zero;	/* Number of nodes */
  int bufid, numt, tid;

  int ave_seq_len;
  int max_sql;
  int ntbuff, nseq, m_seqnm;
  int iln, ocont, maxt;
  long loffset;

  int leng;			/* leng is length of the descriptive line */
  fseek_t qseek,lseek;		/* seek into library of current sequence */
  int qlcont,lcont;			/* continued sequence */
  int n_proc, n_tmp;
  char errstr[120];
  int stats_done =0;			/* flag for z-value processing */
  int tm_best, t_rbest, t_qrbest, t_best, t_n1;
  double e_score, tm_escore, t_rescore, t_qrescore;
  double zscore;			/* tmp value */
  double k_H, k_comp;
  char tmp_str[MAX_FN];
  char pgm_abbr[MAX_SSTR];
  char *bp;
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  void fsigint();
  
  signal(SIGHUP,SIG_IGN);
  if (signal(SIGINT,SIG_IGN) != SIG_IGN) signal(SIGINT,fsigint);
  if (signal(SIGQUIT,SIG_IGN) != SIG_IGN) signal(SIGQUIT,fsigint);
/*  if (signal(SIGSEGV,SIG_IGN) != SIG_IGN) signal(SIGSEGV,fsigint); */

  /* Initialization */


#if defined(UNIX)
  m_msg0.quiet = !isatty(1);
#endif

  /* BFR must be %6 = 0 for TFASTA */
  if ((BFR%6) != 0) {
    fprintf(stderr," BFR size %d not %%6=0 - recompile\n",BFR);
    exit(1);
  }

#ifdef MPI_SRC
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&tid);
  if (tid > 0) {
    workcomp(tid); 
    MPI_Finalize();
    exit(0);
  }
#endif

  printf("#");
  for (i=0; i<argc; i++) {
    if (strchr(argv[i],' ')) printf(" \"%s\"",argv[i]);
    else printf(" %s",argv[i]);
  }
  printf("\n");

#ifdef MPI_SRC
  MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
  if (nnodes <= 1) {
    fprintf(stderr," nnodes = %d; no workers available\n",nnodes);
    exit(1);
  }
  else fprintf(stderr," have %d nodes\n",nnodes);

  tot_speed = nnodes*100;
#endif

  h_init (&pst,&m_msg0, pgm_abbr);

  initenv (argc, argv, &m_msg0, &pst, &aa00);

#ifdef PVM_SRC
  strncpy (workerpgm, WORKERPGM,sizeof(workerpgm)-1);
  strncat(workerpgm, pgm_abbr, sizeof(workerpgm)-strlen(workerpgm)-1);
  workerpgm[sizeof(workerpgm)-1] = '\0';
#endif
  
  strncpy(q_sqnam,"aa",sizeof(q_sqnam));
  m_msg0.quiet = 1;
  if (m_msg0.qdnaseq != SEQT_UNK && 
      (m_msg0.qdnaseq == SEQT_DNA || m_msg0.qdnaseq == SEQT_RNA))
    strncpy(q_sqnam,"nt",sizeof(q_sqnam));

  m_msg0.pstat_void = NULL;
  m_msg0.hist.hist_a = NULL;

  fprintf (stderr, "Pcomp library processor\n");
  fprintf (stderr, "Using %s\n", prog_func);
  
  tstart = tscan = s_time();
  tdstart = time(NULL);
  

#ifdef PVM_SRC
  if ((hosttid=pvm_mytid())<0) {
    pvm_perror("initialization");
    fprintf(stderr,"can't initialize %s\n", argv[0]);
    pvm_exit();
    exit(1);
  }
  
  pvm_config(&nnodes,&narch,&hostp);
  fprintf(stderr,"nnodes: %d, narch: %d\n",nnodes, narch);
  max_nodes = nnodes;

#ifdef DEBUG
  pvm_catchout(stderr);
#endif

/*  if (nnodes < 2 ) nnodes = 4; */
  if (max_workers > 0  && nnodes > max_workers) {
    nnodes = max_workers+FIRSTNODE;
    fprintf(stderr," workers reset from %d to %d\n",
	    max_nodes,nnodes-FIRSTNODE);
  }
  else max_workers = nnodes;
  
  strncpy(nodefile,pgmdir,sizeof(nodefile)-1);
  strncat(nodefile,workerpgm,sizeof(nodefile)-strlen(nodefile)-1);
  nodefile[sizeof(nodefile)-1] = '\0';

  if (worker_1 > 0) {
    /* remap configuration to specific nodes */
    for (i=FIRSTNODE, j=worker_1; i<nnodes && j<=worker_n; i++,j++)
      node_id[i]=j;
    nnodes = i;
    max_workers = i-FIRSTNODE;
    fprintf(stderr," workers remapped from %d to %d\n",
	    max_nodes,nnodes-FIRSTNODE);
    max_nodes = nnodes;
  }
  else {
    for (i=0; i< nnodes; i++) node_map[i]=node_id[i] = i;
  }

  if (nnodes < max_nodes) {
    hostp++;	/* bump over host name for spawn */
    rand_nodes(node_map,nnodes,max_nodes-1);
    for (i=FIRSTNODE; i<nnodes; i++) {
      numt+=pvm_spawn(nodefile,NULL,PvmTaskHost,hostp[node_map[i]].hi_name,
		      1,&pinums[i]);
    }
  }
  else {
    /* i counts through nodes (machines) */
    /* j counts through processes (multiple processes/node) */
    /* node map maps the process (virtual node) to a physical node (machine) */

    for (i=j=FIRSTNODE; i<nnodes && j < MAXWRKR; i++) {
      n_proc = hostp[node_id[i]].hi_speed%100;
      if (n_proc == 0) n_proc = 1;
      if (n_proc > max_workers) n_proc = max_workers;

      n_tmp =pvm_spawn(nodefile,NULL,PvmTaskHost,hostp[node_id[i]].hi_name,
		       n_proc,&pinums[j]);
      if (n_tmp < n_proc)
	fprintf(stderr," spawn problem: %d\n", pinums[j]);
      if (n_tmp > 0) {
	for (k=j; k < j+n_tmp; k++) node_map[k]=node_id[i];
	j += n_tmp;
      }
    }
    nnodes = numt = j;
  }

  if (numt < nnodes) {
    if (numt <= 0) {
      pvm_perror("");
      pvm_exit();
      exit(1);
    }
    nnodes = numt;
  }

  for (tot_speed=0,i=FIRSTNODE; i<nnodes; i++) {
    if (pinums[i]<0) {
      fprintf(stderr," tids %d %8o\n",i,pinums[i]);
      pvm_perror("");
      pvm_exit();
      exit(1);
    }
    else {
      h_speed = hostp[node_map[tidtonode(pinums[i])]].hi_speed;
      if (h_speed <= 0) h_speed = 100;
      fprintf(stderr," tids %d %8o %s %5d\n",i,pinums[i],
	      hostp[node_map[tidtonode(pinums[i])]].hi_name,
	      h_speed);
      tot_speed +=(hostp[node_map[tidtonode(pinums[i])]].hi_speed);
    }
  }

  strncpy(worknode,nodefile,sizeof(worknode));
  fprintf (stderr, "%3d worker programs loaded from %s\n",
	   nnodes-FIRSTNODE,worknode);
#endif  

  /* need to allocate two aa0 arrays so that the old is saved for alignments */

  /* Allocate space for the query sequence */
  if ((aa00 = (unsigned char *) malloc ((MAXTST + SEQ_PAD + 1)* sizeof (char))) == NULL)
    s_abort ("Unable to allocate query sequence", "");
  
  if ((aa01 = (unsigned char *) malloc ((MAXTST + SEQ_PAD + 1) * sizeof (char))) == NULL)
    s_abort ("Unable to allocate query sequence", "");
  
  fputs(iprompt0,stdout);
  fprintf(stdout," %s%s\n",verstr,refstr);

  /* Query library */
  if (m_msg0.tname[0] == '\0') {
      if (m_msg0.quiet == 1) s_abort("query sequence undefined","");
      
      fprintf(stderr, "Pvcomplib [%s]\n",mp_verstr);
    l1:	fputs (iprompt1, stdout);
      fflush  (stdout);
      if (fgets (m_msg0.tname, 80, stdin) == NULL)
	s_abort ("Unable to read query library name","");
      if ((bp=strchr(m_msg0.tname,'\n'))!=NULL) *bp='\0';
      if (m_msg0.tname[0] == '\0') goto l1;
    }
  
  /* Open query library */
  if ((q_file_p=
       openlib(m_msg0.tname, m_msg0.qdnaseq,qascii,!m_msg0.quiet,NULL))==NULL) {
      s_abort(" cannot open library ",m_msg0.tname);
    }
  /*
  else {
    printf ("searching %s library\n",m_msg0.tname);
  }
  */

  ntt.entries = qtt.entries = 0;
  ntt.carry = qtt.carry = 0;
  ntt.length = qtt.length = 0l;

  /* Fetch first sequence */
  qlcont = 0;
  while (qlib < m_msg0.ql_start) {	/* skip through query sequences */
    pst.n0 = qm_msg0.n0 = m_msg0.n0 = 
      QGETLIB (aa00, MAXTST, q_bline, sizeof(q_bline), &qseek, &qlcont,
	       q_file_p,&m_msg0.sq0off);

    strncpy(qm_msg0.libstr,q_bline,sizeof(qm_msg0.libstr)-20);
    qm_msg0.libstr[sizeof(qm_msg0.libstr)-21]='\0';
    if ((bp=strchr(qm_msg0.libstr,' '))!=NULL) *bp='\0';

    /* if annotations are included in sequence, remove them */
    if (m_msg0.ann_flg) {
      pst.n0 = qm_msg0.n0 = m_msg0.n0 = 
	ann_scan(aa00, m_msg0.n0, &m_msg0, m_msg0.qdnaseq);
#ifdef DEBUG
      fprintf(stderr,"m_msp0->/aa0a is: %o/%o\n",&m_msg0,m_msg0.aa0a);
#endif
    }

    if (m_msg0.term_code && 
	!(m_msg0.qdnaseq == SEQT_DNA || m_msg0.qdnaseq==SEQT_RNA) &&
	aa00[m_msg0.n0-1]!='*') {
      aa00[m_msg0.n0++]='*';
      aa00[m_msg0.n0]=0;
      pst.n0 = qm_msg0.n0 = m_msg0.n0;
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
      for (id=0,is=sstart; is<min(m_msg0.n0,sstop); ) aa00[id++]=aa00[is++];
      aa00[id]=0;
      pst.n0 = qm_msg0.n0 = m_msg0.n0 = min(m_msg0.n0,sstop)-sstart;
      if (m_msg0.sq0off==1) m_msg0.sq0off = sstart+1;
    }

    qlib++;

    if (m_msg0.n0 <= 0)
      s_abort ("Unable to fetch sequence from library: ", m_msg0.tname);
  }
  qtt.entries=1;
  qm_msg0.slist = 0;

  /* now have correct query sequence - check sequence type and reset */
  if (m_msg0.qdnaseq == SEQT_UNK) {	/* check for DNA sequence */
    if (m_msg0.n0 > 20 &&
	(float)scanseq(aa00,m_msg0.n0,"ACGTUNacgtun")/(float)m_msg0.n0>0.85) {
      pascii = nascii;
      m_msg0.qdnaseq = SEQT_DNA;
    }
    else {	/* its protein */
      pascii = aascii;
      m_msg0.qdnaseq = SEQT_PROT;
    }

    re_ascii(qascii,pascii);
    init_ascii(pst.ext_sq_set,qascii,m_msg0.qdnaseq);
    m_msg0.n0 = recode(aa00,m_msg0.n0,qascii,pst.nsqx);
  }

  /* for ALTIVEC, must pad with 15 NULL's */
  for (i=0; i<SEQ_PAD+1; i++) {aa00[m_msg0.n0+i]=0;}

  qtt.length = m_msg0.n0;

  if (qlib <= 0) {
    fprintf(stderr," no sequences found in query library\n");
    exit(1);
  }

  resetp (&m_msg0, &pst);

  sprintf(tmp_str," %d %s", qm_msg0.n0, q_sqnam);
  leng = strlen (qm_msg0.libstr);
  if (leng + strlen(tmp_str) >= sizeof(qm_msg0.libstr))
    qm_msg0.libstr[sizeof(qm_msg0.libstr)-strlen(tmp_str)-2] = '\0';
  strncat(&qm_msg0.libstr[0],tmp_str,
	  sizeof(qm_msg0.libstr)-strlen(qm_msg0.libstr)-1);
  qm_msg0.libstr[sizeof(qm_msg0.libstr)-1]='\0';

  qm_msg0.seqnm = qlib-1;
  
  /* Library */

  if (strlen (m_msg0.lname) == 0) {
    if (m_msg0.quiet == 1) s_abort("library name undefined","");
    libchoice(m_msg0.lname, sizeof(m_msg0.lname), &m_msg0);
  }

  libselect(m_msg0.lname, &m_msg0);

  /* Get additional parameters here */
  if (!m_msg0.quiet) query_parm (&m_msg0, &pst);
  
  last_init(&m_msg0, &pst,nnodes-FIRSTNODE);
  memcpy(&m_msg1, &m_msg0, sizeof(m_msg0));
  
  /* m_msg0.maxn needs to be set to MAXLIB or MAXTRN, depending on the
     function - max_tot has the MAXTST + (MAXLIB|MAXTRN) */
  if (m_msg0.maxn <= 0) m_msg0.maxn = m_msg0.max_tot - MAXTST;

  if (m_msg0.maxn < 2 * m_msg0.dupn) m_msg0.maxn = 5*m_msg0.dupn;
  pst.maxlen = m_msg0.maxn;

  m_msg0.loff = m_msg0.dupn;
  m_msg0.maxt3 = m_msg0.maxn-m_msg0.loff;


  /* ******************** */
  /* initial manager code */
  /* ******************** */
  
  outfd = stdout;
  if (m_msg0.outfile[0]!='\0') {
    if ((outfd = fopen(m_msg0.outfile,"w"))==NULL) {
      fprintf(stderr, "cannot open %s for output\n", m_msg0.outfile);
      outfd = stdout;
    }
  }

  /* Label the output */
  printf("Query library %s vs %s library\n", m_msg0.tname, m_msg0.lname);
  
  /* Allocate space for saved scores */
  if ((best = 
       (struct beststr *)malloc((MAXBEST+1)*sizeof(struct beststr)))==NULL)
    s_abort ("Cannot allocate best struct","");
  if ((bptr = 
       (struct beststr **)malloc((MAXBEST+1)*sizeof(struct beststr *)))==NULL)
    s_abort ("Cannot allocate bptr","");
  
  /* Initialize bptr */
  for (nbest = 0; nbest < MAXBEST+1; nbest++)
    bptr[nbest] = &best[nbest];

  best++; bptr++;
  best[-1].score[0]=best[-1].score[1]=best[-1].score[2]=INT_MAX;
  best[-1].zscore = FLT_MAX;
  best[-1].escore = FLT_MIN;
  best_flag = 0;
  
  if ((stats =
       (struct stat_str *)calloc((size_t)MAXSTATS,sizeof(struct stat_str)))
      ==NULL)
    s_abort ("Cannot allocate stats struct","");
  nstats = 0;

  /* Now open the second library, divide it, send sequences to all workers */
  /* Set up buffer for reading the library:

     We will start by using a 2 Mbyte buffer for each worker.  For
     proteins, that means 5,000 sequences of length 400 (average).
     For DNA, that means 2,000 sequences of length 1000.  At the moment,
     those are good averages.
     */

  if (max_buf_cnt <= 0) {
    if (m_msg0.ldnaseq==SEQT_DNA) max_buf_cnt = MAX_NT_BUF;
    else max_buf_cnt = MAX_AA_BUF;
  }

  if (m_msg0.ldnaseq==SEQT_DNA) ave_seq_len = AVE_NT_LEN;
  else ave_seq_len = AVE_AA_LEN;

  /* however - buffer sizes should be a function of the number of
     workers so that all the workers are kept busy.  Assuming a 10,000
     entry library is the smallest we want to schedule, then
     */

  if (max_buf_cnt > 10000/(nnodes-FIRSTNODE)) 
    max_buf_cnt = 10000/(2*(nnodes-FIRSTNODE));

  /* allocate space for sequence buffers */

  m_msg0.pbuf_siz=max_buf_cnt*ave_seq_len;
  if (m_msg0.pbuf_siz < 5*m_msg0.maxn)
    m_msg0.pbuf_siz = 5*m_msg0.maxn;

#ifdef PVM_SRC 
#ifdef ROUTE_DIRECT
  pvm_setopt(PvmRoute,PvmRouteDirect);
#endif
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&nnodes,1,1);
  pvm_pkint(pinums,nnodes,1);
  pvm_pkbyte((char *)&m_msg0,(int)sizeof(m_msg0),1);
  for (node = FIRSTNODE; node<nnodes; node++) 
    if (pvm_send(pinums[node],STARTTYPE0)<0) {
      pvm_perror("pvm_send1");
      pvm_exit();
      exit(1);
    }
#endif
#ifdef MPI_SRC
  for (node = FIRSTNODE; node<nnodes; node++)  {
    MPI_Send(&m_msg0,(int)sizeof(m_msg0),MPI_BYTE,node,STARTTYPE0,
	     MPI_COMM_WORLD);
  }
#endif

  /* now send pst, sascii */
#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkbyte((char *)&pst,(int)sizeof(pst),1);
  pvm_pkbyte((char *)pascii,(int)sizeof(aascii),1);

  for (node = FIRSTNODE; node< nnodes; node++)
    pvm_send(pinums[node],STARTTYPE1);

  /* send pam12 */
  pvm_initsend(PvmDataRaw);
  pvm_pkint(pam12,m_msg0.pamd1*m_msg0.pamd2,1);
  for (node = FIRSTNODE; node< nnodes; node++)
    pvm_send(pinums[node],STARTTYPE2);

  /* send pam12x */
  pvm_initsend(PvmDataRaw);
  pvm_pkint(pam12x,m_msg0.pamd1*m_msg0.pamd2,1);
  for (node = FIRSTNODE; node< nnodes; node++)
    pvm_send(pinums[node],STARTTYPE3);

#endif
#ifdef MPI_SRC
  for (node=FIRSTNODE; node < nnodes; node++) {
    MPI_Send(&pst,(int)sizeof(pst),MPI_BYTE,node,STARTTYPE1,
	     MPI_COMM_WORLD);
    MPI_Send(pascii,(int)sizeof(aascii),MPI_BYTE,node,STARTTYPE1,
	     MPI_COMM_WORLD);
    MPI_Send(pam12,m_msg0.pamd1*m_msg0.pamd2,MPI_INT,node,STARTTYPE2,
	     MPI_COMM_WORLD);
    MPI_Send(pam12x,m_msg0.pamd1*m_msg0.pamd2,MPI_INT,node,STARTTYPE3,
	     MPI_COMM_WORLD);
  }
#endif

  if ((n1_arr =
       (int *)calloc((size_t)(max_buf_cnt+1),sizeof(int)))
      ==NULL) {
    fprintf(stderr," cannot allocate n1_arr %d\n",max_buf_cnt+1);
    s_abort(" cannot allocate n1_arr","");
    exit(1);
  }

  if ((aa1i_arr =
       (int *)calloc((size_t)(max_buf_cnt+1),sizeof(int)))
      ==NULL) {
    fprintf(stderr," cannot allocate aa1i_arr %d\n",max_buf_cnt+1);
    s_abort(" cannot allocate aa1i_arr","");
    exit(1);
  }

  if ((m_seqnm_arr=
       (int *)calloc((size_t)(max_buf_cnt+1),sizeof(int)))
      ==NULL) {
    fprintf(stderr," cannot allocate m_seqnm_arr %d\n",max_buf_cnt+1);
    s_abort(" cannot allocate m_seqnm_arr","");
    exit(1);
  }

  if ((aa1_buf =
       (unsigned char *)calloc((size_t)(m_msg0.pbuf_siz),sizeof(unsigned char)))
      ==NULL) {
    s_abort(" cannot allocate library buffer %d","");
    exit(1);
  }


  /* also allocate space for descriptions.  Assume max of 250,000 sequences/
     worker for now
  */

  /* max_sql is the maxinum number of library sequences that can be stored */
  max_sql = MAXSQL;

  if ((ldes=(struct sql *)calloc(max_sql,sizeof(struct sql)))==NULL) {
    fprintf(stderr," failure to allocate ldes(%d) %ld\n",
	    max_sql,max_sql*sizeof(struct sql));
    s_abort("cannot allocate ldes","");
    exit(1);
  }

  max_bline_b = MAXSQL * (m_msg0.aln.llen+1)/4;
  bline_inc = m_msg0.aln.llen;
  if (m_msg0.markx & MX_M9SUMM) bline_inc += 40;

  i = 4;
  while (i-- > 0) {
    if ((bline_buf=(char *)calloc(max_bline_b,sizeof(char)))!=NULL) break;
    max_bline_b /= 2;
    bline_inc /= 2;
  }
  if (bline_buf == NULL) {
    fprintf(stderr," failure to allocate bline_buf(%d) %d\n",
	    max_sql,max_bline_b);
    s_abort(" cannot allocate bline_buf","");
  }

  bline_bufp = bline_buf;
  bline_buf_mx = bline_buf+max_bline_b;

  /* the code for filling the buffers is copied from comp_thr.c */
  /* the major differences reflect the fact that all library descriptions
     will be kept in memory, indexed by sequence number.

     As a result, one buffer is filled by this loop -
       ldes[] has the descriptive information for every sequence
       this array could potentially be quite large
  */

  /* now open the library and start reading */
  /* get a buffer and fill it up */

  ntbuff = 0;
  m_seqnm = 0;	/* m_seqnm is the number of this library sequence */
  nseq = 0;

  node = FIRSTNODE;

  /* sqs2_buf[0].aa1 = aa1_buf; */
  aa1 = aa1_buf;

  /* iln counts through each library */
  for (iln = 0; iln < m_msg0.nln; iln++) {
    if ((l_file_p=
	 openlib(m_msg0.lbnames[iln], m_msg0.ldnaseq,lascii,!m_msg0.quiet,NULL))==NULL) {
      fprintf(stderr," cannot open library %s\n",m_msg0.lbnames[iln]);
      continue;
    }
    else {
      printf ("searching %s library\n",m_msg0.lbnames[iln]);
    }

    lcont = ocont = 0;
    n1tot_v = n1tot_cnt = 0;
    n1tot_ptr = n1tot_cur = NULL;
    maxt = m_msg0.maxn;
    loffset = 0l;

    /* read sequence directly into buffer */
    aa1ptr = aa1; /* = sqs2_buf[0].aa1; */

    while ((n1= LGETLIB(aa1ptr,maxt,t_bline,sizeof(t_bline),&lseek,&lcont,
			l_file_p,&l_off))>=0) {

      /* skip sequences outside range */
      if (n1 < m_msg0.n1_low || n1 > m_msg0.n1_high) goto loop1;
      
      /* add termination code for proteins, if asked */
      if (m_msg0.term_code && !lcont && 
	  m_msg0.ldnaseq==SEQT_PROT && aa1ptr[n1-1]!=m_msg0.term_code) {
	aa1ptr[n1++]=m_msg0.term_code;
	aa1ptr[n1]=0;
      }

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

      if (bline_bufp + bline_inc > bline_buf_mx) {
	i = 4;
	while (i-- > 0) {
	  if ((bline_buf=(char *)calloc(max_bline_b,sizeof(char)))!=NULL)
	    break;
	  fprintf(stderr," failure to allocate bline_buf(%d) %d\n",
		  max_sql,max_bline_b);
	  max_bline_b /= 2;
	  bline_inc /= 2;
	}
	if (bline_buf != NULL) {
	  bline_bufp = bline_buf;
	  bline_buf_mx = bline_buf+max_bline_b;
	}
	else {
	  s_abort("cannot allocate bline_buf ","");
	  exit(1);
	}
      }

      if (bline_bufp+bline_inc < bline_buf_mx ) {
	strncpy(bline_bufp,t_bline,bline_inc);
	ldes[m_seqnm].bline = bline_bufp;
	bline_bufp[bline_inc]= '\0';
	bline_bufp += bline_inc+1;
      }
      else {
	fprintf(stderr," bline_buf overrun\n");
      }

      ntt.entries++;		/* inc number of sequences */
      ntt.length += n1;	/* update total library length */
      if (ntt.length > LONG_MAX) {ntt.length -= LONG_MAX; ntt.carry++;}

#ifdef DEBUG
      /* This discovers most reasons for core dumps */
      if (pst.debug_lib)
	for (i=0; i<n1; i++)
	  if (aa1[i]>pst.nsq) {
	    fprintf(stderr,
		    "%s residue[%d/%d] %d range (%d) lcont/ocont: %d/%d\n%s\n",
		    qm_msg0.libstr,i,n1,aa1[i],pst.nsq,lcont,ocont,aa1ptr+i);
	    aa1[i]=0;
	    n1=i-1;
	    break;
	  }
#endif

      /* for ALTIVEC, must pad with 15 NULL's */
      for (i=0; i<SEQ_PAD+1; i++) {aa1ptr[n1+i]=0;}

      /* don't count long sequences more than once */
      if (aa1!=aa1ptr) {
	n1 += m_msg0.loff; m_msg0.db.entries--; ntt.entries--;
      }

      if (n1>1) {

	desptr = &ldes[m_seqnm];

	aa1i_arr[nseq] = (int)(aa1-aa1_buf);
	m_seqnm_arr[nseq] = m_seqnm;
	desptr->n1 = n1_arr[nseq] = n1;
	desptr->n1tot_p = n1tot_cur;
	desptr->lseek = lseek;
	desptr->loffset = loffset+l_off;
	desptr->cont = ocont;
	desptr->wrkr = node;
	desptr->nsfnum = nsfnum;
#ifdef SUPERFAMNUM
	if ((desptr->sfnum[0]=sfnum[0])>0 &&
	    (desptr->sfnum[1]=sfnum[1])>0 &&
	    (desptr->sfnum[2]=sfnum[2])>0 &&
	    (desptr->sfnum[3]=sfnum[3])>0 &&
	    (desptr->sfnum[4]=sfnum[4])>0 &&
	    (desptr->sfnum[5]=sfnum[5])>0 &&
	    (desptr->sfnum[6]=sfnum[6])>0 &&
	    (desptr->sfnum[7]=sfnum[7])>0 &&
	    (desptr->sfnum[8]=sfnum[8])>0 &&
	    (desptr->sfnum[9]=sfnum[9])>0) ;
#endif
	m_seqnm++;
	nseq++;

	if (m_seqnm >= max_sql) {
	  max_sql += MAXSQL;
	  if ((ldes=(struct sql *)realloc(ldes,max_sql*sizeof(struct sql)))
	      ==NULL) {
	    fprintf(stderr," failure to realloc ldes(%d) %ld\n",
		    max_sql,max_sql*sizeof(struct sql));
	    s_abort("cannot allocate ldes","");
	    exit(1);
	  }
	}

	/* increment ptrs */
	aa1prev = aa1;

	aa1 += n1+1+SEQ_PAD;
	ntbuff += n1+1+SEQ_PAD;

	/* if the buffer is filled */
	if (nseq >= max_buf_cnt || ntbuff >= m_msg0.pbuf_siz - m_msg0.maxn) {
	  /* provide filled buffer to workers */
#ifdef PVM_SRC
	  pvm_initsend(PvmDataRaw);
	  pvm_pkint(&nseq,1,1);
	  pvm_pkint(&ntbuff,1,1);
	  pvm_pkint(n1_arr,nseq,1);
	  pvm_pkint(aa1i_arr,nseq,1);
	  pvm_pkint(m_seqnm_arr,nseq,1);
	  pvm_send(pinums[node],STARTTYPE4);

	  pvm_initsend(PvmDataRaw);
	  pvm_pkbyte((char *)aa1_buf,ntbuff,1);
	  pvm_send(pinums[node],STARTTYPE5);
#endif
#ifdef MPI_SRC
	  MPI_Send(&nseq,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
	  MPI_Send(&ntbuff,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
	  MPI_Send(n1_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
	  MPI_Send(aa1i_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
	  MPI_Send(m_seqnm_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);

	  MPI_Send(aa1_buf,ntbuff,MPI_BYTE,node,STARTTYPE5,MPI_COMM_WORLD);
#endif
	  nseq =  0;

	  aa1 = aa1_buf;
	  ntbuff = 0;
	  if (++node >= nnodes) node = FIRSTNODE;
	}

      loop1:
	if (lcont) {
	  memcpy(aa1,&aa1prev[n1-m_msg0.loff],m_msg0.loff);
	  aa1ptr = &aa1[m_msg0.loff];
	  ocont = lcont;
	  maxt = m_msg0.maxt3;
	  loffset += n1 - m_msg0.loff;
	}
	else {
	  if (ocont) *n1tot_cur = n1tot_v;
	  n1tot_v = 0;
	  n1tot_cur = NULL;

	  ocont = 0;
	  aa1ptr = aa1;
	  maxt = m_msg0.maxn;
	  loffset = 0l;
	}
      }
    }
  }	/* for (iln < nln) */

  if (nseq > 0) {
#ifdef PVM_SRC
    pvm_initsend(PvmDataRaw);
    pvm_pkint(&nseq,1,1);
    pvm_pkint(&ntbuff,1,1);
    pvm_pkint(n1_arr,nseq,1);
    pvm_pkint(aa1i_arr,nseq,1);
    pvm_pkint(m_seqnm_arr,nseq,1);
    pvm_send(pinums[node],STARTTYPE4);

    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)aa1_buf,ntbuff,1);
    pvm_send(pinums[node],STARTTYPE5);
#endif
#ifdef MPI_SRC
    MPI_Send(&nseq,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(&ntbuff,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(n1_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(aa1i_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(m_seqnm_arr,nseq,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);

    MPI_Send(aa1_buf,ntbuff,MPI_BYTE,node,STARTTYPE5,MPI_COMM_WORLD);
#endif
  }

  /*   fprintf(stderr," all sequences sent\n"); */

  if (ntt.entries <= 0) {
    s_abort("no reference library sequences found\n","");
  }

  zero = 0;
  for (node=FIRSTNODE; node < nnodes; node++) {
#ifdef PVM_SRC
    pvm_initsend(PvmDataRaw);
    pvm_pkint(&zero,1,1);
    pvm_pkint(&zero,1,1);
    pvm_pkint(n1_arr,1,1);
    pvm_pkint(aa1i_arr,1,1);
    pvm_pkint(m_seqnm_arr,1,1);
    pvm_send(pinums[node],STARTTYPE4);
#endif
#ifdef MPI_SRC
    MPI_Send(&zero,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(&zero,1,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(n1_arr,0,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(aa1i_arr,0,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
    MPI_Send(m_seqnm_arr,0,MPI_INT,node,STARTTYPE4, MPI_COMM_WORLD);
#endif
  }

  for (node = FIRSTNODE; node < nnodes; node++) {
#ifdef PVM_SRC
    bufid = pvm_recv(-1,STARTTYPE0);
    pvm_bufinfo(bufid,NULL,NULL,&tid);
    snode = tidtonode(tid);
    pvm_upkint(&lcnt,1,1);
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv(&lcnt,1,MPI_INT,MPI_ANY_SOURCE,STARTTYPE0,
	     MPI_COMM_WORLD,&mpi_status);
    snode= mpi_status.MPI_SOURCE;
#endif
    wlsn [snode-FIRSTNODE] = lcnt;
    fprintf(stderr," %d sequences at %d\n",lcnt,snode);
  }

  /* print out all descriptions */
  /*
  for (node = FIRSTNODE; node < nnodes; node++)
    for (lcnt = 0; lcnt < wlsn[node-FIRSTNODE]; lcnt ++)
      printf("%2d:%3d\t%s\n",node,lcnt,ldes[lcnt].bline);
  */

  /* Calculate cumulative totals and send to workers for a self search */

  clsn [0] = nclib= 0;
  for (node = FIRSTNODE; node < nnodes-1; node++) {
    /* clsn[] is for the next node */
    clsn[node-FIRSTNODE+1] = nclib += wlsn[node-FIRSTNODE];
  }

  if (m_msg0.self)
    for (node = FIRSTNODE; node < nnodes; node++) {
#ifdef PVM_SRC
      pvm_initsend(PvmDataRaw);
      pvm_pkint(&clsn[node-FIRSTNODE],1,1);
      pvm_send(pinums[node],STARTTYPE1);
#endif
#ifdef MPI_SRC
    MPI_Send(&clsn[node-FIRSTNODE],1,MPI_INT,node,STARTTYPE1,MPI_COMM_WORLD);
#endif
      fprintf(stderr,"sending lend: %d to worker %d\n",clsn[node-FIRSTNODE],node);
    }

  last_msg_b[0] = m_msg0.nbr_seq = m_msg1.nbr_seq = ntt.entries;

  qres_bufsize = BFR;
  /* if BFR is too big for this library, reduce it */
  while ( ntt.entries*(m_msg0.nitt1+1)/(2*nnodes) < qres_bufsize) {
    qres_bufsize /= 2;
    if ((qres_bufsize%(m_msg0.nitt1+1))!= 0) {
      qres_bufsize *= (m_msg0.nitt1+1);
      break;
    }
    if (qres_bufsize < 50) break;
  }
  last_msg_b[1] = qres_bufsize;

  fprintf(stderr," using BFR=%d/%d\n",qres_bufsize,BFR);

#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkint(last_msg_b,2,1);
  for (node=FIRSTNODE; node < nnodes; node++)
    pvm_send(pinums[node],STARTTYPE0);
#endif
#ifdef MPI_SRC
  for (node=FIRSTNODE; node < nnodes; node++)
    MPI_Send(last_msg_b,2,MPI_INT,node,STARTTYPE0,MPI_COMM_WORLD);
#endif  

  tscan = tprev = s_time();

/**************************************
  The logic of this section has been simplified to allow multistage
  comparison functions to be used and alignments to be generated.

  	send 1st query to workers
	get next query sequence from host (m_msp1)
    L1:	get results from next-1 search (m_msp0)
	sort the results of the next-1 search
	(possibly) do additional stages of search
	(possibly produce alignments for search
	send next query to workers (m_msp1)
	display result of next-1 search (m_msp0)
	get next query sequence from host (m_msp1)
	goto L1;

As a result of the interleaving, there must be two qm_msg structures,
one for the next-1 sequence (which is required for labeling the
output), and one for the next sequence (which is sent to the workers
while the results are being displayed.  qm_msp0 and qm_msp1 alternate
between these two structures.
***************************************/

/*
  qm_msp0 points to the older qm_msg 
  qm_msp1 points to the newer qm_msg
  the assignment below goes with curtype==ONETYPE
*/
  m_msp0 = &m_msg0;
  m_msp1 = &m_msg1;

  qm_msp0 = &qm_msg0;
  qm_msp1 = &qm_msg1;
  
  aa0p0 = aa00;		/* aa0p0 is the "old" sequence */
  aa0p1 = aa01;		/* aa0p1 is the "new" sequence */

  last_params(aa00,m_msp0->n0,m_msp0,&pst,qm_msp0);

  /* process_hist() is called here to get find_zp(), and some other
     structures initialized that would otherwise not be initialized
     because z-scores are not being calculated */

  if (m_msp0->escore_flg) {
     pst.zsflag_f = process_hist(stats,nstats,*m_msp0,pst,
				 &m_msp0->hist,&m_msp0->pstat_void,0);
     stats_done=1;
  }

  if (m_msp0->qshuffle && qstats==NULL) {
    if ((qstats =
	 (struct stat_str *)calloc(m_msg0.shuff_max+1,sizeof(struct stat_str)))==NULL)
      s_abort ("Cannot allocate qstats struct","");
  }
  nqstats = 0;

/* Send first query sequence to each worker */

  if (m_msg0.dfile[0] && (fdata=fopen(m_msg0.dfile,"w"))!=NULL)
    fprintf(fdata,"%3d>%-50s\n",qlib,qm_msp0->libstr);

#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkbyte((char *)qm_msp0,sizeof(qm_msg0),1);
  if (qm_msp0->n0 > 0) {
    pvm_pkbyte((char *)aa0p0,qm_msp0->n0+1+SEQ_PAD,1);
    if (m_msg0.ann_flg) pvm_pkbyte((char *)m_msp0->aa0a,qm_msp0->n0+1,1);
  }
  for (node = FIRSTNODE; node < nnodes; node++)
    pvm_send(pinums[node],MSEQTYPE);
#endif
#ifdef MPI_SRC
  for (node = FIRSTNODE; node < nnodes; node++) {
    MPI_Send(qm_msp0,sizeof(qm_msg0),MPI_BYTE,node,MSEQTYPE,MPI_COMM_WORLD);
    if (qm_msp0->n0 > 0) {
      MPI_Send(aa0p0,qm_msp0->n0+1+SEQ_PAD,MPI_BYTE,node,
	       MSEQTYPE1,MPI_COMM_WORLD);
      if (m_msg0.ann_flg) {
	if (m_msp0->aa0a == NULL) {
	  fprintf(stderr," m_msp0: %o/%oaa0a is null\n",m_msp0,m_msp0->aa0a);
	}
	MPI_Send(m_msp0->aa0a,qm_msp0->n0+1,MPI_BYTE,node, MSEQTYPE2,MPI_COMM_WORLD);
      }
    }
  }
#endif

  /* Get second query sequence (additional query sequences are read in
     the main loop */

  m_msp1->n0 = qm_msp1->n0 = 
    QGETLIB(aa0p1,MAXTST,q_bline, sizeof(q_bline),&qseek, &qlcont,q_file_p,&m_msp1->sq0off);
  strncpy(qm_msp1->libstr,q_bline,sizeof(qm_msg0.libstr)-20);
  qm_msp1->libstr[sizeof(qm_msg0.libstr)-21]='\0';
  if ((bp=strchr(qm_msp1->libstr,' '))!=NULL) *bp='\0';

  /* if annotations are included in sequence, remove them */
  if (m_msg0.ann_flg) {
    m_msp1->n0 = qm_msp1->n0 = 
      ann_scan(aa0p1,qm_msp1->n0,m_msp1,m_msp1->qdnaseq);
#ifdef DEBUG
    fprintf(stderr,"m_msp1->/aa0a is: %o/%o\n",m_msp1,m_msp1->aa0a);
#endif
  }

  if (qm_msp1->n0 > 0 && m_msg0.term_code && !qlcont && 
      m_msg0.qdnaseq == SEQT_PROT &&
      aa0p1[m_msp1->n0-1]!=m_msg0.term_code) {
    aa0p1[m_msp1->n0++]=m_msg0.term_code;
    aa0p1[m_msp1->n0]=0;
    qm_msp1->n0 = m_msp1->n0;
  }

  /* for ALTIVEC, must pad with 15 NULL's */
  if (m_msp1->n0 > 0) {
    for (i=0; i<SEQ_PAD+1; i++) {aa0p1[m_msp1->n0+i]=0;}
  }

  qm_msp1->slist = 0;
  qm_msp1->seqnm = qlib;

  last_params(aa0p1,m_msp1->n0,m_msp1,&pst,qm_msp1);

  sprintf(tmp_str," - %d %s", qm_msp1->n0, q_sqnam);
  if (strlen(qm_msp1->libstr) + strlen(tmp_str) >= sizeof(qm_msg0.libstr))
    qm_msp1->libstr[sizeof(qm_msg0.libstr)-strlen(tmp_str)-2] = '\0';
  strncat(qm_msp1->libstr,tmp_str,
	  sizeof(qm_msg0.libstr)-strlen(qm_msp1->libstr)-1);
  qm_msp1->libstr[sizeof(qm_msg0.libstr)-1]='\0';

  naa0 = 0;  /* reset node counter */

  /* sit in loop and collect results */
  nbest = nopt = 0;
  zbestcut = -BIGNUM;


  while (1) {

#ifdef PVM_SRC
    bufid = pvm_recv(-1,curtype);
    pvm_bufinfo(bufid,NULL,NULL,&tid);
    pvm_upkbyte((char *)&bestr[0],sizeof(struct comstr)*(qres_bufsize+1),1);
    snode = tidtonode(tid);
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv(bestr,sizeof(struct comstr)*(qres_bufsize+1),
	     MPI_BYTE,MPI_ANY_SOURCE,curtype,MPI_COMM_WORLD,&mpi_status);
    snode = mpi_status.MPI_SOURCE;
#endif

    nres = bestr[qres_bufsize].seqnm & ~FINISHED;

#ifdef DEBUG
    fprintf(stderr,"%d results from %d\n",nres,snode);
#endif

    if (bestr[qres_bufsize].seqnm&FINISHED) {	/* a worker is finished */
      naa0++;

      /* fast_flag == 1 => send new sequences immediately */
      fast_flag = ((m_msp0->stages==1) && !(m_msp0->markx & MX_M9SUMM) &&
		   (m_msp0->ashow == 0) && (m_msp0->last_calc_flg==0));
      /* send a new query sequence if no more processing required */
      if (fast_flag) {
#ifdef PVM_SRC
	pvm_initsend(PvmDataRaw);
	pvm_pkbyte((char *)qm_msp1,sizeof(qm_msg1),1);
	if (qm_msp1->n0 != -1) {
	  pvm_pkbyte((char *)aa0p1,qm_msp1->n0+1+SEQ_PAD,1);
	  if (m_msp1->ann_flg) pvm_pkbyte((char *)m_msp1->aa0a,qm_msp1->n0+1,1);
	}
	pvm_send(tid,MSEQTYPE);
#endif
#ifdef MPI_SRC
	MPI_Send(qm_msp1,sizeof(qm_msg1),MPI_BYTE,snode,MSEQTYPE,MPI_COMM_WORLD);
	if (qm_msp1->n0 != -1) {
	  MPI_Send(aa0p1,qm_msp1->n0+1+SEQ_PAD,MPI_BYTE,snode,MSEQTYPE1,MPI_COMM_WORLD);
	  if (m_msp1->ann_flg)
	    MPI_Send(m_msp1->aa0a,qm_msp1->n0+1,MPI_BYTE,snode,MSEQTYPE2,MPI_COMM_WORLD);
	}
#endif
      }
    }

#ifdef DEBUG
    if (pst.debug_lib)
      fprintf(stderr," unpacking %d from %d; nbest %d\n",nres,snode,nbest);
#endif

    /* this section is now more complex because can get groups of
       sequence results; e.g. forward and reverse frame */

    t_best = t_rbest = t_qrbest = -1;
    tm_escore = t_rescore = t_qrescore = FLT_MAX;
    for (ires = 0; ires < nres; ires++) {
      desptr = &ldes[bestr[ires].m_seqnm];

      /* save raw results */
      if (fdata) {
	strncpy(tlibstr,desptr->bline,10);
	if ((bp=strchr(tlibstr,' '))!=NULL) *bp='\0';
	fprintf(fdata,"%-10s\t%4d\t%4d\t%d\t%4d\t%4d\t%4d\t%8ld\n",
		tlibstr,desptr->sfnum[0],desptr->n1,bestr[ires].frame,
		bestr[ires].score[0],bestr[ires].score[1],bestr[ires].score[2],
		desptr->lseek);
      }

      i_score = bestr[ires].score[pst.score_ix];
      e_score = bestr[ires].escore;
      k_comp = bestr[ires].comp;
      k_H = bestr[ires].H;

      t_n1 = desptr->n1;
      if (i_score > t_best) {tm_best = t_best = i_score;}
      if (e_score < tm_escore) tm_escore = e_score;

      if (m_msp0->qshuffle) {
	if (bestr[ires].qr_score > t_qrbest) 
	  t_qrbest = bestr[ires].qr_score;
	if (bestr[ires].qr_escore < t_qrescore) 
	  t_qrescore = bestr[ires].qr_escore;

	if (bestr[ires].frame==m_msp0->nitt1 && 
	    nqstats < m_msp0->shuff_max &&
	    bestr[ires].qr_score >= 0) {
	  qstats[nqstats].n1 = t_n1;	/* save the best score */
	  qstats[nqstats].comp =  bestr[ires].comp;
	  qstats[nqstats].H = bestr[ires].H;
	  qstats[nqstats].escore = t_qrescore;
	  qstats[nqstats++].score = t_qrbest;
	  t_qrbest = -1;	/* reset t_qrbest, t_qrescore */
	  t_qrescore = FLT_MAX;
	}
      }

      if (pst.zsflag >= 10 && bestr[ires].r_score > t_rbest) {
	t_rbest = bestr[ires].r_score;
	t_rescore = bestr[ires].r_escore;
      }

      if (nstats < MAXSTATS) {
	if (bestr[ires].frame == m_msg0.nitt1) {
	  stats[nstats].n1 = t_n1;
	  stats[nstats].comp = k_comp;
	  stats[nstats].H = k_H;

	  if (pst.zsflag > 10) {
	    tm_best = t_rbest;
	    tm_escore = t_rescore;
	    t_rbest = -1;
	    t_rescore = FLT_MAX;
	  }
	  stats[nstats].escore = tm_escore;
	  stats[nstats++].score = tm_best;
	  tm_escore = FLT_MAX;
	  t_best = -1;
	}
      }
      else if (pst.zsflag >=0) {	/* nstats >= MAXSTATS, zsflag >=0 */
	if (!stats_done ) {
	  pst.n0 = qm_msp0->n0;
	  pst.zsflag_f = process_hist(stats,nstats,*m_msp0,pst,
				      &m_msp0->hist, &m_msp0->pstat_void,0);
	  stats_done = 1;
	  kstats = nstats;
	  for (i=0; i<nbest; i++) {
	    bptr[i]->zscore = (*find_zp)(bptr[i]->score[pst.score_ix],
					 bptr[i]->escore,bptr[i]->n1,
					 bptr[i]->comp, m_msp0->pstat_void);
	  }
	}
#ifdef SAMP_STATS
	if (!m_msp0->escore_flg) {
	  jstats = nrand(kstats++);
	  if (jstats < MAXSTATS) {
	    stats[jstats].n1 = t_n1;	/* save the best score */
	    stats[jstats].comp =  k_comp;
	    stats[jstats].H = k_H;
	    if (pst.zsflag >=10) t_best = t_rbest;
	    stats[jstats].score = t_best;
	  }
	}
#endif
      }

      if (stats_done) {
	zscore=(*find_zp)(i_score,e_score,desptr->n1,k_comp,
			  m_msp0->pstat_void);
	if (bestr[ires].frame == m_msg0.nitt1) {
	  addhistz((*find_zp)(tm_best,tm_escore,t_n1,k_comp,
			      m_msp0->pstat_void),
		   &(m_msp0->hist));
	  t_best = t_rbest = -1;
	}

      }
      else zscore = (double) i_score;

      if (zscore > zbestcut) {
	if (nbest>=MAXBEST) {
	  selectbestz(bptr, nbest-MAXBEST/4-1, nbest); 
	  nbest -= MAXBEST/4;
	  zbestcut = bptr[nbest-1]->zscore;
	  best_flag = 0;
	}
	/* if zbestcut == -BIGNUM, bptr[] has not been reinitialized */
	else if (best_flag) bptr[nbest]=&best[nbest];

	bptr[nbest]->m_seqnm  = bestr[ires].m_seqnm ;
	bptr[nbest]->seqnm  = bestr[ires].seqnm;
	bptr[nbest]->score[0] = bestr[ires].score[0];
	bptr[nbest]->score[1] = bestr[ires].score[1];
	bptr[nbest]->score[2] = bestr[ires].score[2];
	bptr[nbest]->escore = bestr[ires].escore;
	bptr[nbest]->segnum = bestr[ires].segnum;
	bptr[nbest]->seglen = bestr[ires].seglen;
	bptr[nbest]->comp = bestr[ires].comp;
	bptr[nbest]->H = bestr[ires].H;
	bptr[nbest]->zscore = zscore;
	bptr[nbest]->wrkr   = snode;
	bptr[nbest]->desptr = desptr;
	bptr[nbest]->lseek = desptr->lseek; /* needed for identifying alternate
					    strand scores from same sequence */
	bptr[nbest]->n1 = desptr->n1;
	bptr[nbest]->frame = bestr[ires].frame;

	/*	this was used when -m 9 info was calculated in 1st scan */
	/*
	bptr[nbest]->sw_score = bestr[ires].sw_score;
	if (bestr[ires].sw_score > -1) {
	  nopt++;
	  bptr[nbest]->a_len = bestr[ires].a_len;
	  bptr[nbest]->percent = bestr[ires].percent;
	  bptr[nbest]->gpercent = bestr[ires].gpercent;
	  bptr[nbest]->min0 = bestr[ires].min0;
	  bptr[nbest]->min1 = bestr[ires].min1;
	  bptr[nbest]->max0 = bestr[ires].max0;
	  bptr[nbest]->max1 = bestr[ires].max1;
	  bptr[nbest]->ngap_q = bestr[ires].ngap_q;
	  bptr[nbest]->ngap_l = bestr[ires].ngap_l;
	}
	else {
	  bptr[nbest]->percent = -1.0;
	  bptr[nbest]->min0 = bptr[nbest]->min1 = bptr[nbest]->max0 = 
	    bptr[nbest]->max1 = 0;
	}
	*/

	nbest++;
      }
    }	/* for loop */
    if (naa0 < nnodes-FIRSTNODE) continue;

    gstring2[0]='\0';

    /* get gstring2,3 - algorithm/parameter description */
#ifdef PVM_SRC
    bufid = pvm_recv(pinums[FIRSTNODE],PARAMTYPE);
    pvm_upkbyte(gstring2,sizeof(gstring2),1);
    pvm_upkbyte(gstring3,sizeof(gstring3),1);
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv(gstring2,sizeof(gstring2),MPI_BYTE,FIRSTNODE,PARAMTYPE,
	     MPI_COMM_WORLD,&mpi_status);
    MPI_Recv(gstring3,sizeof(gstring3),MPI_BYTE,FIRSTNODE,PARAMTYPE,
	     MPI_COMM_WORLD,&mpi_status);
#endif

/* ********************** */
/* analyze the results    */
/* ********************** */
    
    if (!stats_done) {
      if (nbest < 20 || pst.zsflag <= 0) {
	pst.zsflag_f = -1;
      }
      else {
	pst.n0 = qm_msp0->n0;
	pst.zsflag_f = process_hist(stats,nstats,*m_msp0,pst,
				    &m_msp0->hist, &m_msp0->pstat_void,stats_done);

	for (i=0; i<nbest; i++)
	  bptr[i]->zscore = (*find_zp)(bptr[i]->score[pst.score_ix],
				       bptr[i]->escore, bptr[i]->n1,
				       bptr[i]->comp, m_msp0->pstat_void);
      }
    }

    m_msp0->db.entries = ntt.entries;
    m_msp0->db.length = ntt.length;
    m_msp0->db.carry = ntt.carry;

    if (pst.zdb_size < 1) pst.zdb_size = ntt.entries;

    if (!qm_msp0->qshuffle) {
      last_stats(aa0p0, m_msp0->n0,
		 stats,nstats, bptr,nbest, *m_msp0, pst, 
		 &m_msp0->hist, &m_msp0->pstat_void);
    }
    else {
      last_stats(aa0p0, m_msp0->n0,
		 qstats,nqstats, bptr,nbest, *m_msp0, pst, 
		 &m_msp0->hist, &m_msp0->pstat_void);
    }

    if (m_msp0->last_calc_flg) {
      nbest = last_calc(bptr,nbest, *m_msp0, &pst,qm_msp0,
			m_msp0->pstat_void);
    }

    sortbeste(bptr,nbest);
    scale_scores(bptr,nbest,m_msp0->db,pst,m_msp0->pstat_void);

    if (pst.zsflag >= 0 && bptr[0]->escore >= m_msg0.e_cut) goto no_results;

    /*      else sortorder(bptr,nbest,wlsn,nnodes); */

/* if more than one stage or markx==9, calculate opt scores or do alignment */
/* send results to workers as available */

    if (m_msg0.stages > 1 || m_msg0.markx & MX_M9SUMM) {

      /* to determine how many sequences to re-align (either for
	 do_opt() or calc_id() we need to modify m_msg.mshow to get
	 the correct number of alignments */

      if (m_msg0.mshow_flg != 1 && pst.zsflag >= 0) {
	for (i=0; i<nbest && bptr[i]->escore< m_msg0.e_cut; i++) {}
	m_msg0.mshow = i;
      }

      /* allocate space for a_struct info */
      if (m_msg0.markx & MX_M9SUMM && m_msg0.mshow > 0) {
	if ((aln_d_base=(struct a_struct *)
	     calloc((size_t)m_msg0.mshow,sizeof(struct a_struct)))==NULL) {
	  fprintf(stderr," cannot allocate a_struct %d\n", m_msg0.mshow);
	  exit(1);
	}

	for (is = 0; is < m_msg0.mshow; is++ ) {
	  bptr[is]->aln_d = &aln_d_base[is];
	}
      }

      do_stage2(bptr,m_msg0.mshow, *m_msp0, DO_OPT_FLG, qm_msp0);
    }

  no_results:
    tdone = s_time();
    tddone = time(NULL);

    /* changed from >> to >>> because qm_msp0->libstr is missing '>' */
    fprintf (outfd, "%3d>>>%s\n", qlib,qm_msp0->libstr);

    /* make certain that m_msp0->n0, libstr are current */
    m_msp0->n0 = qm_msp0->n0;
    /*    strncpy(m_msp0->libstr,qm_msp0->libstr,sizeof(m_msg0.libstr)); */

    prhist (outfd,*m_msp0,pst,m_msp0->hist,nstats,m_msp0->db,gstring2);

    if (bptr[0]->escore < m_msg0.e_cut) {

      showbest (outfd, bptr, nbest, qlib, m_msp0,pst,ntt,gstring2);

      if (m_msg0.markx & MX_M9SUMM) {
	fprintf(outfd,"\n>>>%s#%d %s%s, %d %s vs %s library\n",
		m_msg0.tname,qlib,qm_msp0->libstr,
		(m_msg0.revcomp ? "-":"\0"), qm_msp0->n0, m_msg0.sqnam,
		m_msg0.lname);
      }
      else if (m_msg0.markx & MX_M10FORM) {
	if ((bp=strchr(qm_msp0->libstr,' '))!=NULL) *bp = '\0';
	fprintf(outfd,"\n>>>%s#%d %s%s, %d %s vs %s library\n",
		m_msg0.tname,qlib,qm_msp0->libstr,
		(m_msg0.revcomp ? "-":"\0"), qm_msp0->n0, m_msg0.sqnam,
		m_msg0.lname);
	if (bp!=NULL) *bp=' ';
	fprintf(outfd,"; mp_name: %s\n",argv[0]);
	fprintf(outfd,"; mp_ver: %s\n",mp_verstr);
	fprintf(outfd,"; mp_argv:");
	for (i=0; i<argc; i++)
	  fprintf(outfd," %s",argv[i]);
	fputc('\n',outfd);
	fputs(gstring3,outfd);
	fputs(hstring1,outfd);
      }

      /* ashow is -1 if not set, -d 0 indicates no alignments, > 0 if set */
      /* if ashow is -1, m_msg.nshow (set by e_cut above) sets limit
	 in showalign */
      
      if (m_msp0->ashow != 0) {
	/* showalign needs m_msp->qtitle, so fill it in */
	strncpy(m_msp0->qtitle,qm_msp0->libstr,MAX_FN-1);
	m_msp0->qtitle[MAX_FN-1]='\0';
	showalign (outfd, bptr, nbest, qlib, *m_msp0, pst, gstring2);
      }
    }
    else {
      if (m_msg0.markx & (MX_M9SUMM + MX_M10FORM)) {
	fprintf(outfd,"\n>>>%s#%d %s%s, %d %s vs %s library\n",
		m_msg0.tname,qlib,qm_msp0->libstr,(m_msg0.revcomp ? "-":"\0"), qm_msg0.n0, m_msg0.sqnam,
		m_msg0.lname);
	fprintf(outfd,">>>!!! No sequences with E() < %f\n",m_msg0.e_cut);
      }
      else fprintf(outfd,"!! No sequences with E() < %f\n",m_msg0.e_cut);
    }

    if (! (m_msg0.markx & (MX_M9SUMM + MX_M10FORM))) {
      fprintf(outfd,"/** search time: ");
      ptime(outfd,tdone-tprev);
      fprintf(outfd," **/\n");
      tprev = tdone;
    }
    else if (m_msg0.markx & MX_M9SUMM) {
      if (aln_d_base != NULL) {
	free((void *)aln_d_base);
	aln_d_base = NULL;
      }
      fprintf(outfd,">>>***\n");
      fprintf(outfd,"/** %s **/\n",gstring2);
      fprintf(outfd,"/** %s **/\n",m_msp0->hist.stat_info);
      fprintf(outfd,">>><<<\n");
    }
    else if (m_msg0.markx & MX_M10FORM) {
      fprintf(outfd,">>><<<\n");
    }
    fflush(outfd);
    
/* *********************** */
/* end of analysis/display */
/* *********************** */


/* *********************** */
/* start the next search   */                                           
/* *********************** */

    if (fdata) {		/* label the results file */
      fprintf(fdata,"/** %s **/\n",gstring2);
      fprintf(fdata,"%3d>%-50s\n",qlib-1,qm_msp1->libstr);
      fflush(fdata);
    }
    
    if (m_msp1->escore_flg) {	/* re-initialize some stats stuff before search */
      pst.zsflag_f = process_hist(stats,nstats,*m_msp1,pst,
				  &m_msp1->hist,&m_msp1->pstat_void,0);
      stats_done=1;
    }
    else stats_done = 0;

    /* set up qstats if necessary - different queries have different qshuffle */
    if (m_msp1->qshuffle && qstats==NULL) {
      if ((qstats =
	   (struct stat_str *)calloc(m_msg0.shuff_max+1,sizeof(struct stat_str)))==NULL)
	s_abort ("Cannot allocate qstats struct","");
    }

    nqstats = nstats = 0;

    /* send new qm_msp, sequence */
    if (!fast_flag) {
#ifdef PVM_SRC
      pvm_initsend(PvmDataRaw);
      pvm_pkbyte((char *)qm_msp1,sizeof(qm_msg1),1);
      if (qm_msp1->n0 != -1) {
	pvm_pkbyte((char *)aa0p1,qm_msp1->n0+1+SEQ_PAD,1);
	if (m_msp1->ann_flg) {
	  pvm_pkbyte((char *)m_msp1->aa0a,qm_msp1->n0+1,1);
	}	  
      }
      for (node = FIRSTNODE; node < nnodes; node++)
	pvm_send(pinums[node],MSEQTYPE);
#endif
#ifdef MPI_SRC
      for (node=FIRSTNODE; node < nnodes; node++) {
	MPI_Send(qm_msp1,sizeof(qm_msg1),MPI_BYTE,node,MSEQTYPE,
		 MPI_COMM_WORLD);
	if (qm_msp1->n0 != -1) {
	  MPI_Send(aa0p1,qm_msp1->n0+1+SEQ_PAD,MPI_BYTE,node,MSEQTYPE1,MPI_COMM_WORLD);
	  if (m_msp1->ann_flg)
	    MPI_Send(m_msp1->aa0a,qm_msp1->n0+1,MPI_BYTE,snode,MSEQTYPE2,MPI_COMM_WORLD);
	}
      }
#endif
    }

    qlib++;
    if (qm_msp1->n0 != -1) {
      qtt.entries++;
      qtt.length += qm_msp1->n0;
    }
    else goto done;
    
/* ******************************** */
/* flip m_msg, qm_msg, aa0 pointers */
/* ******************************** */

    naa0 = 0;
    best_flag = 1;
    nbest = nopt = 0;
    zbestcut = -BIGNUM;
    if (curtype == ONETYPE) {
      curtype = TWOTYPE;
      qm_msp0 = &qm_msg1;
      qm_msp1 = &qm_msg0;
      m_msp0 = &m_msg1;
      m_msp1 = &m_msg0;
      aa0p0 = aa01;
      aa0p1 = aa00;
    }
    else  {
      curtype = ONETYPE;
      qm_msp0 = &qm_msg0;
      qm_msp1 = &qm_msg1;
      m_msp0 = &m_msg0;
      m_msp1 = &m_msg1;
      aa0p0 = aa00;
      aa0p1 = aa01;
    }


/* **********************************************************/
/* all library sequences are done get next library sequence */
/* **********************************************************/

    m_msp1->n0 = qm_msp1->n0 = 
      QGETLIB(aa0p1,MAXTST,q_bline, sizeof(q_bline),&qseek, &qlcont,q_file_p,&m_msp1->sq0off);
    strncpy(qm_msp1->libstr,q_bline,sizeof(qm_msg0.libstr)-20);
    qm_msp1->libstr[sizeof(qm_msg0.libstr)-21]='\0';

    if ((qlib+1) >= m_msg0.ql_stop) { qm_msp1->n0 = m_msp1->n0 = -1;}

    if (qm_msp1->n0 > 0 && m_msg0.term_code && !qlcont &&
	m_msg0.qdnaseq==SEQT_PROT &&
	aa0p1[m_msp1->n0-1]!=m_msg0.term_code) {
      aa0p1[m_msp1->n0++]=m_msg0.term_code;
      aa0p1[m_msp1->n0]=0;
      qm_msp1->n0 = m_msp1->n0;
    }

    /* for ALTIVEC, must pad with 15 NULL's */
    if (m_msg0.n0 > 0) {
       for (i=0; i<SEQ_PAD+1; i++) {aa00[m_msg0.n0+i]=0;}
    }

    qm_msp1->slist = 0;
    /*
    leng = strlen (qm_msp1->libstr);
    sprintf (&(qm_msp1->libstr[leng]), " %d %s", qm_msp1->n0, q_sqnam);
    */
    sprintf(tmp_str," %d %s", qm_msp1->n0, q_sqnam);
    if (strlen(qm_msp1->libstr) + strlen(tmp_str) >= sizeof(qm_msg0.libstr))
      qm_msp1->libstr[sizeof(qm_msg0.libstr)-strlen(tmp_str)-2] = '\0';
    strncat(qm_msp1->libstr,tmp_str,
	    sizeof(qm_msg0.libstr)-strlen(qm_msp1->libstr)-1);
    qm_msp1->libstr[sizeof(qm_msg0.libstr)-1]='\0';

    qm_msp1->seqnm = qlib;

    last_params(aa0p1,m_msp1->n0,m_msp1,&pst,qm_msp1);

  }	    /* while loop */
  
  /* ******************** */
  /* end of library while */
  /* ******************** */

 done:
  tdone = s_time();
  if (m_msg0.markx & (MX_M9SUMM + MX_M10FORM)) fputs(">>>///\n",outfd);
  printsum(outfd);
  if (outfd!=stdout) printsum(stdout);
  printsum(stderr);
#ifdef PVM_SRC
  pvm_exit();
#endif
#ifdef MPI_SRC
  MPI_Finalize();
#endif

  exit(0);
}   /* End of main program */

void
printsum(FILE *fd)
{
  double db_tt;
  char tstr1[26], tstr2[26];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */
  if (qtt.carry==0) {
    fprintf(fd, "\n%ld residues in %d query   sequences\n", qtt.length, qtt.entries);
  }
  else {
    db_tt = (double)qtt.carry*(double)LONG_MAX + (double)qtt.length;
    fprintf(fd, "\n%.0g residues in %d query   sequences\n", db_tt, qtt.entries);
  }

  if (ntt.carry==0) {
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt.length, ntt.entries);
  }
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.6f residues in %ld library sequences\n", db_tt, ntt.entries);
  }

  fprintf(fd," %d processors (%d workers) were used\n",
	  nnodes+-FIRSTNODE+1,nnodes-FIRSTNODE);
  fprintf(fd," Pvcomplib [%s]\n start: %s done: %s\n",mp_verstr,tstr1,tstr2);
  fprintf(fd," Loading time: ");
  ptime(fd, tscan - tstart);
  fprintf (fd," Scan time: ");
  ptime (fd, tdone - tscan);
  fprintf (fd,"\n");
  fprintf (fd, "\nFunction used was %s [%s]\n", prog_func,verstr);
}

void fsigint()
{
  int i;

  tdone = s_time();
  tddone = time(NULL);

  if (outfd!=stdout) fprintf(outfd,"/*** interrupted ***/\n");
  fprintf(stderr,"/*** interrupted ***/\n");

  printsum(stdout);
  if (outfd!=stdout) printsum(outfd);

#ifdef PVM_SRC
  for (i=FIRSTNODE; i<nnodes; i++) pvm_kill(pinums[i]);
  pvm_exit();
#endif
#ifdef MPI_SRC
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif
  exit(1);
}
