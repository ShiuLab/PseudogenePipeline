/*	initfa.c	*/

/* $Name: fa_34_26_5 $ - $Id: initfa.c,v 1.148 2007/04/26 18:40:58 wrp Exp $ */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */

/* init??.c files provide function specific initializations */

/* h_init()	- called from comp_lib.c, comp_thr.c to initialize pstruct ppst
   		  which includes the alphabet, and pam matrix

   alloc_pam()	- allocate pam matrix space
   init_pam2()	- convert from 1D to 2D pam

   init_pamx()	- convert from 1D to 2D pam

   f_initenv()	- set up mngmsg and pstruct defaults
   f_getopt()	- read fasta specific command line options
   f_getarg()	- read ktup

   resetp()	- reset the parameters, scoring matrix for DNA-DNA/DNA-prot

   query_parm()	- ask for ktup
   last_init()	- some things must be done last

   f_initpam()	- set some parameters based on the pam matrix
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "defs.h"
#include "structs.h"
#include "param.h"

#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#define XTERNAL
#include "upam.h"
#include "uascii.h"
#undef XTERNAL

#define MAXWINDOW 32

int initpam(char *, struct pstruct *);
void init_pam2 (struct pstruct *ppst);
void extend_pssm(unsigned char *aa0, int n0, struct pstruct *ppst);
void build_xascii(int *qascii, char *save_str);
void ann_ascii(int *qascii, char *ann_arr);
void re_ascii(int *qascii, int *pascii);
extern int nrand(int);

/*  at some point, all the defaults should be driven from this table */
/*
#pgm	q_seq	l_seq	p_seq	matrix	g_open	g_ext	fr_shft	e_cut	ktup
#	-n/-p		-s	-e	-f	-h/-j	-E	argv[3]
fasta	prot(0)	prot(0)	prot(0)	bl50	-10	-2	-	10.0	2
fasta	dna(1)	dna(1)	dna(1)	+5/-4	-14	-4	-	2.0	6
ssearch	prot(0)	prot(0)	prot(0)	bl50	-10	-2	-	10.0	-
ssearch	dna(1)	dna(1)	dna(1)	+5/-4	-14	-4	-	2.0	-
fastx	dna(1)	prot(0)	prot(0)	BL50	-12	-2	-20	5.0	2
fasty	dna(1)	prot(0)	prot(0)	BL50	-12	-2	-20/-24	5.0	2
tfastx	dna(1)	prot(0)	prot(0)	BL50	-14	-2	-20	5.0	2
tfasty	dna(1)	prot(0)	prot(0)	BL50	-14	-2	-20/-24	5.0	2
fasts	prot(0)	prot(0)	prot(0)	MD20-MS	-	-	-	5.0	-
fasts	dna(1)	dna(1)	dna(1)	+2/-4	-	-	-	5.0	1
tfasts	prot(0)	dna(1)	prot(0)	MD10-MS	-	-	-	2.0	1
fastf	prot(0)	prot(0)	prot(0)	MD20	-	-	-	2.0	1
tfastf	prot(0)	dna(1)	prot(0)	MD10	-	-	-	1.0	1
fastm	prot(0)	prot(0)	prot(0)	MD20	-	-	-	5.0	1
fastm	dna(1)	dna(1)	dna(1)	+2/-4	-	-	-	2.0	1
tfastm	prot(0)	dna(1)	prot(0)	MD10	-	-	-	2.0	1
*/

struct pgm_def_str {
  int pgm_id;
  char *prog_func;
  char *pgm_abbr;
  char *iprompt0;
  char *ref_str;
  int PgmDID;
  char *smstr;
  int g_open_mod;
  int gshift;
  int hshift;
  int e_cut;
  int ktup;
};

char *ref_str_a[]={
  "\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448\n",
  "\nPlease cite:\n T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197; \n W.R. Pearson (1991) Genomics 11:635-650\n",
 "\nPlease cite:\n Pearson et al, Genomics (1997) 46:24-36\n",
 "\nPlease cite:\n Mackey et al. Mol. Cell. Proteomics  (2002) 1:139-147\n",
  "\nPlease cite:\n W.R. Pearson (1996) Meth. Enzymol. 266:227-258\n"
};

#define FA_PID	1
#define SS_PID	2
#define FX_PID	3
#define FY_PID	4
#define FS_PID	5
#define FF_PID	6
#define FM_PID	7
#define RSS_PID	8
#define RFX_PID	9
#define SSS_PID 10	/* old (slow) non-PG Smith-Waterman */
#define TFA_PID	FA_PID+10
#define TFX_PID	FX_PID+10
#define TFY_PID	FY_PID+10
#define TFS_PID	FS_PID+10
#define TFF_PID	FF_PID+10
#define TFM_PID FM_PID+10

struct pgm_def_str
pgm_def_arr[20] = {
  {0, "", "", "", NULL, 400, "", 0, 0, 0, 1.0, 0 },  /* 0 */
  {FA_PID, "FASTA", "fa",
   "FASTA searches a protein or DNA sequence data bank",
   NULL, 401, "BL50", 0, 0, 0, 10.0, 2}, /* 1 - FASTA */
  {SS_PID, "SSEARCH","gsw","SSEARCH searches a sequence data bank",
   NULL, 404, "BL50", 0, 0, 0, 10.0, 0}, /* 2 - SSEARCH */
  {FX_PID, "FASTX","fx",
   "FASTX compares a DNA sequence to a protein sequence data bank",
   NULL, 405, "BL50", -2, -20, 0, 5.0, 2}, /* 3 - FASTX */
  {FY_PID, "FASTY", "fy",
   "FASTY compares a DNA sequence to a protein sequence data bank",
   NULL, 405, "BL50", -2, -20, -24, 5.0, 2}, /* 4 - FASTY */
  {FS_PID, "FASTS", "fs",
   "FASTS compares linked peptides to a protein data bank",
   NULL, 400, "MD20-MS", 0, 0, 0, 5.0, 1}, /* 5 - FASTS */
  {FF_PID, "FASTF", "ff",
   "FASTF compares mixed peptides to a protein databank",
   NULL, 400, "MD20", 0, 0, 0, 2.0, 1 }, /* 6 - FASTF */
  {FM_PID, "FASTM", "fm",
   "FASTM compares ordered peptides to a protein data bank",
     NULL, 400, "MD20", 0, 0, 0, 5.0, 1 }, /* 7 - FASTM */
  {RSS_PID, "PRSS", "rss",
   "PRSS evaluates statistical signficance using Smith-Waterman",
   NULL, 401, "BL50", 0, 0, 0, 1000.0, 0 }, /* 8 - PRSS */
  {RFX_PID,"PRFX", "rfx",
   "PRFX evaluates statistical signficance using FASTX",
   NULL, 401, "BL50", -2, -20, -24, 1000.0, 2 }, /* 9 - PRFX */
  {SSS_PID, "OSEARCH","ssw","OSEARCH searches a sequence data bank",
   NULL, 404, "BL50", 0, 0, 0, 10.0, 0}, /* 2 - OSEARCH */
  {TFA_PID, "TFASTA", "tfa",
   "TFASTA compares a protein  to a translated DNA data bank",
   NULL, 402, "BL50", -2, 0, 0, 5.0, 2 },
  {0, "", "", "", NULL, 400, "", 0, 0, 0, 1.0, 0 },  /* 0 */
  {TFX_PID, "TFASTX", "tfx",
   "TFASTX compares a protein to a translated DNA data bank",
   NULL, 406, "BL50", -2, -20, 0, 2.0, 2},
  {TFY_PID, "TFASTY", "tfy",
   "TFASTY compares a protein to a translated DNA data bank",
   NULL, 406, "BL50", -2, -20, -24, 2.0, 2},
  {TFS_PID, "TFASTS", "tfs",
   "TFASTS compares linked peptides to a translated DNA data bank",
   NULL, 400, "MD10-MS", 0, 0, 0, 2.0, 2 },
  {TFF_PID, "TFASTF", "tff",
   "TFASTF compares mixed peptides to a protein databank",
   NULL, 400, "MD10", 0, 0, 0, 1.0, 1 },
  {TFM_PID, "TFASTM", "tfm",
   "TFASTM compares ordered peptides to a translated DNA databank",
   NULL, 400, "MD10", 0, 0, 0, 1.0, 1 }
};

struct msg_def_str {
  int pgm_id;
  int q_seqt;
  int l_seqt;
  int p_seqt;
  int sw_flag;
  int stages;
  int qframe;
  int nframe;
  int nrelv, srelv, arelv;
  char *f_id0, *f_id1, *label;
};

/* pgm_id    q_seqt     l_seqt   p_seqt sw_f st qf nf nrv srv arv s_ix */
struct msg_def_str msg_def_arr[20] = {
  {0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, "", "", ""},	/* ID=0 */
  {FA_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 1, 3,
   "fa","sw", "opt"},
  {SS_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "sw","sw", "s-w"},
  {FX_PID, SEQT_DNA, SEQT_PROT, SEQT_PROT, 1, 1, 2, -1, 3, 1, 3,
   "fx","sx", "opt"},
  {FY_PID, SEQT_DNA, SEQT_PROT, SEQT_PROT, 1, 1, 2, -1, 3, 1, 3,
   "fy","sy", "opt"},
  {FS_PID, SEQT_UNK, SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "fs","fs", "initn init1"},
  {FF_PID, SEQT_PROT,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "ff","ff", "initn init1"},
  {FM_PID, SEQT_PROT,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 3, 2, 3,
   "fm","fm","initn init1"},
  {RSS_PID, SEQT_UNK,SEQT_PROT, SEQT_PROT, 0, 1, 1, -1, 1, 1, 1,
   "rss","sw","s-w"},
  {RFX_PID, SEQT_DNA,SEQT_PROT, SEQT_PROT, 0, 1, 2, -1, 3, 1, 3,
   "rfx","sx","opt"},
  {SSS_PID, SEQT_UNK,SEQT_PROT, SEQT_PROT, 1, 1, 1, -1, 1, 1, 1,
   "sw","sw", "s-w"},
  {TFA_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 0, 1, 1, 6, 3, 1, 3,
   "tfa","fa","initn init1"},
  {0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, "", "", ""},	/* ID=12 */
  {TFX_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 2, 3, 2, 3,
   "tfx","sx","initn opt"},
  {TFY_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 2, 3, 2, 3,
   "tfy","sy","initn opt"},
  {TFS_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tfs","fs","initn init1"},
  {TFF_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tff","ff","initn init1"},
  {TFM_PID, SEQT_PROT,SEQT_DNA, SEQT_PROT, 1, 1, 1, 6, 3, 2, 3,
   "tfm","fm","initn init1"}
};

int
get_pgm_id() {

  int rval=0;

#ifdef FASTA
#ifndef TFAST
  pgm_def_arr[FA_PID].ref_str = ref_str_a[0];
  rval=FA_PID;
#else
  pgm_def_arr[TFA_PID].ref_str = ref_str_a[0];
  rval=TFA_PID;
#endif
#endif

#ifdef FASTX
#ifndef TFAST
#ifndef PRSS
  pgm_def_arr[FX_PID].ref_str = ref_str_a[2];
  rval=FX_PID;
#else
  pgm_def_arr[RFX_PID].ref_str = ref_str_a[2];
  rval=RFX_PID;
#endif
#else
  pgm_def_arr[TFX_PID].ref_str = ref_str_a[2];
  rval=TFX_PID;
#endif
#endif

#ifdef FASTY
#ifndef TFAST
  pgm_def_arr[FY_PID].ref_str = ref_str_a[2];
  rval=FY_PID;
#else
  pgm_def_arr[TFY_PID].ref_str = ref_str_a[2];
  rval=TFY_PID;
#endif
#endif

#ifdef FASTS
#ifndef TFAST
  pgm_def_arr[FS_PID].ref_str = ref_str_a[3];
  rval=FS_PID;
#else
  pgm_def_arr[TFS_PID].ref_str = ref_str_a[3];
  rval=TFS_PID;
#endif
#endif

#ifdef FASTF
#ifndef TFAST
  pgm_def_arr[FF_PID].ref_str = ref_str_a[3];
  rval=FF_PID;
#else
  pgm_def_arr[TFF_PID].ref_str = ref_str_a[3];
  rval=TFF_PID;
#endif
#endif

#ifdef FASTM
#ifndef TFAST
  pgm_def_arr[FM_PID].ref_str = ref_str_a[3];
  rval=FM_PID;
#else
  pgm_def_arr[TFM_PID].ref_str = ref_str_a[3];
  rval=TFM_PID;
#endif
#endif

#ifdef SSEARCH
  pgm_def_arr[SS_PID].ref_str = ref_str_a[1];
  rval=SS_PID;
#endif

#ifdef OSEARCH
  pgm_def_arr[SSS_PID].ref_str = ref_str_a[1];
  rval=SSS_PID;
#endif

#ifdef PRSS
#ifndef FASTX
  pgm_def_arr[RSS_PID].ref_str = ref_str_a[4];
  rval=RSS_PID;
#endif
#endif

  return rval;
}

char *iprompt1=" test sequence file name: ";
char *iprompt2=" database file name: ";

char *verstr="version 34.26.5 April 26, 2007";

char   *s_optstr = "13Ac:f:g:h:j:k:nopP:r:s:St:Ux:y:";

static int mktup=2;
static int ktup_set = 0;
static int gap_set=0;
static int del_set=0;
static int mshuff_set = 0;
static int prot2dna = 0;

extern int max_workers;

extern void s_abort(char *, char *);
extern void init_ascii(int ext_sq, int *sascii, int dnaseq);
extern int standard_pam(char *smstr, struct pstruct *ppst,
			int del_set, int gap_set);
extern void mk_n_pam(int *arr,int siz, int mat, int mis);
extern int karlin(int , int, double *, double *, double *);
extern void init_karlin_a(struct pstruct *, double *, double **);
extern int do_karlin_a(int **, struct pstruct *, double *,
		       double *, double *, double *, double *);

#if defined(TFAST) || defined(FASTX) || defined(FASTY)
extern void aainit(int tr_type, int debug);
#endif

char *iprompt0, *prog_func, *refstr;


/* Sets defaults assuming a protein sequence */
void h_init (struct pstruct *ppst, struct mngmsg *m_msp, char *pgm_abbr)
{
  struct pgm_def_str pgm_def;
  int i, pgm_id;

  ppst->pgm_id  = pgm_id =  get_pgm_id();
  pgm_def = pgm_def_arr[pgm_id];

  /* check that pgm_def_arr[] is valid */
  if (pgm_def.pgm_id != pgm_id) {
    fprintf(stderr,
	    "**pgm_def integrity failure: def.pgm_id %d != pgm_id %d**\n",
	    pgm_def.pgm_id, pgm_id);
    exit(1);
  }

  /* check that msg_def_arr[] is valid */
  if (msg_def_arr[pgm_id].pgm_id != pgm_id) {
    fprintf(stderr,
	    "**msg_def integrity failure: def.pgm_id %d != pgm_id %d**\n",
	    msg_def_arr[pgm_id].pgm_id, pgm_id);
    exit(1);
  }

  strncpy(pgm_abbr,pgm_def.pgm_abbr,MAX_SSTR);
  iprompt0 = pgm_def.iprompt0;
  refstr = pgm_def.ref_str;
  prog_func = pgm_def.prog_func;

  /* MAXTOT = MAXTST + MAXLIB for everything except TFAST,
     where it is MAXTST + MAXTRN */
  m_msp->max_tot = MAXTOT;

  /* set up DNA query sequence if required*/
  if (msg_def_arr[pgm_id].q_seqt == SEQT_DNA) {
    memcpy(qascii,nascii,sizeof(qascii));
    m_msp->qdnaseq = SEQT_DNA;
  }
  else { 	/* when SEQT_UNK, start with protein */
    memcpy(qascii,aascii,sizeof(qascii));
    m_msp->qdnaseq = msg_def_arr[pgm_id].q_seqt;
  }

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  qascii[','] = ESS;
  /* also initialize aascii, nascii for databases */
  qascii['*'] = NA;
#endif

  /* initialize a pam matrix */
  strncpy(ppst->pamfile,pgm_def.smstr,MAX_FN);
  standard_pam(ppst->pamfile,ppst,del_set,gap_set);
  ppst->have_pam2 = 0;

  /* this is always protein by default */
  ppst->nsq = naa;
  ppst->nsqx = naax;
  for (i=0; i<=ppst->nsqx; i++) {
    ppst->sq[i] = aa[i];
    ppst->hsq[i] = haa[i];
    ppst->sqx[i]=aax[i];	/* sq = aa */
    ppst->hsqx[i]=haax[i];	/* hsq = haa */
  }
  ppst->sq[ppst->nsqx+1] = ppst->sqx[ppst->nsqx+1] = '\0';

  /* set up the c_nt[] mapping */

#if defined(FASTS) || defined(FASTF) || defined(FASTM)
  ppst->c_nt[ESS] = ESS;
#endif
  ppst->c_nt[0]=0;
  for (i=1; i<=nnt; i++) {
    ppst->c_nt[i]=gc_nt[i];
    ppst->c_nt[i+nnt]=gc_nt[i]+nnt;
  }
}

/*
 * alloc_pam(): allocates memory for the 2D pam matrix as well
 * as for the integer array used to transmit the pam matrix
 */
void
alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;
  char err_str[128];

  if ((ppst->pam2[0] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((ppst->pam2[1] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((d2p = pam12 = (int *) calloc (d1 * d2, sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
   }

   for (i = 0; i < d1; i++, d2p += d2)
      ppst->pam2[0][i] = d2p;

   if ((d2p=pam12x= (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2d pam matrix: %d",d2);
     s_abort (err_str,"");
   }

   for (i = 0;  i < d1; i++, d2p += d2)
      ppst->pam2[1][i] = d2p;

   ppst->have_pam2 = 1;
}

/*
 *  init_pam2(struct pstruct pst): Converts 1-D pam matrix to 2-D
 */
void
init_pam2 (struct pstruct *ppst) {
  int     i, j, k, nsq;

  nsq = ppst->nsq;

  ppst->pam2[0][0][0] = -BIGNUM;
  ppst->pam_h = -1; ppst->pam_l = 1;

  k = 0;
  for (i = 1; i <= nsq; i++) {
    ppst->pam2[0][0][i] = ppst->pam2[0][i][0] = -BIGNUM;
    for (j = 1; j <= i; j++) {
      ppst->pam2[0][j][i] = ppst->pam2[0][i][j] = pam[k++] - ppst->pamoff;
      if (ppst->pam_l > ppst->pam2[0][i][j]) ppst->pam_l =ppst->pam2[0][i][j];
      if (ppst->pam_h < ppst->pam2[0][i][j]) ppst->pam_h =ppst->pam2[0][i][j];
    }
  }
}

void
init_pamx (struct pstruct *ppst) {
  int     i, j, k, nsq, pam_xx, pam_xm;
  int sa_x, sa_t, tmp;

  nsq = ppst->nsq;

  ppst->nt_align = (ppst->dnaseq== SEQT_DNA || ppst->dnaseq == SEQT_RNA);

  if (ppst->nt_align) {
    sa_x = pascii['N'];
    sa_t = sa_x;
  }
  else {
    sa_x = pascii['X'];
    sa_t = pascii['*'];
  }

  if (ppst->dnaseq == SEQT_RNA) {
    tmp = ppst->pam2[0][nascii['G']][nascii['G']] - 1;
    ppst->pam2[0][nascii['A']][nascii['G']] = 
      ppst->pam2[0][nascii['C']][nascii['T']] = 
      ppst->pam2[0][nascii['C']][nascii['U']] = tmp;
  }

  if (ppst->pam_x_set) {
    for (i=1; i<=nsq; i++) {
      ppst->pam2[0][sa_x][i] = ppst->pam2[0][i][sa_x]=ppst->pam_xm;
      ppst->pam2[0][sa_t][i] = ppst->pam2[0][i][sa_t]=ppst->pam_xm;
    }
    ppst->pam2[0][sa_x][sa_x]=ppst->pam_xx;
    ppst->pam2[0][sa_t][sa_t]=ppst->pam_xx;
  }
  else {
    ppst->pam_xx = ppst->pam2[0][sa_x][sa_x];
    ppst->pam_xm = ppst->pam2[0][1][sa_x];
  }

  pam_xx = ppst->pam_xx;
  pam_xm = ppst->pam_xm;

  if (ppst->ext_sq_set) {	/* using extended alphabet */
    /* fill in pam2[1] matrix */
    ppst->pam2[1][0][0] = -BIGNUM;
    /* fill in additional parts of the matrix */
    for (i = 1; i <= nsq; i++) {

      /* -BIGNUM to all matches vs 0 */
      ppst->pam2[0][0][i+nsq] = ppst->pam2[0][i+nsq][0] = 
	ppst->pam2[1][0][i+nsq] = ppst->pam2[1][i+nsq][0] = 
	ppst->pam2[1][0][i] = ppst->pam2[1][i][0] = -BIGNUM;

      for (j = 1; j <= nsq; j++) {

	/* replicate pam2[0] to i+nsq, j+nsq */
	ppst->pam2[0][i+nsq][j] = ppst->pam2[0][i][j+nsq] =
	  ppst->pam2[0][i+nsq][j+nsq] = ppst->pam2[1][i][j] =
	  ppst->pam2[0][i][j];

	/* set the high portion of pam2[1] to the corresponding value
	   of pam2[1][sa_x][j] */

	ppst->pam2[1][i+nsq][j] = ppst->pam2[1][i][j+nsq]=
	  ppst->pam2[1][i+nsq][j+nsq]=ppst->pam2[0][sa_x][j];
      }
    }
  }
}

/*  function specific initializations */
void
f_initenv (struct mngmsg *m_msp, struct pstruct *ppst, unsigned char **aa0) {
  struct msg_def_str m_msg_def;
  int pgm_id;

  pgm_id = ppst->pgm_id;
  m_msg_def = msg_def_arr[pgm_id];

  m_msp->last_calc_flg=0;

  strncpy(m_msp->f_id0,m_msg_def.f_id0,sizeof(m_msp->f_id0));
  strncpy(m_msp->f_id1,m_msg_def.f_id1,sizeof(m_msp->f_id1));
  strncpy (m_msp->label, m_msg_def.label, sizeof(m_msp->label));

#ifndef SSEARCH
  strncpy (m_msp->alab[0],"initn",20);
  strncpy (m_msp->alab[1],"init1",20);
  strncpy (m_msp->alab[2],"opt",20);
#else
  strncpy (m_msp->alab[0],"s-w opt",20);
#endif

  ppst->gdelval += pgm_def_arr[pgm_id].g_open_mod;
  ppst->sw_flag = m_msg_def.sw_flag;
  m_msp->e_cut=pgm_def_arr[pgm_id].e_cut;

  ppst->score_ix = 0;
  ppst->histint = 2;
  m_msp->qframe = m_msg_def.qframe;
  ppst->sw_flag = m_msg_def.sw_flag;
  m_msp->nframe = m_msg_def.nframe;
  m_msp->nrelv = m_msg_def.nrelv;
  m_msp->srelv = m_msg_def.srelv;
  m_msp->arelv = m_msg_def.arelv;
  m_msp->stages = m_msg_def.stages;
#if defined(PRSS)
  m_msp->shuff_wid = 0;
  m_msp->shuff_max = 200;
#endif

  /* see param.h for the definition of all these */

  m_msp->qshuffle = 0;
  m_msp->nm0 = 1;
  m_msp->escore_flg = 0;

  /* pam information */
  ppst->pam_pssm = 0;
#if defined(FASTS) || defined(FASTF) || defined(FASTM)
   ppst->pam_xx = ppst->pam_xm = 0;
#else
  ppst->pam_xx = 1;  /* set >0 to use pam['X']['X'] value */
  ppst->pam_xm = -1;  /* set >0 to use pam['X']['A-Z'] value */
#endif
  ppst->pam_x_set = 0;
  ppst->pam_set = 0;
  ppst->pam_pssm = 0;
  ppst->p_d_set = 0;
  ppst->pamoff = 0;
  ppst->ext_sq_set = 0;

  if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
    mktup = 2;
    ppst->param_u.fa.bestscale = 300;
    ppst->param_u.fa.bestoff = 36;
    ppst->param_u.fa.bkfact = 6;
    ppst->param_u.fa.scfact = 3;
    ppst->param_u.fa.bktup = 2;
    ppst->param_u.fa.ktup = 0;
    ppst->param_u.fa.bestmax = 50;
    ppst->param_u.fa.pamfact = 1;
    ppst->param_u.fa.altflag = 0;
    ppst->param_u.fa.optflag = 1;
    ppst->param_u.fa.iniflag = 0;
    ppst->param_u.fa.optcut = 0;
    ppst->param_u.fa.optcut_set = 0;
    ppst->param_u.fa.cgap = 0;
    ppst->param_u.fa.optwid = MAXWINDOW;
  }

}

/*  switches for fasta only */

static int shift_set=0;
static int subs_set=0;
static int sw_flag_set=0;
static int nframe_set=0;
static int wid_set=0;

void
f_getopt (char copt, char *optarg,
	  struct mngmsg *m_msg, struct pstruct *ppst)
{
  int pgm_id;
  char *bp;

  pgm_id = ppst->pgm_id;

  switch (copt) {
  case '1': 
    if (pgm_def_arr[pgm_id].ktup > 0) {
      ppst->param_u.fa.iniflag=1;
    }
     break;
  case '3':
    nframe_set = 1;
    if (pgm_id == TFA_PID) {
      m_msg->nframe = 3; break;
    }
    else {
      m_msg->nframe = 1;	/* for TFASTXY */
      m_msg->qframe = 1;  /* for FASTA, FASTX */
    }
    break;
  case 'A':
    ppst->sw_flag= 1;
    sw_flag_set = 1;
    break;
  case 'c':
    if (pgm_def_arr[pgm_id].ktup > 0) {
      sscanf (optarg, "%d", &ppst->param_u.fa.optcut);
      ppst->param_u.fa.optcut_set = 1;
    }
    break;
  case 'f':
    sscanf (optarg, "%d", &ppst->gdelval);
    if (ppst->gdelval > 0) ppst->gdelval = -ppst->gdelval;
    del_set = 1;
    break;
  case 'g':
    sscanf (optarg, "%d", &ppst->ggapval);
    if (ppst->ggapval > 0) ppst->ggapval = -ppst->ggapval;
    gap_set = 1;
    break;
  case 'h':
    sscanf (optarg, "%d", &ppst->gshift);
    if (ppst->gshift > 0) ppst->gshift = -ppst->gshift;
    shift_set = 1;
    break;
  case 'j':
    sscanf (optarg, "%d", &ppst->gsubs);
    subs_set = 1;
    break;
  case 'k':
    sscanf (optarg, "%d", &m_msg->shuff_max);
    mshuff_set = 1;
    break;
  case 'n':
    m_msg->qdnaseq = SEQT_DNA;
    re_ascii(qascii,nascii);
    strncpy(m_msg->sqnam,"nt",4);
    prot2dna = 1;
    break;
  case 'o':
    if (pgm_def_arr[pgm_id].ktup > 0) {
      ppst->param_u.fa.optflag = 0;
      msg_def_arr[pgm_id].nrelv = m_msg->nrelv = 2;
    }
    break;
  case 'p':
    m_msg->qdnaseq = SEQT_PROT;
    ppst->dnaseq = SEQT_PROT;
    strncpy(m_msg->sqnam,"aa",4);
    break;
  case 'P':
    strncpy(ppst->pgpfile,optarg,MAX_FN);
    if ((bp=strchr(ppst->pgpfile,' '))!=NULL) {
      *bp='\0';
      ppst->pgpfile_type = atoi(bp+1);
    }
    else ppst->pgpfile_type = 0;
    ppst->pgpfile[MAX_FN-1]='\0';
    ppst->pam_pssm = 1;
    break;
  case 'r':
    sscanf(optarg,"%d/%d",&ppst->p_d_mat,&ppst->p_d_mis);
    if (ppst->p_d_mat > 0 && ppst->p_d_mis < 0) {
      ppst->p_d_set = 1;
      strncpy(ppst->pamfile,optarg,40);
    }
    break;
  case 's':
    strncpy (ppst->pamfile, optarg, 120);
    ppst->pamfile[120-1]='\0';
    if (!standard_pam(ppst->pamfile,ppst,del_set, gap_set)) {
      initpam (ppst->pamfile, ppst);
    }
    ppst->pam_set=1;
    break;
  case 'S':	/* turn on extended alphabet for seg */
    ppst->ext_sq_set = 1;
    break;
  case 't':
    if (tolower(optarg[0])=='t') {
      m_msg->term_code = aascii['*']; optarg++;
    }
    if (*optarg) {sscanf (optarg, "%d", &ppst->tr_type);}
    break;
  case 'U':
    m_msg->qdnaseq = SEQT_RNA;
    memcpy(qascii,nascii,sizeof(qascii));
    strncpy(m_msg->sqnam,"nt",4);
    nt[nascii['T']]='U';
    prot2dna=1;
    break;
  case 'x':
    if (strchr(optarg,',')!=NULL) {
      sscanf (optarg,"%d,%d",&ppst->pam_xx, &ppst->pam_xm);
    }
    else {
      sscanf (optarg,"%d",&ppst->pam_xx);
      ppst->pam_xm = ppst->pam_xx;
    }
    ppst->pam_x_set=1;
    break;
  case 'y':
    if (pgm_def_arr[pgm_id].ktup > 0) {
      sscanf (optarg, "%d", &ppst->param_u.fa.optwid);
      wid_set = 1;
    }
    break;
  }
}

void
f_lastenv (struct mngmsg *m_msg, struct pstruct *ppst)
{
  char save_str[MAX_SSTR];

#if !defined(FASTM) && !defined(FASTS) && !defined(FASTF)
  strncpy(save_str,"*",sizeof(save_str));
#else
  strncpy(save_str,",",sizeof(save_str));
#endif

  if (m_msg->qdnaseq == SEQT_UNK) {
    build_xascii(qascii,save_str);
    if (m_msg->ann_flg) ann_ascii(qascii,m_msg->ann_arr);
  }  

/* this check allows lc DNA sequence queries with FASTX */
#if defined(FASTA) && !defined(FASTS) && !defined(FASTM) && !defined(FASTF)
  else
   init_ascii(ppst->ext_sq_set,qascii,m_msg->qdnaseq);
#endif
}

void
f_getarg (int argc, char **argv, int optind,
	  struct mngmsg *m_msg, struct pstruct *ppst)
{

  if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
    if (argc - optind >= 4) {
      sscanf (argv[optind + 3], "%d", &ppst->param_u.fa.ktup);
      ktup_set = 1;
    }
    else
      ppst->param_u.fa.ktup = -ppst->param_u.fa.bktup;
  }
  
  if (ppst->pgm_id == RSS_PID && argc - optind > 3) {
    sscanf (argv[optind + 3], "%d", &m_msg->shuff_max);
  }

  if (ppst->pgm_id == RFX_PID && argc - optind > 4) {
    sscanf (argv[optind + 4], "%d", &m_msg->shuff_max);
  }
}

/* fills in the query ascii mapping from the parameter
   ascii mapping.
*/

void
re_ascii(int *qascii, int *pascii) {
  int i;

  for (i=0; i < 128; i++) {
    if (qascii[i] > '@' || qascii[i] < ESS) {
      qascii[i] = pascii[i];
    }
  }
}


/* recode has become function specific to accommodate FASTS/M */
/* modified 28-Dec-2004 to ensure that all mapped characters
   are valid */
int
recode(unsigned char *seq, int n, int *qascii, int nsqx) {
  int i,j;
  char save_c;

#if defined(FASTS) || defined(FASTM)
  qascii[',']=ESS;
#endif

  for (i=0; i < n; i++) {
    save_c = seq[i];
    if (seq[i] > '@') seq[i] = qascii[seq[i]];
    if (seq[i] > nsqx && seq[i]!=ESS) {
      fprintf(stderr, "*** Warning - unrecognized residue at %d:%c - %2d\n",
	      i,save_c,save_c);
      seq[i] = qascii['X'];
    }
  }
  seq[i]=EOSEQ;
  return i;
}

/* here we have the query sequence, all the command line options,
   but we need to set various parameter options based on the type
   of the query sequence (m_msg->qdnaseq = 0:protein/1:DNA) and
   the function (FASTA/FASTX/TFASTA)
*/

/* this resetp is for conventional a FASTA/TFASTXYZ search */
void
resetp (struct mngmsg *m_msg, struct pstruct *ppst) {
  int i, pgm_id;

  pgm_id = ppst->pgm_id;

#if defined(TFAST)
  if (m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA) {
    fprintf(stderr," %s compares a protein to a translated\n\
DNA sequence library.  Do not use a DNA query/scoring matrix.\n",prog_func);
    exit(1);
  }
#else
#if (defined(FASTX) || defined(FASTY))
  if (!(m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA)) {
    fprintf(stderr," FASTX/Y compares a DNA sequence to a protein database\n");
    fprintf(stderr," Use a DNA query\n");
    exit(1);
  }
#endif
#endif

/* this code changes parameters for programs (FA_PID, SS_PID, FS_PID,
   RSS_PID) that can examine either protein (initial state) or DNA 
   Modified May, 2006 to reset e_cut for DNA comparisons.
*/

  if (msg_def_arr[pgm_id].q_seqt == SEQT_UNK) {
    if (m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA) {
      msg_def_arr[pgm_id].q_seqt = m_msg->qdnaseq;
      msg_def_arr[pgm_id].p_seqt = SEQT_DNA;
      msg_def_arr[pgm_id].l_seqt = SEQT_DNA;
      if (m_msg->qdnaseq == SEQT_DNA) msg_def_arr[pgm_id].qframe = 2;
      pgm_def_arr[pgm_id].e_cut /= 5.0;
    }
    else {
      msg_def_arr[pgm_id].q_seqt = SEQT_PROT;
    }
  }

  ppst->dnaseq = msg_def_arr[pgm_id].p_seqt;
  if (!sw_flag_set) ppst->sw_flag = msg_def_arr[pgm_id].sw_flag;
  if (!m_msg->e_cut_set) m_msg->e_cut=pgm_def_arr[pgm_id].e_cut;

  if (ppst->dnaseq == SEQT_DNA && m_msg->qdnaseq==SEQT_RNA) {
    ppst->dnaseq = SEQT_RNA;
    ppst->nt_align = 1;
  }
  if (ppst->dnaseq==SEQT_DNA) pascii = &nascii[0];
  else if (ppst->dnaseq==SEQT_RNA) {
    pascii = &nascii[0];
    ppst->sq[nascii['T']] = 'U';
  }
  else pascii = &aascii[0];
  m_msg->ldnaseq = msg_def_arr[pgm_id].l_seqt;
  if (m_msg->ldnaseq & SEQT_DNA) {
    memcpy(lascii,nascii,sizeof(lascii));
#ifndef TFAST
#ifdef DNALIB_LC
   init_ascii(ppst->ext_sq_set,lascii,m_msg->ldnaseq);
#endif
#else
  /* no init_ascii() because we translate lower case library sequences */
#endif
  }
  else {
    memcpy(lascii,aascii,sizeof(lascii));	/* initialize lib mapping */

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
    lascii['*'] = NA;
#endif
    init_ascii(ppst->ext_sq_set,lascii,m_msg->ldnaseq);
  }

  if (!nframe_set) {
    m_msg->qframe = msg_def_arr[pgm_id].qframe;
    m_msg->nframe = msg_def_arr[pgm_id].nframe;
  }

  /* the possibilities:
  	     -i  -3	qframe	revcomp
   FA_D/FX   -    -        2       0	
   FA_D/FX   +    -        2       1	
   FA_D/FX   -    +        1       0	
   FA_D/FX   +    +        2       1	
  */

  if (m_msg->qdnaseq == SEQT_DNA) {
    m_msg->nframe = 1;
    if (m_msg->qframe == 1 && m_msg->revcomp==1) {
      m_msg->qframe = m_msg->revcomp+1;
    }
  }
  else if (m_msg->qdnaseq == SEQT_RNA) {
    m_msg->qframe = m_msg->revcomp+1;
    m_msg->nframe = 1;
  }

  /* change settings for DNA search */
  if (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA) {
    ppst->histint = 4;

    if (!del_set) {
#ifdef OLD_FASTA_GAP
      ppst->gdelval = -16;	/* def. del penalty */
#else
      ppst->gdelval = -12;	/* def. open penalty */
#endif
    }
    if (!gap_set) ppst->ggapval = -4;	/* def. gap penalty */

    if (pgm_def_arr[pgm_id].ktup > 0) {
      /* these parameters are used to scale optcut, they should be replaced
	 by statistically based parameters */
      if (!wid_set) ppst->param_u.fa.optwid = 16;
      ppst->param_u.fa.bestscale = 80;
      ppst->param_u.fa.bkfact = 5;
      ppst->param_u.fa.scfact = 1;
      ppst->param_u.fa.bktup = 6;
      ppst->param_u.fa.bestmax = 80;
      ppst->param_u.fa.bestoff = 45;

      if (!sw_flag_set) {
	ppst->sw_flag = 0;
	strncpy(m_msg->f_id1,"bs",sizeof(m_msg->f_id1));
      }

      /* largest ktup */
      mktup = 6;
      
      if (ppst->param_u.fa.pamfact >= 0) ppst->param_u.fa.pamfact = 0;
      if (ppst->param_u.fa.ktup < 0)
	ppst->param_u.fa.ktup = -ppst->param_u.fa.bktup;
    }

    ppst->nsq = nnt;
    ppst->nsqx = nntx;
    for (i=0; i<=ppst->nsqx; i++) {
      ppst->hsq[i] = hnt[i];
      ppst->sq[i] = nt[i];
      ppst->hsqx[i] = hntx[i];
      ppst->sqx[i] = ntx[i];
    }
    ppst->sq[ppst->nsqx+1] = ppst->sqx[ppst->nsqx+1] = '\0';

    if (!ppst->pam_set) {
      if (ppst->p_d_set)
	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
#if !defined(FASTS) && !defined(FASTM)
      else if (ppst->pamfile[0]=='\0' || strncmp(ppst->pamfile,"BL50",4)==0) {
	strncpy (ppst->pamfile, "+5/-4", sizeof(ppst->pamfile));
      }
#else
      else if (strncmp(ppst->pamfile,"MD20",4)==0) {
	strncpy (ppst->pamfile, "+2/-2", sizeof(ppst->pamfile));
	ppst->p_d_mat = +2;
	ppst->p_d_mis = -2;
      	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
      }
#endif
      pam = npam;
    }

    strncpy (m_msg->sqnam, "nt",sizeof(m_msg->sqnam));
    strncpy (m_msg->sqtype, "DNA",sizeof(m_msg->sqtype));
  }	/* end DNA reset */

  else {  /* other parameters for protein comparison */
    if (pgm_def_arr[pgm_id].ktup > 0) {
      if (!wid_set) {
	if (ppst->param_u.fa.ktup==1) ppst->param_u.fa.optwid = 32;
	else ppst->param_u.fa.optwid = 16;
      }
    }
    if (!del_set) {ppst->gdelval += pgm_def_arr[pgm_id].g_open_mod;}
    if (!shift_set) {ppst->gshift = pgm_def_arr[pgm_id].gshift;}
    if (!subs_set) {ppst->gsubs = pgm_def_arr[pgm_id].hshift;}
  }

}

/* query_parm()	this function asks for any additional parameters
	that have not been provided.  Could be null. */
void
query_parm (struct mngmsg *m_msp, struct pstruct *ppst)
{
   char    qline[40];

   if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
     if (ppst->param_u.fa.ktup < 0)
       ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

     if (ppst->param_u.fa.ktup == 0) {
       printf (" ktup? (1 to %d) [%d] ", mktup, ppst->param_u.fa.bktup);
       if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
       else sscanf(qline,"%d",&ppst->param_u.fa.ktup);
     }
     if (ppst->param_u.fa.ktup == 0)
       ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
     else ktup_set = 1;
   }

#if defined(PRSS)
   if (m_msp->shuff_max < 10) m_msp->shuff_max = 200;

   if (!mshuff_set) {
     printf(" number of shuffles [%d]? ",m_msp->shuff_max);
     fflush(stdout);
     if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
     else sscanf(qline,"%d",&m_msp->shuff_max);
   }

   if (ppst->zs_win == 0) {
     printf (" local (window) (w) or uniform (u) shuffle [u]? ");
     if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
     else if (qline[0]=='w' || qline[0]=='W') {
       m_msp->shuff_wid = 20;
       printf(" local shuffle window size [%d]? ",m_msp->shuff_wid);
       if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
       else sscanf(qline,"%d",&m_msp->shuff_wid);
     }
   }
#endif
}

/* last_init() cannot look at aa0, n0, because it is only run once,
   it is not run before each new aa0 search */
void
last_init (struct mngmsg *m_msg, struct pstruct *ppst
#ifdef PCOMPLIB
	   ,int nnodes
#endif
	   )
{
  int ix_l, ix_i, i, pgm_id;
  double *kar_p;
  double aa0_f[MAXSQ];

  pgm_id = ppst->pgm_id;

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  m_msg->nohist = 1;
  m_msg->shuff_max = 2000;
#ifndef PCOMPLIB
  ppst->shuff_node = m_msg->shuff_max/max_workers;
#else
  ppst->shuff_node = m_msg->shuff_max/nnodes;
#endif
#endif

  if (m_msg->aln.llen < 1) {
    m_msg->aln.llen = 60;
  }

#ifndef PCOMPLIB
#if defined(FASTX) || defined(FASTY) || defined(TFAST)
  /* set up translation tables: faatran.c */
  aainit(ppst->tr_type,ppst->debug_lib);
#endif
#endif

/* a sanity check */
#if !defined(TFAST)
   if (m_msg->revcomp && m_msg->qdnaseq!=SEQT_DNA && m_msg->qdnaseq!=SEQT_RNA) {
     fprintf(stderr," cannot reverse complement protein\n");
     m_msg->revcomp = 0;
   }
#endif

   if (pgm_def_arr[pgm_id].ktup > 0) {

     if (ppst->param_u.fa.ktup < 0)
       ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

     if (ppst->param_u.fa.ktup < 1 || ppst->param_u.fa.ktup > mktup) {
       fprintf(stderr," warning ktup = %d out of range [1..%d], reset to %d\n",
	       ppst->param_u.fa.ktup, mktup, ppst->param_u.fa.bktup);
       ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
     }
   }

   if (pgm_id == TFA_PID) {
     m_msg->revcomp *= 3;
     if (m_msg->nframe == 3) m_msg->nframe += m_msg->revcomp;
   }
   else if (pgm_id == TFX_PID || pgm_id == TFY_PID) {
     if (m_msg->nframe == 1) m_msg->nframe += m_msg->revcomp;
   }

#if !defined(TFAST)
  /* for fasta/fastx searches, itt iterates the the query strand */
  m_msg->nitt1 = m_msg->qframe-1;
#else
  /* for tfasta/tfastxy searches, itt iterates the library frames */
  m_msg->nitt1 = m_msg->nframe-1;
#endif

  if (pgm_def_arr[pgm_id].ktup > 0) {
    if (ppst->param_u.fa.ktup>=2 && !wid_set) {
      ppst->param_u.fa.optwid=16;
      switch (pgm_id) {
      case FA_PID:
	m_msg->thr_fact = 32;
	break;
      case FX_PID:
      case FY_PID:
	m_msg->thr_fact = 16;
	break;
      case TFA_PID:
      case TFX_PID:
      case TFY_PID:
	m_msg->thr_fact = 8;
	break;
      default:
	m_msg->thr_fact = 4;
      }
    }
    else { m_msg->thr_fact = 4;}
  }
  else m_msg->thr_fact = 4;

#if defined(PRSS)
   if (m_msg->shuff_max < 10) m_msg->shuff_max = 200;
   if (ppst->zsflag < 10) ppst->zsflag += 10;
   if (ppst->zs_win > 0) {
     m_msg->shuff_wid = ppst->zs_win;
   }
#endif

   if (pgm_def_arr[ppst->pgm_id].ktup > 0) {
     if (ppst->param_u.fa.iniflag) {
       ppst->score_ix = 1;
       strncpy (m_msg->label, "initn init1", sizeof(m_msg->label));
     }
     else if (ppst->param_u.fa.optflag) {
       ppst->score_ix = 2;
       m_msg->stages = 1;
     }
   }

   if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
     init_pam2(ppst);
   }
   init_pamx(ppst);

   if (ppst->pam_ms) {
     if (m_msg->qdnaseq == SEQT_PROT) {
       /* code to make 'L'/'I' identical scores */
       ix_l = pascii['L'];
       ix_i = pascii['I'];
       ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
	 ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
	 (ppst->pam2[0][ix_l][ix_l]+ppst->pam2[0][ix_i][ix_i]+1)/2;
       for (i=1; i<=ppst->nsq; i++) {
	 ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	   (ppst->pam2[0][i][ix_l]+ppst->pam2[0][i][ix_i]+1)/2;
	 ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	   (ppst->pam2[0][ix_i][i]+ppst->pam2[0][ix_l][i]+1)/2;
       }

       /* code to make 'Q'/'K' identical scores */
       if (!shift_set) {
	 ix_l = pascii['Q'];
	 ix_i = pascii['K'];
	 ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
	   ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
	   (ppst->pam2[0][ix_l][ix_l]+ppst->pam2[0][ix_i][ix_i]+1)/2;
	 for (i=1; i<=ppst->nsq; i++) {
	   ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	     (ppst->pam2[0][i][ix_l]+ppst->pam2[0][i][ix_i]+1)/2;
	   ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	     (ppst->pam2[0][ix_i][i]+ppst->pam2[0][ix_l][i]+1)/2;
	 }
       }
     }
   }

   /*
   print_pam(ppst);
   */

   /* once we have a complete pam matrix, we can calculate Lambda and K 
      for "average" sequences */
   kar_p = NULL;
   init_karlin_a(ppst, aa0_f, &kar_p);
   do_karlin_a(ppst->pam2[0], ppst, aa0_f,
	       kar_p, &m_msg->Lambda, &m_msg->K, &m_msg->H);
   free(kar_p);

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
   if (ppst->ext_sq_set) {
     fprintf(stderr," -S not available on [t]fast[fs]\n");
     ppst->ext_sq_set = 0;

     /* reset sascii to ignore -S, map lc */
     init_ascii(0,lascii,0);
   }
#endif
}

/* this function is left over from the older FASTA format scoring
   matrices that allowed additional parameters (bktup, bkfact) to be
   set in the scoring matrix.  It is no longer used.  A modern version
   would set parameters based on lambda and K.
*/
/*
void
f_initpam (line, ppst)
char   *line;
struct pstruct *ppst;
{
   if (sscanf (line, " %d %d %d %d %d %d %d", &ppst->param_u.fa.scfact,
	       &ppst->param_u.fa.bestoff, &ppst->param_u.fa.bestscale,
	       &ppst->param_u.fa.bkfact, &ppst->param_u.fa.bktup,
	       &ppst->param_u.fa.bestmax, &ppst->histint) != 7)
   {
      printf ("  bestcut parameters - bad format\n");
      exit (1);
   }
}
*/

/* alloc_pam2 creates a profile structure */
int **
alloc_pam2p(int len, int nsq) {
  int i;
  int **pam2p;

  if ((pam2p = (int **)calloc(len,sizeof(int *)))==NULL) {
    fprintf(stderr," Cannot allocate pam2p: %d\n",len);
    return NULL;
  }

  if((pam2p[0] = (int *)calloc((nsq+1)*len,sizeof(int)))==NULL) {
    fprintf(stderr, "Cannot allocate pam2p[0]: %d\n", (nsq+1)*len);
    free(pam2p);
    return NULL;
  }

  for (i=1; i<len; i++) {
    pam2p[i] = pam2p[0] + (i*(nsq+1));
  }

  return pam2p;
}

void free_pam2p(int **pam2p) {
  if (pam2p) {
    free(pam2p[0]);
    free(pam2p);
  }
}

/* sortbest has now become comparison function specific so that we can use
   a different comparison for fasts/f 
*/
#if !defined(FASTS) && !defined (FASTF) && !defined(FASTM)
#ifndef PCOMPLIB
void
qshuffle() {}
#endif

int
last_calc(unsigned char *aa0, unsigned char *aa1, int maxn,
	  struct beststr **bestp_arr, int nbest,
	  struct mngmsg *m_msg, struct pstruct *pst,
	  void **f_str, void *rs_str)
{
  return nbest;
}

void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->score[irelv] >= bptr[j + gap]->score[irelv]) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

void show_aux(FILE *fp, struct beststr *bptr) {}
void header_aux(FILE *fp) {}

#else
void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->escore < bptr[j + gap]->escore) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

#if defined(FASTS) || defined(FASTM)

#ifndef PCOMPLIB
/* this shuffle is for FASTS  */
/* convert ',' -> '\0', shuffle each of the substrings */
void
qshuffle(unsigned char *aa0, int n0, int nm0) {

  unsigned char **aa0start, *aap, tmp;
  int i,j,k, ns;

  if ((aa0start=(unsigned char **)calloc(nm0+1,
					 sizeof(unsigned char *)))==NULL) {
    fprintf(stderr,"cannot calloc for qshuffle %d\n",nm0);
    exit(1);
  }

  aa0start[0]=aa0;
  for (k=1,i=0; i<n0; i++) {
    if (aa0[i]==EOSEQ || aa0[i]==ESS) {
      aa0[i]='\0';
      aa0start[k++] = &aa0[i+1];
    }
  }  

  /* aa0start has the beginning of each substring */
  for (k=0; k<nm0; k++) {
    aap=aa0start[k];
    ns = strlen((char *)aap);
    for (i=ns; i>1; i--) {
      j = nrand(i);
      tmp = aap[j];
      aap[j] = aap[i-1];
      aap[i-1] = tmp;
    }
    aap[ns] = 0;
  }

  for (k=1; k<nm0; k++) {
/*  aap = aa0start[k];
    while (*aap) fputc(pst.sq[*aap++],stderr);
    fputc('\n',stderr);
*/
    aa0start[k][-1]=ESS;
  }

  free(aa0start);
}
#endif
#endif

#ifdef FASTF
#ifndef PCOMPLIB
void qshuffle(unsigned char *aa0, int n0, int nm0) {

  int i, j, k, nmpos;
  unsigned char tmp;
  int nmoff;
  
  nmoff = (n0 - nm0 - 1)/nm0 + 1;

  for (i = nmoff-1 ; i > 0 ; i--) {

    /* j = nrand(i); if (i == j) continue;*/       /* shuffle columns */ 
    j = (nmoff -1 ) - i; 
    if (i <= j) break; /* reverse columns */

    /* swap all i'th column residues for all j'th column residues */
    for(nmpos = 0, k = 0 ; k < nm0 ; k++, nmpos += nmoff+1 ) {
      tmp = aa0[nmpos + i];
      aa0[nmpos + i] = aa0[nmpos + j];
      aa0[nmpos + j] = tmp;
    }
  }
}
#endif
#endif


/* show additional best_str values */
void show_aux(FILE *fp, struct beststr *bptr) {
  fprintf(fp," %2d %3d",bptr->segnum,bptr->seglen);
}

void header_aux(FILE *fp) {
  fprintf(fp, " sn  sl");
}
#endif

void
fill_pam(int **pam2p, int n0, int nsq, double **freq2d, double scale) {
  int i, j;
  double freq;

  /* fprintf(stderr, "scale: %g\n", scale); */
  
  /* now fill in the pam matrix: */
  for (i = 0 ; i < n0 ; i++) {
    for (j = 1 ; j <=20 ; j++) {
      freq = scale * freq2d[i][j-1];
      if ( freq < 0.0) freq -= 0.5;
      else freq += 0.5;
      pam2p[i][j] = (int)(freq);
    }
  }
}

double
get_lambda(int **pam2p, int n0, int nsq, unsigned char *query) {
  double lambda, H;
  double *pr, tot, sum;
  int i, ioff, j, min, max;

  /* get min and max scores */
  min = BIGNUM;
  max = -BIGNUM;
  if(pam2p[0][1] == -BIGNUM) {
    ioff = 1;
    n0++;
  } else {
    ioff = 0;
  }

  for (i = ioff ; i < n0 ; i++) {
    for (j = 1; j <= nsq ; j++) {
      if (min > pam2p[i][j])
	min = pam2p[i][j];
      if (max < pam2p[i][j])
	max = pam2p[i][j];
    }
  }

  /*  fprintf(stderr, "min: %d\tmax:%d\n", min, max); */
  
  if ((pr = (double *) calloc(max - min + 1, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for score probabilities: %d\n", max - min + 1);
    exit(1);
  }

  tot = (double) rrtotal * (double) rrtotal * (double) n0;
  for (i = ioff ; i < n0 ; i++) {
    for (j = 1; j <= nsq ; j++) {
      pr[pam2p[i][j] - min] +=
	(double) ((double) rrcounts[aascii[query[i]]] * (double) rrcounts[j]) / tot;
    }
  }

  sum = 0.0;
  for(i = 0 ; i <= max-min ; i++) { 
    sum += pr[i];
    /*     fprintf(stderr, "%3d: %g %g\n", i+min, pr[i], sum); */
  }
  /*   fprintf(stderr, "sum: %g\n", sum); */

  for(i = 0 ; i <= max-min ; i++) { pr[i] /= sum; }

  if (!karlin(min, max, pr, &lambda, &H)) {
    fprintf(stderr, "Karlin lambda estimation failed\n");
  }

  /*   fprintf(stderr, "lambda: %g\n", lambda); */
  free(pr);

  return lambda;
}

/*
   *aa0 - query sequence
   n0   - length
   pamscale - scaling for pam matrix - provided by apam.c, either
              0.346574 = ln(2)/2 (P120, BL62) or
	      0.231049 = ln(2)/3 (P250, BL50) 
*/

void
scale_pssm(int **pssm2p, double **freq2d,
	   unsigned char *query, int n0,
	   int **pam2, double pamscale);

static unsigned char ustandard_aa[] ="\0ARNDCQEGHILKMFPSTWYV";

void
read_pssm(unsigned char *aa0, int n0, int nsq,
	  double pamscale, 
	  FILE *fp, int pgpf_type, struct pstruct *ppst) {
  int i, j, len, k;
  int qi, rj;	/* qi - index query; rj - index residues (1-20) */
  int **pam2p;
  int first, too_high;
  unsigned char *query, ctmp;
  char dline[512];
  double freq, **freq2d, lambda, new_lambda;
  double scale, scale_high, scale_low;

  pam2p = ppst->pam2p[0];

  if (pgpf_type == 0) {

    if(1 != fread(&len, sizeof(int), 1, fp)) {
      fprintf(stderr, "error reading from checkpoint file: %d\n", len);
      exit(1);
    }

    if(len != n0) {
      fprintf(stderr, "profile length (%d) and query length (%d) don't match!\n",
	      len,n0);
      exit(1);
    }

    /* read over query sequence stored in BLAST profile */
    if(NULL == (query = (unsigned char *) calloc(len+2, sizeof(char)))) {
      fprintf(stderr, "Couldn't allocate memory for query!\n");
      exit(1);
    }

    if(len != fread(query, sizeof(char), len, fp)) {
      fprintf(stderr, "Couldn't read query sequence from profile: %s\n", query);
      exit(1);
    }
  }
  else if (pgpf_type == 1) {

    if ((fgets(dline,sizeof(dline),fp) == NULL)  ||
	(1 != sscanf(dline, "%d",&len))) {
      fprintf(stderr, "error reading from checkpoint file: %d\n", len);
      exit(1);
    }

    if(len != n0) {
      fprintf(stderr, "profile length (%d) and query length (%d) don't match!\n",
	      len,n0);
      exit(1);
    }

    /* read over query sequence stored in BLAST profile */
    if(NULL == (query = (unsigned char *) calloc(len+2, sizeof(char)))) {
      fprintf(stderr, "Couldn't allocate memory for query!\n");
      exit(1);
    }

    if (fgets((char *)query,len+2,fp)==NULL) {
      fprintf(stderr, "Couldn't read query sequence from profile: %s\n", query);
      exit(1);
    }
  }  
  else {
    fprintf(stderr," Unrecognized PSSM file type: %d\n",pgpf_type);
    exit(1);
  }

  /* currently we don't do anything with query; ideally, we should
     check to see that it actually matches aa0 ... */

  /* quick 2d array alloc: */
  if((freq2d = (double **) calloc(n0, sizeof(double *))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  if((freq2d[0] = (double *) calloc(n0 * 20, sizeof(double))) == NULL) {
    fprintf(stderr, "Couldn't allocate memory for frequencies!\n");
    exit(1);
  }

  /* a little pointer arithmetic to fill out 2d array: */
  for (i = 1 ; i < n0 ; i++) {
    freq2d[i] = freq2d[i-1] + 20;
  }

  if (pgpf_type == 0) {
    for (qi = 0 ; qi < n0 ; qi++) {
      for (rj = 0 ; rj < 20 ; rj++) {
	if(1 != fread(&freq, sizeof(double), 1, fp)) {
	  fprintf(stderr, "Error while reading frequencies!\n");
	  exit(1);
	}
	freq2d[qi][rj] = freq;
      }
    }
  }
  else {
    for (qi = 0 ; qi < n0 ; qi++) {
      if ((fgets(dline,sizeof(dline),fp) ==NULL) ||
      (k = sscanf(dline,"%c %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		 &ctmp, &freq2d[qi][0], &freq2d[qi][1], &freq2d[qi][2], &freq2d[qi][3], &freq2d[qi][4], 
		 &freq2d[qi][5], &freq2d[qi][6], &freq2d[qi][7], &freq2d[qi][8], &freq2d[qi][9],
		 &freq2d[qi][10], &freq2d[qi][11], &freq2d[qi][12], &freq2d[qi][13], &freq2d[qi][14],
		      &freq2d[qi][15], &freq2d[qi][16], &freq2d[qi][17], &freq2d[qi][18], &freq2d[qi][19]))<1) {
	fprintf(stderr, "Error while reading frequencies: %d read!\n",k);
	exit(1);
      }
      for (rj=0; rj<20; rj++) { freq2d[qi][rj] /= 10.0; }	/* reverse scaling */
    }
  }

  scale_pssm(ppst->pam2p[0], freq2d, query, n0, ppst->pam2[0],pamscale);

  free(freq2d[0]);
  free(freq2d);

  free(query);
}

void
scale_pssm(int **pssm2p, double **freq2d, unsigned char *query, int n0, int **pam2, double pamscale) {
  int i, qi, rj;
  double freq, new_lambda, lambda;
  int first, too_high;
  double scale, scale_high, scale_low;

  for (qi = 0 ; qi < n0 ; qi++) {
    for (rj = 0 ; rj < 20 ; rj++) {
      if (freq2d[qi][rj] > 1e-20) {
	freq = log(freq2d[qi][rj] /((double) (rrcounts[rj+1])/(double) rrtotal));
	freq /= pamscale; /* this gets us close to originial pam scores */
	freq2d[qi][rj] = freq;
      }
      else {		
	/* when blastpgp decides to leave something out, it puts 0's in all the frequencies
	   in the binary checkpoint file.  In the ascii version, however, it uses BLOSUM62
	   values.  I will put in scoring matrix values as well */

	freq2d[qi][rj] = pam2[aascii[query[qi]]][rj+1];
      }
    }
  }

  /* now figure out the right scale */
  scale = 1.0;
  lambda = get_lambda(pam2, 20, 20, ustandard_aa);

  /* should be near 1.0 because of our initial scaling by ppst->pamscale */
  /* fprintf(stderr, "real_lambda: %g\n", lambda); */

  /* get initial high/low scale values: */
  first = 1;
  while (1) {
    fill_pam(pssm2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pssm2p, n0, 20, query); 

    if (new_lambda > lambda) {
      if (first) {
	first = 0;
	scale = scale_high = 1.0 + 0.05;
	scale_low = 1.0;
	too_high = 1;
      } else {
	if (!too_high) break;
	scale = (scale_high += scale_high - 1.0);
      }
    } else if (new_lambda > 0) {
      if (first) {
	first = 0;
	scale_high = 1.0;
	scale = scale_low = 1.0 - 0.05;
	too_high = 0;
      } else {
	if (too_high) break;
	scale = (scale_low += scale_low - 1.0);
      }
    } else {
      fprintf(stderr, "new_lambda (%g) <= 0; matrix has positive average score", new_lambda);
      exit(1);
    }
  }

  /* now do binary search between low and high */
  for (i = 0 ; i < 10 ; i++) {
    scale = 0.5 * (scale_high + scale_low);
    fill_pam(pssm2p, n0, 20, freq2d, scale);
    new_lambda = get_lambda(pssm2p, n0, 20, query);
    
    if (new_lambda > lambda) scale_low = scale;
    else scale_high = scale;
  }

  scale = 0.5 * (scale_high + scale_low);
  fill_pam(pssm2p, n0, 20, freq2d, scale);

  /*
  fprintf(stderr, "final scale: %g\n", scale);

  for (qi = 0 ; qi < n0 ; qi++) {
    fprintf(stderr, "%4d %c:  ", qi+1, query[qi]);
    for (rj = 1 ; rj <= 20 ; rj++) {
      fprintf(stderr, "%4d", pssm2p[qi][rj]);
    }
    fprintf(stderr, "\n");
  }
  */
}

#if defined(SSEARCH) || (defined(PRSS) && !defined(FASTX))
int
parse_pssm_asn_fa(FILE *afd, int *n_rows, int *n_cols,
		  unsigned char **query, double ***freqs,
		  char *matrix, int *gap_open, int *gap_extend,
		  double *lambda);

/* the ASN.1 pssm includes information about the scoring matrix used
   (though not the gap penalty in the current version PSSM:2) The PSSM
   scoring matrix and gap penalties should become the default if they
   have not been set explicitly.
*/

int
read_asn_pssm(unsigned char *aa0, int n0, int nsq,
	      double pamscale, FILE *fp, struct pstruct *ppst) {

  int i, j, len, k;
  int qi, rj;	/* qi - index query; rj - index residues (1-20) */
  int **pam2p;
  int first, too_high;
  unsigned char *query, ctmp;
  char dline[512];
  char matrix[MAX_SSTR];
  double psi2_lambda;
  double freq, **freq2d, lambda, new_lambda;
  double scale, scale_high, scale_low;
  int gap_open, gap_extend;
  int n_rows, n_cols;

  pam2p = ppst->pam2p[0];

  if (parse_pssm_asn_fa(fp, &n_rows, &n_cols, &query, &freq2d,
			matrix, &gap_open, &gap_extend, &psi2_lambda)<=0) {
    return -1;
  }

  if (!gap_set) {
    if (gap_open) {
      if (gap_open > 0) {gap_open = -gap_open;}
      ppst->gdelval = gap_open;
    }
    else if (strncmp(matrix,"BLOSUM62",8)==0) {
      ppst->gdelval = -11;
    }
    gap_set = 1;
  }
  if (!del_set) {
    if (gap_extend) {
      if (gap_extend > 0) {gap_extend = -gap_extend;}
      ppst->ggapval = gap_extend;
    }
    else if (strncmp(matrix,"BLOSUM62",8)==0) {
      ppst->ggapval = -1;
    }
    del_set = 1;
  }

  if (strncmp(matrix, "BLOSUM62", 8)== 0 && !ppst->pam_set) {
    strncpy(ppst->pamfile, "BL62", 120);
    standard_pam(ppst->pamfile,ppst,del_set, gap_set);
    if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
    }
    init_pam2(ppst);
    ppst->pam_set = 1;
  }

  if (n_cols < n0) { 
    fprintf(stderr, " query length: %d != n_cols: %d\n",n0, n_cols);
    exit(1);
  }

  scale_pssm(ppst->pam2p[0], freq2d, query, n0, ppst->pam2[0],pamscale);

  free(freq2d[0]);
  free(freq2d);

  free(query);
  return 1;
}
#endif

void
last_params(unsigned char *aa0, int n0, 
	    struct mngmsg *m_msg,
	    struct pstruct *ppst
#ifdef PCOMPLIB
	    , struct qmng_str *qm_msg
#endif
	    ) {
  int i, nsq;
  FILE *fp;

  if (n0 < 0) { return;}

  ppst->n0 = m_msg->n0;

  if (ppst->ext_sq_set) { nsq = ppst->nsqx; }
  else {nsq = ppst->nsq;}

/* currently, profiles are only available for SSEARCH, PRSS */
#if defined(SSEARCH) || defined(PRSS)

  ppst->pam2p[0] = alloc_pam2p(n0,nsq);
  ppst->pam2p[1] = alloc_pam2p(n0,nsq);

  if (ppst->pam_pssm) {
    if ((ppst->pgpfile_type == 0) && (fp=fopen(ppst->pgpfile,"rb"))) {
      read_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, 0, ppst);
      extend_pssm(aa0, n0, ppst);
    }
    else if ((ppst->pgpfile_type == 1) && (fp=fopen(ppst->pgpfile,"r"))) {
      read_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, 1, ppst);
      extend_pssm(aa0, n0, ppst);
    }
#if defined(SSEARCH) || (defined(PRSS) && !defined(FASTX))
    else if ((ppst->pgpfile_type == 2) && (fp=fopen(ppst->pgpfile,"rb"))) {
      if (read_asn_pssm(aa0, n0, ppst->nsq, ppst->pamscale, fp, ppst)>0) {
	extend_pssm(aa0, n0, ppst);
      }
      else {
	fprintf(stderr," Could not parse PSSM file: %s\n",ppst->pgpfile);
	ppst->pam_pssm = 0;
	return;
      }
    }
#endif
    else {
      fprintf(stderr," Could not open PSSM file: %s\n",ppst->pgpfile);
      ppst->pam_pssm = 0;
      return;
    }
  }
#endif

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  m_msg->nm0 = 1;
  for (i=0; i<n0; i++)
    if (aa0[i]==EOSEQ || aa0[i]==ESS) m_msg->nm0++;

/*
  for FASTS, we can do statistics in one of two different ways
  if there are <= 10 query fragments, then we calculate probabilistic
  scores for every library sequence.  If there are > 10 fragments, this
  takes much too long and too much memory, so we use the old fashioned
  raw score only z-score normalized method initially, and then calculate
  the probabilistic scores for the best hits.  To scale those scores, we
  also need a set of random probabilistic scores.  So we do the qshuffle
  to get them.

  For FASTF, precalculating probabilities is prohibitively expensive,
  so we never do it; FASTF always acts like FASTS with nfrags>10.

*/

#if defined(FASTS) || defined(FASTM)
  if (m_msg->nm0 > 10) m_msg->escore_flg = 0;
  else m_msg->escore_flg = 1;
#endif

  if (m_msg->escore_flg && (ppst->zsflag&1)) {
    m_msg->last_calc_flg = 0;
    m_msg->qshuffle = 0;
  }
  else {	/* need random query, second set of 2000 scores */
    m_msg->last_calc_flg = 1;
    m_msg->qshuffle = 1;
  }
#else
  m_msg->last_calc_flg = 0;
  m_msg->qshuffle = 0;
  m_msg->escore_flg = 0;
  m_msg->nm0 = 1;
#endif

/* adjust the ktup if appropriate */  

  if (!ktup_set && pgm_def_arr[ppst->pgm_id].ktup > 0) {
    if (m_msg->qdnaseq == SEQT_PROT) {
      ppst->param_u.fa.ktup = pgm_def_arr[ppst->pgm_id].ktup;
#if defined(FASTS) || defined(FASTM)
      if (n0 > 100) ppst->param_u.fa.ktup = 2;
#endif
      if (n0 < 40) ppst->param_u.fa.ktup = 1;
    }
    else if (m_msg->qdnaseq == SEQT_DNA || m_msg->qdnaseq == SEQT_RNA) {
      if (n0 < 20) ppst->param_u.fa.ktup = 1;
#if defined(FASTS) || defined(FASTM)
      /* with the current (April 12 2005) dropfs2.c - ktup cannot be > 2 */
      else ppst->param_u.fa.ktup = 2;
#else
      else if (n0 < 50) ppst->param_u.fa.ktup = 2;
      else if (n0 < 100)  ppst->param_u.fa.ktup = 3;
#endif
    }
  }

#ifdef PCOMPLIB
  qm_msg->nm0 = m_msg->nm0;
  qm_msg->escore_flg = m_msg->escore_flg;
  qm_msg->qshuffle = m_msg->qshuffle;
  qm_msg->pam_pssm = 0;
#endif
}

/* given a good profile in ppst->pam2p[0], make an extended profile
   in ppst->pam2p[1]
*/
void
extend_pssm(unsigned char *aa0, int n0, struct pstruct *ppst) {

  int i, j, nsq;
  int sa_x, sa_t, sa_b, sa_z;
  int **pam2p0, **pam2p1;

  nsq = ppst->nsq;

  pam2p0 = ppst->pam2p[0];
  pam2p1 = ppst->pam2p[1];

  sa_x = pascii['X'];
  sa_t = pascii['*'];
  sa_b = pascii['B'];
  sa_z = pascii['Z'];

  /* fill in boundaries, B, Z, *, X */
  for (i=0; i<n0; i++) {
    pam2p0[i][0] = -BIGNUM;
    pam2p0[i][sa_b] = (int)
      (((float)pam2p0[i][pascii['N']]+(float)pam2p0[i][pascii['D']]+0.5)/2.0);
    pam2p0[i][sa_z] = (int)
      (((float)pam2p0[i][pascii['Q']]+(float)pam2p0[i][pascii['E']]+0.5)/2.0);
    pam2p0[i][sa_x] = ppst->pam_xm;
    pam2p0[i][sa_t] = ppst->pam_xm;
  }

  /* copy pam2p0 into pam2p1 */
  for (i=0; i<n0; i++) {
    pam2p1[i][0] = -BIGNUM;
    for (j=1; j<=ppst->nsq; j++) {
      pam2p1[i][j] = pam2p0[i][j];
    }
  }

  /* then fill in extended characters, if necessary */
  if (ppst->ext_sq_set) {
    for (i=0; i<n0; i++) {
      for (j=1; j<=ppst->nsq; j++) {
	pam2p0[i][nsq+j] = pam2p0[i][j];
	pam2p1[i][nsq+j] = ppst->pam_xm;
      }
    }
  }
}
