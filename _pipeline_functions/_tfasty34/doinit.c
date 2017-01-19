/*	doinit.c	general and function-specific initializations */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: doinit.c,v 1.62 2007/01/08 15:38:46 wrp Exp $ */

/* this file performs general initializations of search parameters

   In addition, it calls several functions in init??.c that provide
   program-specific initializations:

   f_initenv()	- called from initenv()
   f_getopt()	- called from initenv() during a getopt() scan
   f_getarg()	- called from initenv() after the getopt() scan

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "param.h"
#include "upam.h"	/* required for 'U' option change of nascii */

#include "structs.h"

#define XTERNAL
#include "uascii.h"
#undef XTERNAL

extern char *s_optstr;
extern int optind;		/* used by getopt() */

#ifdef PCOMPLIB
#define PARALLEL
#include "p_mw.h"
extern char pgmdir[];
extern char managepgm[];
extern char workerpgm[];
extern int max_buf_cnt;
#define MAX_WORKERS MAXWRKR
#endif

char prog_name[MAX_FN];

extern void f_initenv(struct mngmsg *, struct pstruct *, unsigned char **);
extern void f_lastenv(struct mngmsg *, struct pstruct *);
extern void f_getopt(char, char *, struct mngmsg *, struct pstruct *);
extern void f_getarg(int, char **, int, struct mngmsg *, struct pstruct *);
void ann_ascii(int *qascii, char *ann_arr);
int set_markx(int markx, int val);

int optcnt;
int max_workers=MAX_WORKERS;
#ifdef PCOMPLIB
int worker_1=0;
int worker_n=0;
#endif
extern char *optarg;

/* initenv ()  initializes the environment */
void initenv (int argc, char **argv, struct mngmsg *m_msg, 
		 struct pstruct *ppst, unsigned char **aa0)
{
   char   *cptr, ctmp;
   int     copt, itmp;

   /* options for all search functions */
   char   *g_optstr = "ab:BC:d:DE:F:HiIJ:K:l:Lm:M:N:O:QqR:T:v:V:w:W:X:z:Z:";
   char    optstring[MAX_STR];

/*  these initializations will be used by all functions */

   /* prog_name[] is only used for error messages */
   strncpy(prog_name,argv[0],sizeof(prog_name));
   prog_name[sizeof(prog_name)-1]='\0';

#ifdef PARALLEL
   if ((cptr = getenv ("MANAGEPGM")) != NULL) strncpy (managepgm, cptr, 120);
   if ((cptr = getenv ("WORKERPGM")) != NULL) strncpy (workerpgm, cptr, 120);
   if ((cptr = getenv ("PGMDIR")) != NULL) strncpy (pgmdir, cptr, 120);
#endif

   m_msg->ltitle[0] = '\0';

   if ((cptr=getenv("FASTLIBS"))!=NULL) {
     strncpy(m_msg->flstr,cptr,MAX_FN);
     m_msg->flstr[MAX_FN-1] = '\0';
   }
   else m_msg->flstr[0]='\0';

   m_msg->hist.hist_a = NULL;
   m_msg->outfile[0] = '\0';
   m_msg->ldnaseq = SEQT_PROT;	/* library is protein */
   m_msg->n1_low = 0;
   m_msg->n1_high = BIGNUM;
   m_msg->ql_start = 1;	/* start with first query sequence */
   m_msg->ql_stop = BIGNUM;	/* end with the last query sequence */

   m_msg->pamd1 = MAXSQ;
   m_msg->pamd2 = MAXSQ;

   m_msg->term_code = 0;
   ppst->tr_type = 0;
   ppst->debug_lib = 0;
   m_msg->nshow = 20;
#if defined(PCOMPLIB)
   m_msg->nohist = 1;
   m_msg->mshow = 20;
#else
   m_msg->nohist = 0;
   m_msg->mshow = 50;
#endif
   m_msg->ashow = -1;
   m_msg->nmlen = DEF_NMLEN;
   m_msg->z_bits = 1;
   m_msg->mshow_flg = 0;
   m_msg->aln.llen = 0;
   m_msg->aln.llcntx = 30;
   m_msg->aln.llcntx_flg = 0;
   m_msg->e_cut = 10.0;
   m_msg->e_low = 0.0;
   m_msg->e_cut_set = 0;
   m_msg->revcomp = 0;
   m_msg->self = 0;
   m_msg->long_info = 0;
   m_msg->maxn = 0;
   m_msg->dupn = SEQDUP;
   m_msg->dfile[0] = '\0';
   m_msg->tname[0] = '\0';
   m_msg->lname[0] = '\0';
   m_msg->show_code = 0;
   m_msg->aln.showall = 0;
   m_msg->markx = 0;
   m_msg->sq0off = m_msg->sq1off = 1;
   strncpy(m_msg->sqnam,"aa",4);
   strncpy(m_msg->sqtype,"protein",10);
   m_msg->ann_flg = 0;
   m_msg->ann_arr[0] = '\0';
   m_msg->aa0a = NULL;
   
   ppst->zsflag = ppst->zsflag_f = 1;
   ppst->zs_win = 0;

   ppst->zdb_size = -1;
   ppst->dnaseq = SEQT_PROT;	/* default is protein */
   ppst->nt_align = 0;

   f_initenv (m_msg, ppst, aa0);

   strncpy (optstring, g_optstr, sizeof (optstring));
   strncat (optstring, s_optstr, sizeof (optstring));

   while ((copt = getopt (argc, argv, optstring)) != EOF)
   {
      if (strchr (g_optstr, copt) != NULL)
      {
	switch (copt) {  /* switches for all options */
	case 'a': m_msg->aln.showall = 1; break;
	case 'b':
	  if (optarg[0] == '$') {
	    m_msg->mshow = -1;
	    m_msg->e_cut = 10000000.0;
	    break;
	  }
	  else sscanf (optarg, "%d", &m_msg->mshow);
	  m_msg->e_cut = 10000000.0;
	  m_msg->e_cut_set = 1;
	  m_msg->mshow_flg = 1;
	  break;
	case 'B': m_msg->z_bits = 0; break;
	case 'C': sscanf(optarg,"%d",&m_msg->nmlen);
	  if (m_msg->nmlen > MAX_UID-1) m_msg->nmlen = MAX_UID-1;
	  break;
	case 'd': sscanf(optarg,"%d",&m_msg->ashow);
	  if (m_msg->ashow > m_msg->mshow) m_msg->mshow=m_msg->ashow;
	  /* m_msg->ashow_flg = 1; (ashow_flg not in structs.h, not used)*/
	  break;
	case 'D': ppst->debug_lib = 1;
	  break;
	case 'E':
	  sscanf(optarg,"%lf",&m_msg->e_cut);
	  m_msg->e_cut_set = 1;
	  break;
	case 'F':
	  sscanf(optarg,"%lg",&m_msg->e_low);
	  m_msg->e_cut_set = 1;
	  break;
	case 'H':
#if defined(PCOMPLIB)
	  m_msg->nohist = 0; break;
#else
	  m_msg->nohist = 1; break;
#endif
	case 'i':
	  m_msg->revcomp = 1; break;
#ifdef PARALLEL
	case 'I':
	  m_msg->self = 1; break;
	case 'J':
	  if (optarg[0]==':') {
	    m_msg->ql_start = 0;
	    sscanf(optarg,":%d",&m_msg->ql_stop);
	    m_msg->ql_stop++;
	  }
	  else if (!strchr(optarg,':')) {
	    m_msg->ql_stop = BIGNUM;
	    sscanf(optarg,"%d",&m_msg->ql_start);
	  }
	  else {
	    sscanf(optarg,"%d:%d",&m_msg->ql_start,&m_msg->ql_stop);
	    m_msg->ql_stop++;
	  }
	  break;
	case 'K':
	  sscanf(optarg,"%d",&max_buf_cnt);
	  break;
#endif
	case 'l':
	  strncpy(m_msg->flstr,optarg,MAX_FN);
	  m_msg->flstr[MAX_FN-1]='\0';
	  break;
	case 'L':
	  m_msg->long_info = 1; break;
	case 'm':
	  sscanf(optarg,"%d%c",&itmp,&ctmp);
	  if (itmp==9 && ctmp=='c') {
	    m_msg->show_code = SHOW_CODE_ALIGN;
	  }
	  else if (itmp==9 && ctmp=='i') {
	    m_msg->show_code = SHOW_CODE_ID;
	  }
	  if (itmp > 6 && itmp != 10 && itmp != 9) itmp = 0;
	  m_msg->markx = set_markx(m_msg->markx,itmp);
	  break;
	case 'M':
	  sscanf(optarg,"%d-%d",&m_msg->n1_low,&m_msg->n1_high);
	  if (m_msg->n1_low < 0) {
	    m_msg->n1_high = -m_msg->n1_low;
	    m_msg->n1_low = 0;
	  }
	  if (m_msg->n1_high == 0) m_msg->n1_high = BIGNUM;
	  if (m_msg->n1_low > m_msg->n1_high) {
	    fprintf(stderr," low cutoff %d greater than high %d\n",
		    m_msg->n1_low, m_msg->n1_high);
	    m_msg->n1_low = 0;
	    m_msg->n1_high = BIGNUM;
	  }
	  break;
	case 'N':
	  sscanf(optarg,"%d",&m_msg->maxn);
	  break;
	case 'p':
	  m_msg->qdnaseq = SEQT_PROT;
	  ppst->dnaseq = SEQT_PROT;
	  strncpy(m_msg->sqnam,"aa",4);
	  break;
	case 'O':
	  strncpy(m_msg->outfile,optarg,MAX_FN);
	  m_msg->outfile[MAX_FN-1]='\0';
	  break;
	case 'q':
	case 'Q':
	  m_msg->quiet = 1;
	  break;
	case 'R':
	  strncpy (m_msg->dfile, optarg, MAX_FN);
	  m_msg->dfile[MAX_FN-1]='\0';
	  break;
	case 'T':
#ifdef PCOMPLIB
	  if (strchr(optarg,'-') != NULL) {
	    sscanf(optarg,"%d-%d",&worker_1,&worker_n);
	    if (worker_1 > worker_n) {
	      worker_1 = worker_n = 0;
	    }
	  }
	  else 
#endif
	    sscanf (optarg, "%d", &max_workers);
	  if (max_workers < 0) max_workers=1;
	  break;
	case 'v':
	  sscanf (optarg,"%d",&ppst->zs_win);
	  break;
	case 'V':
	  strncpy(m_msg->ann_arr+1,optarg,MAX_FN-2);
	  m_msg->ann_arr[0]='\0';
	  m_msg->ann_arr[MAX_FN-2]='\0';
	  m_msg->ann_flg = 1;
	  ann_ascii(qascii, m_msg->ann_arr);
	  break;
/*
	case 'V':
	  fprintf(stderr," -V option not currently supported in parallel\n");
	  break;
*/
	case 'w':
	  sscanf (optarg,"%d",&m_msg->aln.llen);
	  if (m_msg->aln.llen < 10) m_msg->aln.llen = 10;
	  if (m_msg->aln.llen > 200) m_msg->aln.llen = 200;
	  if (!m_msg->aln.llcntx_flg) m_msg->aln.llcntx = m_msg->aln.llen/2;
	  break;
	case 'W':
	  sscanf (optarg,"%d",&m_msg->aln.llcntx);
	  m_msg->aln.llcntx_flg = 1;
	  break;
	case 'X':
	  sscanf (optarg,"%ld %ld",&m_msg->sq0off,&m_msg->sq1off); break;
	case 'z':
	  sscanf(optarg,"%d",&ppst->zsflag);
	  break;
	case 'Z':
	  sscanf(optarg,"%ld",&ppst->zdb_size);
	  break;
	}
      }
      else if (strchr (s_optstr, copt))
	 f_getopt (copt, optarg, m_msg, ppst);
   }
   optind--;

   f_lastenv (m_msg, ppst);

   if (argc - optind < 3) return;
   m_msg->tnamesize = sizeof (m_msg->tname);
   if (argc - optind > 1) strncpy (m_msg->tname, argv[optind + 1],MAX_FN);
   if (argc - optind > 2) { strncpy(m_msg->lname, argv[optind + 2],MAX_FN); }
   f_getarg (argc, argv, optind, m_msg, ppst);
}

int
ann_scan(unsigned char *aa0, int n0, struct mngmsg *m_msg, int seqtype)
{
  unsigned char *aa0p, *aa0d, *aa0ad;
  int n_n0;

  /* count how many "real" residues */

  if (seqtype==SEQT_UNK) {
    for (n_n0=0, aa0p = aa0; aa0p < aa0+n0; aa0p++) {
      if (*aa0p > '@' || *aa0p == ESS ) n_n0++;
    }
  }
  else {
    for (n_n0=0, aa0p = aa0; aa0p < aa0+n0; aa0p++) {
      if (*aa0p < NANN ) n_n0++;
    }
  }

  aa0d = aa0;
  /* n_n0 has the real sequence length */
  if ((m_msg->aa0a = calloc(n_n0+2, sizeof(char)))==NULL) {
    fprintf(stderr," cannot allocate annotation sequence: %d\n",n_n0);
    m_msg->ann_flg = 0;
    if (seqtype==SEQT_UNK) {
      for (aa0p = aa0; aa0p < aa0+n0; aa0p++) {
	if (*aa0p > '@' || *aa0p == ESS) {*aa0d++ = *aa0p;}
      }
    }
    else {
      for (aa0p = aa0; aa0p < aa0+n0; aa0p++) {
	if (*aa0p < NANN) {*aa0d++ = *aa0p;}
      }
    }
      *aa0d = '\0';
    return n_n0;
  }

  aa0ad = m_msg->aa0a;
  if (seqtype==SEQT_UNK) {
    for (aa0p = aa0; aa0p<aa0+n0; aa0p++) {
      if (*aa0p > '@' || *aa0p == ESS) {*aa0d++ = *aa0p; *aa0ad++='\0';}
      else if (aa0ad > m_msg->aa0a) { aa0ad[-1] = *aa0p - NANN;}
    }
  }
  else {
    for (aa0p = aa0; aa0p<aa0+n0; aa0p++) {
      if (*aa0p < NANN) {*aa0d++ = *aa0p; *aa0ad++='\0';}
      else if (aa0ad > m_msg->aa0a) { aa0ad[-1] = *aa0p - NANN;}
    }
  }
  *aa0ad = *aa0d = '\0';
  return n_n0;
}

void
ann_ascii(int *qascii, char *ann_arr)
{
  char *ann_p;
  int ann_ix = NANN+1;

  ann_arr[0] = ' ';
  if (strchr(ann_arr+1,'*')) {qascii['*'] = NA;}

  for (ann_p = ann_arr+1; *ann_p; ann_p++) {
    if (qascii[*ann_p] == NA) { qascii[*ann_p] = ann_ix++;}
  }
}

int 
set_markx(int markx, int val) {

  if (val < 3) {
    return markx | (MX_ATYPE & val);
  }
  else if (val == 3) {
    markx |= (MX_ATYPE + MX_ASEP);
  }
  else if (val == 4) {
    markx |= (MX_ATYPE + MX_AMAP);
  }
  else if (val == 5) {
    markx |= MX_AMAP;
  }
  else if (val == 6) {
    markx |= (MX_HTML) ;
  }
  else if (val == 9) {
    markx |= MX_M9SUMM;
  }
  else if (val == 10) {
    markx |= MX_M10FORM;
  }

  return markx;
}
