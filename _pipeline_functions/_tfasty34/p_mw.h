/* Concurrent read version */

/* $Name: fa_34_26_5 $ - $Id: p_mw.h,v 1.17 2006/04/12 18:00:02 wrp Exp $ */

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long fseek_t;
#else
typedef off_t fseek_t;
#endif
#endif

struct beststr {
  int n1;		/* sequence number */
  int score[3];		/* score */
  int rscore;	/* score from shuffled sequence */
  int sw_score;	/* optimal score from alignment */
  double comp;	/* karlin 1/lambda comp.parameter */
  double H;	/* karlin H information content */
  double zscore;
  double escore;
  double r_escore;
  int segnum;
  int seglen;
  int  lib;
  fseek_t lseek;
  int cont;
  int frame;
  int m_seqnm;
  int seqnm;
  int wrkr;
  struct sql *desptr;
  struct a_struct *aln_d;
  char *aln_code;
  int aln_code_n;
  float percent, gpercent;
};

struct stat_str {
  int score;
  int n1;
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};

/* this structure passes library sequences to the worker threads
   and returns scores */

#include "w_mw.h"

/*
struct pbuf_head {
  int buf_cnt;
  unsigned char *start;
  struct sqs2 *buf;
};
*/
