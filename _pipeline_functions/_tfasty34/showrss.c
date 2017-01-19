
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: showrss.c,v 1.12 2006/04/12 18:00:02 wrp Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#include "structs.h"
#include "param.h"

extern double
zs_to_E(double zs, int n1, int isdna, long entries,struct db_str db);
extern double zs_to_bit(double zs, int n0, int n1);
extern double zs_to_p(double zs);

extern double (*find_zp)(int score, double escore, int length, double comp, void *);

extern char *prog_func;

void showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	       struct beststr **bptr, int nbest, int qlib, struct mngmsg *m_msg,
	       struct pstruct pst, struct db_str db,
	       char *gstring2, void **f_str)
{
  double zs;
  int score;
  char *rlabel;
  struct beststr *bbp;

  if ((rlabel=strrchr(m_msg->label,' '))==NULL) rlabel = m_msg->label;

  fprintf(fp,"\n %s - %d shuffles; ",prog_func,m_msg->shuff_max);
  if (m_msg->shuff_wid > 0)
    fprintf(fp," window shuffle, window size: %d\n",m_msg->shuff_wid);
  else
    fprintf(fp," uniform shuffle\n");

  bbp = bptr[0];

  fprintf(fp," unshuffled %s score: %d;  bits(s=%d|n_l=%d): %4.1f p(%d) < %g\n",
	  rlabel,bbp->score[0],bbp->score[0], bbp->n1,
	  zs_to_bit(bbp->zscore,m_msg->n0,bbp->n1),bbp->score[0],zs_to_p(bbp->zscore));

  fprintf(fp,"For %ld sequences, a score >= %d is expected %4.4g times\n\n", 
	  pst.zdb_size,bbp->score[0],zs_to_E(bbp->zscore,bbp->n1,0l,pst.zdb_size,db)); 
}

void showalign (FILE *fp, unsigned char *aa0, unsigned char *aa1, int maxn,
		struct beststr **bptr, int nbest,int qlib, struct mngmsg m_msg,
		struct pstruct pst, void *f_str, char *gstring2)
{
}

void
aancpy(char *to, char *from, int count,
       struct pstruct pst)
{
  char *tp;

  tp=to;
  while (count-- && *from) {
    if (*from <= pst.nsq) *tp++ = pst.sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}
