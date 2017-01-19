/* Concurrent read version */

/* $Name: fa_34_26_5 $ - $Id: mw.h,v 1.20 2006/03/20 17:38:15 wrp Exp $ */

#include <sys/types.h>

#include "aln_structs.h"

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long fseek_t;
#else
typedef off_t fseek_t;
#endif
#endif

struct beststr {
  int n1;		/* sequence length */
  int *n1tot_p;		/* pointer (or NULL) to long sequence length */
  int score[3];		/* score */
  int sw_score;		/* do_walign() score */
  double comp;
  double H;
  double zscore;
  double escore;
  int segnum;
  int seglen;
  struct lmf_str *m_file_p;
  fseek_t lseek;
  char libstr[MAX_UID];
  int cont;
  int frame;
  int nsfnum;
  int sfnum[10];
  long loffset;
  struct a_struct aln_d;	/* these values are used by -m9 */
  struct a_res_str a_res;	/* need only a_res, not a_res[2], because different frames
				   for the same sequence are stored separately */
  int have_ares;
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


