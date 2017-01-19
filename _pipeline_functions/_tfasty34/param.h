/* $Name: fa_34_26_5 $ - $Id: param.h,v 1.41 2007/04/26 18:37:19 wrp Exp $ */


#ifndef P_STRUCT
#define P_STRUCT

#define MAXSQ 50


/* Concurrent read version */

struct fastr {
  int ktup;
  int cgap;
  int pgap;
  int pamfact;
  int scfact;
  int bestoff;
  int bestscale;
  int bkfact;
  int bktup;
  int bestmax;
  int altflag;
  int optflag;
  int iniflag;
  int optcut;
  int optcut_set;
  int optwid;
};

struct prostr {
    int gopen;
    int gextend;
    int width;
};

struct pstruct		/* parameters */
{
  int n0;	/* length of query sequence, used for statistics */
  int gdelval;	/* value gap open (-10) */
  int ggapval;	/* value for additional residues in gap (-2) */
  int gshift;	/* frameshift for fastx, fasty */
  int gsubs;	/* nt substitution in fasty */
  int p_d_mat;	/* dna match penalty */
  int p_d_mis;	/* dna mismatch penalty */
  int p_d_set;	/* using match/mismatch */
  int score_ix;	/* index to sorted score */
  int zsflag;	/* use scalebest() */
  int zsflag_f;	/* use scalebest() */
  int zs_win;
  int histint;		/* histogram interval */
  char sq[MAXSQ+1];
  int hsq[MAXSQ+1];
  int nsq;		/* length of normal sq */
  int ext_sq_set;	/* flag for using extended alphabet */
  char sqx[MAXSQ];
  int hsqx[MAXSQ+1];
  int c_nt[MAXSQ+1];
  int nsqx;	/* length of extended sq */
  int dnaseq;	/* -1 = not set (protein); 0 = protein; 1 = DNA; 2 = other, 3 RNA */
  int nt_align;	/* DNA/RNA alignment = 1 */
  int debug_lib;
  int tr_type;	/* codon table */
  int sw_flag;
  char pamfile[120];	/* pam file type */
  char pgpfile[120];
  int pgpfile_type;
  float pamscale;
  int pam_pssm;
  int pam_set;
  int have_pam2;
  int **pam2[2];
  int **pam2p[2];
  int pamoff;	/* offset for pam values */
  int pam_l, pam_h, pam_xx, pam_xm;	/* lowest, highest pam value */
  int pam_x_set;
  int pam_ms;		/* use a Mass Spec pam matrix */
  int maxlen;
  long zdb_size; /* force database size */
  int pgm_id;
  union {
    struct fastr fa;
    struct prostr pr;
  } param_u;
  int pseudocts;
  int shuff_node;
};

/* Result structure - do not remove */
struct rstruct
{
  int score[3];
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};

#ifndef PCOMPLIB
struct thr_str {
  int worker;
  void *status;
  int max_work_buf;
  int qframe;
  struct pstruct *ppst;
  int qshuffle;
  unsigned char *aa0;
  int n0;
  int nm0;
  int max_tot;
};

#include <sys/types.h>

/* this structure passes library sequences to the worker threads
   and returns scores */

struct buf_str {
  int n1;
  int *n1tot_p;
  unsigned char *aa1b;
#ifndef USE_FSEEKO
  long lseek;
#else
  off_t lseek;
#endif
  struct lmf_str *m_file_p;
  int cont;
  int qframe;
  int frame;
  int nsfnum;
  int sfnum[10];
  char libstr[20];	/* set to MAX_UID */
  struct rstruct rst;
  int r_score, qr_score;
  double r_escore, qr_escore;
};

struct buf_head {
  int buf_cnt;
  int have_results;
  unsigned char *start;
  struct buf_str *buf;
};

#endif

#endif  /* PSTRUCT */

#include "aln_structs.h"
