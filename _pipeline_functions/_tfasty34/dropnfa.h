
/* global definitions shared by dropnfa.c and altivec.c */

#ifndef MAXSAV
#define MAXSAV 10
#endif



struct dstruct		/* diagonal structure for saving current run */
{			
   int     score;	/* hash score of current match */
   int     start;	/* start of current match */
   int     stop;	/* end of current match */
   struct savestr *dmax;   /* location in vmax[] where best score data saved */
};

struct savestr
{
   int     score;		/* pam score with segment optimization */
   int     score0;		/* pam score of best single segment */
   int     gscore;		/* score from global match */
   int     dp;			/* diagonal of match */
   int     start;		/* start of match in lib seq */
   int     stop;		/* end of match in lib seq */
};

struct bdstr { int CC, DD, CP, DP;};

struct f_struct {
  struct dstruct *diag;
  struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
  struct savestr *vptr[MAXSAV];
  struct savestr *lowmax;
  int ndo;
  int noff;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  int *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#ifdef TFASTA
  unsigned char *aa1x;
  int n10;
#endif
  struct bdstr *bss;
  struct swstr *ss;
  struct swstr *f_ss, *r_ss;
  int *waa0;
  int *waa1;
  int *res;
  int max_res;
  double aa0_f[MAXSQ];
  double *kar_p;
 
#ifdef FA_ALTIVEC
  int vec_len;
  vecInt **vec_matrix;
  vector signed ALTIVEC_SIZE *vec_HH;
  vector signed ALTIVEC_SIZE *vec_EE;

  int vec_len2;
  vecInt2 **vec_matrix2;
  vector signed ALTIVEC_SIZE2 *vec_HH2;
  vector signed ALTIVEC_SIZE2 *vec_EE2;
#endif
};

static int
FLOCAL_ALIGN(const unsigned char *A, const unsigned char *B,
	     int M, int N, int low, int up,
	     int **W, int G,int H, int MW,
	     struct f_struct *f_str);
