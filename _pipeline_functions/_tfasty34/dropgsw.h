
/* global definitions shared by dropgsw.c and altivec.c */

/* definitions for SW */

struct f_struct {
  struct swstr *ss;
  struct swstr *f_ss, *r_ss;
  int *waa_s, *waa_a;
  int **pam2p[2];
  int *res;
  double aa0_f[MAXSQ];
  double *kar_p;
#if defined(SW_ALTIVEC) || defined(SW_SSE2)
  unsigned char      bias;
  unsigned short *   word_score;
  unsigned char *    byte_score;
  void *             workspace;
  int                alphabet_size;
  void *             word_score_memory;
  void *             byte_score_memory;
  void *             workspace_memory;
  int                try_8bit;
  int                done_8bit;
  int                done_16bit;
#endif
};

