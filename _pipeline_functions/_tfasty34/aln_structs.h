
#ifndef A_STRUCT
#define A_STRUCT

struct a_struct {
  int smin0;		/* coordinate of display start in seqc0 */
  int smin1;		/* coordinate of display start in seqc1 */
  int amin0, amax0;	/* coordinate of alignment start in seqc0 */
  int amin1, amax1;	/* coordinate of alignment start in seqc1 */

  int llen;
  int llcntx, llcntx_flg, showall;

  int qlrev, qlfact;
  int llrev, llfact, llmult;
  int frame;

  int a_len;			/* consensus alignment length */
  int nident, nsim, ngap_q, ngap_l, nfs;	/* number of identities, gaps in q, l */
  long d_start0,d_stop0;
  long d_start1,d_stop1;
};

struct a_res_str {
  int min0, max0;	/* boundaries of alignment in aa0 */
  int min1, max1;	/* boundaries of alignment in aa1 */
  int *res;
  int nres;
};
#endif
