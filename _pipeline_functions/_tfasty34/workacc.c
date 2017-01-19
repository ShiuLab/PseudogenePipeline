
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: workacc.c,v 1.19 2006/02/07 17:58:19 wrp Exp $ */

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "param.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

char err_str[128];

/* Initialization - set up defaults - assume protein sequence */
void w_init ()
{
  pascii=aascii;
}

#ifndef MPI_SRC
/* allocate memory for pam matrix - identical to version in initfa/sw.c */
alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;
  char err_str[128];

  if ((ppst->pam2[0] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     return -1;
  }

  if ((ppst->pam2[1] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     return -1;
  }

  if ((d2p = pam12 = (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     return -1;
   }
   for (i = 0; i < d1; i++, d2p += d2) ppst->pam2[0][i] = d2p;

   if ((d2p=pam12x= (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2d pam matrix: %d",d2);
     return -1;
   }
   for (i = 0;  i < d1; i++, d2p += d2) ppst->pam2[1][i] = d2p;

   return 1;
}

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

void
aancpy(char *to, char *from, int count, struct pstruct pst)
{
  char *tp, *sq;
  int nsq;

  tp=to;

  if (pst.ext_sq_set) {
    nsq = pst.nsqx;
    sq = pst.sqx;
  }
  else {
    nsq = pst.nsq;
    sq = pst.sq;
  }

  while (count-- && *from) {
    if (*from <= nsq) *tp++ = sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}
#endif

/* copies from from to to shuffling */

void
shuffle(unsigned char *from, unsigned char *to, int n)
{
  int i,j; unsigned char tmp;

  if (from != to) memcpy((void *)to,(void *)from,(size_t)n);

  for (i=n; i>0; i--) {
    j = nrand(i);
    tmp = to[j];
    to[j] = to[i-1];
    to[i-1] = tmp;
  }
  to[n] = 0;
}

/* this shuffle is for FASTS */
/* convert ',' -> '\0', shuffle each of the substrings */
qshuffle(unsigned char *aa0, int n0, int nm0)
{
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

/* copies from from to from shuffling */
void
wshuffle(unsigned char *from, unsigned char *to, int n, int wsiz, int *ieven)
{
  int i,j, k, mm; 
  unsigned char tmp, *top;

  memcpy((void *)to,(void *)from,n);
	
  mm = n%wsiz;

  if (*ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>0; i--) {
	j = nrand(i);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>0; i--) {
      j = nrand(i);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    *ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>0; i--) {
	j = nrand(i);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[0];
    for (i=mm; i>0; i--) {
      j = nrand(i);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    *ieven = 1;
  }
  to[n] = 0;
}

void initseq(char **seqc0, char **seqc0a, char **seqc1, char **seqca, int seqsiz)	/* initialize arrays */
{
  *seqc0=(char *)calloc((size_t)(seqsiz+1)*4,sizeof(char));
  *seqc0a= *seqc0+seqsiz+1;
  *seqc1= *seqc0a+seqsiz+1;
  *seqca= *seqc1+seqsiz+1;
  if (*seqc0==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
}

void freeseq(char **seqc0, char **seqc1, char **seqca)
{
  free(*seqc0);
}

#define ESS 49

void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
  seq[n]=0;
}
