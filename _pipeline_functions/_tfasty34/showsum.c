/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: showsum.c,v 1.21 2006/06/22 15:00:51 wrp Exp $ */

/* 10 December 1999 --

   code modified to reflect the fact that there may be two scores for
   each sequence - e.g. forward and reverse strand - and only one of them
   - presumably the best - is a related score.
*/

/* showsum.c should report statistics for success -

   given the sorted results

   (1) find the highest scoring unrelated sequence: unf_score0
   	find the number of related sequences missed: relm_num0
   (2) find the 0.5% highest scoring unrelated sequence: unf_score05
   	find the number of related sequences missed: relm_num05
   (3) find the score where the number of related sequences
   	missed and the number of unrelated sequences found
	matches; report the score and the number: equ_score, equ_num;

The query sequence library number will be put in qsfnum.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#include "structs.h"

#ifndef SFCHAR
#define SFCHAR ':'
#define NSFCHAR '!'
#endif

#ifdef PCOMPLIB
#define BSFNUM(i) bptr[i]->desptr->sfnum
#define QSFNUM qsfnum
#define NQSFNUM qsfnum_n
#else
#define BSFNUM(i) bptr[i]->sfnum
#define QSFNUM m_msg->qsfnum
#define NQSFNUM m_msg->qsfnum_n
#endif

#define MAX_BLINE 200

double E_to_zs(double, long);
double zs_to_E(double,int,int,long,struct db_str db);
double zs_to_bit(double,int,int);
#ifdef PVM_SRC
void sf_sort(int *s, int n);
#endif
void lnum_sort(struct beststr **s, int n);

void showbest (FILE *fp, 
#ifndef PCOMPLIB
	       unsigned char **aa0, unsigned char *aa1, int maxn,
#endif
	       struct beststr **bptr,int nbest,
	       int qlib, struct mngmsg *m_msg, struct pstruct pst,
	       struct db_str db,
	       char *gstring2
#ifndef PCOMPLIB
	       ,void *f_str
#endif
	       )
{
  int i, j, k, rel_tot;
  int irelv;

  int unf_num0, relm_num0;
  int unf_num01,relm_num01;
  int unf_num02, relm_num02;
  int unf_num05, relm_num05;
  int unf_num100, relm_num100;
  int equ_num, rel_3_num, rel_1_num;

  double unf_score0, unf_score01, unf_score02 ,unf_score05;
  double unf_score100, equ_score, rel_3_score, rel_1_score;
  double unf_score0_b, unf_score01_b, unf_score02_b ,unf_score05_b;
  double unf_score100_b, equ_score_b, rel_3_score_b, rel_1_score_b;
  char *bp;

#ifdef PCOMPLIB
  int qsfnum[10],qsfnum_n[10],isf,nsf,nsf_n;
  char  *bp1, *bpn, *tp;
  char sfstr[MAX_FN];
#endif

#ifdef PCOMPLIB
  /*	not done here because done in pvcomplib.c */
  if ((bp=strchr(m_msg->qtitle,SFCHAR))!=NULL) {
    strncpy(sfstr,bp+1,sizeof(sfstr));
    sfstr[sizeof(sfstr)-1]='\0';
    if ((bp1=strchr(sfstr,SFCHAR)) != NULL) { /* look for second | */
      if ((bpn=strchr(sfstr,NSFCHAR))!=NULL) *bpn = '\0';
      *bp1='\0';
      tp = strtok(sfstr," \t");
      qsfnum[0]=atoi(tp);
      isf = 1;
      while ((tp=strtok(NULL," \t"))!=NULL) {
	qsfnum[isf++] = atoi(tp);
	if (isf >= 10) {
	  fprintf(stderr," error - too many superfamilies: %d\n %s\n",
		  isf,m_msg->qtitle);
	  break;
	}
      }
      qsfnum[nsf=isf]=0;
      sf_sort(qsfnum,nsf);

      /* now get negatives */
      qsfnum_n[0]= nsf_n = 0;
      if (bpn != NULL) {
	tp = strtok(bpn+1," \t");
	qsfnum_n[0]=atoi(tp);
	isf = 1;
	while ((tp=strtok(NULL," \t"))!=NULL) {
	  qsfnum_n[isf++] = atoi(tp);
	  if (isf >= 10) {
	    fprintf(stderr,
		    " error - too many negative superfamilies: %d\n %s\n",
		    isf,m_msg->qtitle);
	    break;
	  }
	}
	qsfnum[nsf_n=isf]=0;
	sf_sort(qsfnum_n,nsf_n);
      }
    }
    else {	/* only one sfnum */
      sscanf(bp+1,"%d",qsfnum);
      qsfnum[1]=0;
      qsfnum_n[0]= nsf_n = 0;
    }
  }
  else {
    fprintf(stderr," no query superfamily number\n %s\n",m_msg->qtitle);
    return;
  }
#endif

  if (m_msg->qframe > 1 || m_msg->nframe > 1) {

    /* this code is included for cases where there are several scores -
       forward and reverse, or six in the case of tfastf33s, for each
       sequence

       lnum_sort sorts the library by lseek position, which will be
       the same for the same sequence
    */

    lnum_sort(bptr,nbest);

  /* merge, saving the best score */
    i = j = 0;

    /* i has the source position we are currently examining
       k has the adjacent alternative scores ( k > i) 
       j has the destination 
    */

    while (i<nbest) {
      for (k=i+1; k < nbest && bptr[i]->lseek == bptr[k]->lseek; k++) {
	if (bptr[i]->zscore < bptr[k]->zscore) bptr[i] = bptr[k];
      }
      bptr[j++]=bptr[i];
      i = k;
    }

    if (j != m_msg->nbr_seq) {
      fprintf(stderr,"*** warning ***, nbest (%d/%d) != nbr_seq (%d)\n",
	      j,nbest,m_msg->nbr_seq);
      fprintf(stdout,"*** warning ***, nbest (%d/%d) != nbr_seq (%d)\n",
	      j,nbest,m_msg->nbr_seq);
    }
    nbest = j;

    if (pst.zsflag >=0) sortbeste(bptr, nbest);
    else sortbest(bptr,nbest,pst.score_ix);
  }

/* fprintf(stderr," %1d label is %s (%s)\n",irelv,labptr,label); */

/* get the query superfamily */
  
  for (i=0; i<nbest; i++) {
    /*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)==0 && sfn_cmp(BSFNUM(i),NQSFNUM)==0) {
      unf_num0=i;
      unf_score0=bptr[i]->zscore;
      unf_score0_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);
      break;
    }
  }

  if (i>=nbest) {
    fprintf(stderr," %s: %d\n error - no unrelated sequences\n",
	    m_msg->qtitle,QSFNUM[0]);
    return;
  }
  
  for (i=rel_tot=relm_num0=0; i<nbest; i++) {
    /*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)>0 ) {
      rel_tot++;			/* total related */
      if (bptr[i]->zscore <= unf_score0) relm_num0++;
#ifdef DEBUG      
      if (pst.debug_lib)
	fprintf(stderr,"%d\t%l\t%.1f\n",i,bptr[i]->lseek,bptr[i]->zscore);
#endif
    }
  }
  
  /* relm_num0, unf_num0, unf_score0 done */
  
  /* now calculate number missed at various expectation value cutoffs */
  /* calculate z-score cutoff for E()=0.01, 0.02, 0.05 */

  unf_score01 = E_to_zs(0.01,db.entries);
  unf_score02 = E_to_zs(0.02,db.entries);
  unf_score05 = E_to_zs(0.05,db.entries);
  unf_score100 = E_to_zs(1.00,db.entries);

  /* relm_num01, unf_num01, unf_score01 done */
  
  for (i=unf_num01=0,relm_num01=rel_tot;
       i<nbest && bptr[i]->zscore >= unf_score01; i++) {
/*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)==0) {
      if (sfn_cmp(BSFNUM(i),NQSFNUM)==0) unf_num01++;
    }
    else relm_num01--;
  }
  unf_score01_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);

  for (i=unf_num02=0,relm_num02=rel_tot;
       i<nbest && bptr[i]->zscore >= unf_score02; i++) {
/*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)==0) {
      if (sfn_cmp(BSFNUM(i),NQSFNUM)==0) unf_num02++;
    }
    else relm_num02--;
  }
  unf_score02_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);
      
  for (i=unf_num05=0,relm_num05=rel_tot;
       i<nbest && bptr[i]->zscore >= unf_score05; i++) {
/*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)==0) {
      if (sfn_cmp(BSFNUM(i),NQSFNUM)==0) unf_num05++;
    }
    else relm_num05--;
  }
  unf_score05_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);
      
  for (i=unf_num100=0,relm_num100=rel_tot;
       i<nbest && bptr[i]->zscore >= unf_score100; i++) {
/*     if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)==0) {
      if (sfn_cmp(BSFNUM(i),NQSFNUM)==0) unf_num100++;
    }
    else relm_num100--;
  }
  unf_score100_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);
      
  /* the final criterion finds the score and the number of sequences
     where the number of unrelated sequences found == the number of
     related sequences missed. */
  
  equ_num=0;
  i = 0; j=nbest-1;

/* j is counting up the list of scores (actually down the array) from
  the lowest scoring related sequence

  i is counting down the list of scores (actually up the array)
  from the highest scoring unrelated sequence */

  for (i=0, j=nbest-1; j>=0 && i<nbest; i++,j--) {
    /* i++ while sequences are related, stop at next unrelated */
    while (i<nbest && (sfn_cmp(BSFNUM(i),QSFNUM) || sfn_cmp(BSFNUM(i),NQSFNUM))) i++; 
    /* j-- while sequences are unrelated, stop at next related */
    while (j>=0 && ( sfn_cmp(BSFNUM(j),QSFNUM)==0)) j--;
    /*
      fprintf(stderr,"i: %3d %3d %4d; j: %3d %3d %4d\n",i,bptr[i]->zscore,
      BSFNUM(i),j,bptr[j]->zscore,BSFNUM(j));
      */
    /* if unrelated [i] score <= related [j] score, quit */
    if (bptr[i]->zscore <= bptr[j]->zscore) break;
    equ_num++;
  }
  
  equ_score = 0.0;
  if (i>=nbest || j<0) {
#ifndef PCOMPLIB
    if (pst.debug_lib) 
#endif
      fprintf(stderr," i (%3d), j (%3d) off end\n %s\n", i, j,m_msg->qtitle);
    equ_num = rel_tot+1; equ_score = 0.0;
  }
  else {
    equ_score=bptr[i]->zscore;
    equ_score_b =zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);
  }
  
  /* get the lowest scoring related */
  for (i=0,rel_1_num=rel_tot-1; i<nbest && rel_1_num > 0; i++) {
/*    if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)) rel_1_num--;
  }
  rel_1_num = i;
  rel_1_score = bptr[i]->zscore;
  rel_1_score_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);

  /* get the 3rd lowest scoring related */
  for (i=0,rel_3_num=rel_tot-3; i<nbest && rel_3_num > 0; i++) {
/*     if (sfn_cmp(BSFNUM(i),NQSFNUM)) continue; */
    if (sfn_cmp(BSFNUM(i),QSFNUM)) rel_3_num--;
  }
  rel_3_num = i;
  rel_3_score = bptr[i]->zscore;
  rel_3_score_b=zs_to_bit(bptr[i]->zscore,m_msg->n0,bptr[i]->n1);

  fprintf(fp,"%3d>%s - %d (%d/%d)\n",
	  qlib,m_msg->qtitle, QSFNUM[0],rel_tot,nbest);
  fprintf(fp," 0.0 criterion- relm: %3d pos: %3d score: %5.1f exp: %6.4g\n",
	  relm_num0, unf_num0+1, unf_score0_b,
	  zs_to_E(unf_score0,m_msg->n0,pst.dnaseq,pst.zdb_size,db));
  fprintf(fp," 0.01 criterion- relm: %3d unf: %3d score: %5.1f exp: %6.4g\n",
	  relm_num01, unf_num01, unf_score01_b,
	  zs_to_E(unf_score01,m_msg->n0,pst.dnaseq,pst.zdb_size,db));
  fprintf(fp," 0.02 criterion- relm: %3d unf: %3d score: %5.1f exp: %6.4g\n",
	  relm_num02, unf_num02, unf_score02_b,
	  zs_to_E(unf_score02,m_msg->n0,pst.dnaseq,pst.zdb_size,db));
  fprintf(fp," 0.05 criterion- relm: %3d unf: %3d score: %5.1f exp: %6.4g\n",
	  relm_num05, unf_num05, unf_score05_b,
	  zs_to_E(unf_score05,m_msg->n0,pst.dnaseq,pst.zdb_size,db));
  fprintf(fp," 1.00 criterion- relm: %3d unf: %3d score: %5.1f exp: %6.4g\n",
	  relm_num100, unf_num100, unf_score100_b,
	  zs_to_E(unf_score100,m_msg->n0,pst.dnaseq,pst.zdb_size,db));

  fprintf(fp," equ num: %3d score: %5.1f exp: %6.4g\n",equ_num,equ_score_b,
	  zs_to_E(equ_score,m_msg->n0,pst.dnaseq,pst.zdb_size,db));

  fprintf(fp," rel[-1]: %3d score: %5.1f exp: %6.4g\n",rel_1_num+1,rel_1_score_b,
	  zs_to_E(rel_1_score,m_msg->n0,pst.dnaseq,pst.zdb_size,db));
  fprintf(fp," rel[-3]: %3d score: %5.1f exp: %6.4g\n",rel_3_num+1,rel_3_score_b,
	  zs_to_E(rel_3_score,m_msg->n0,pst.dnaseq,pst.zdb_size,db));

  /* 
  fprintf(fp,"/ ** %s ** /\n",gstring2);
  fflush(fp);
  */
  m_msg->nshow = m_msg->ashow;
}

#ifdef PCOMPLIB
void showalign()
{}

#if !defined(MPI_SRC) && !defined(PCOMPLIB)
void
sf_sort(int *s, int n)
{
  int gap, i, j;
  int itmp;
	
  for (i=0; i<n-1; i++)
    if (s[i]>s[i+1]) goto l2;
  return;

l2:
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if (s[j] <= s[j+gap]) break;
	itmp = s[j];
	s[j]=s[j+gap];
	s[j+gap]=itmp;
      }
}

#endif
#endif

void
lnum_sort(struct beststr **s, int n)
{
  int gap, i, j;
  struct beststr *btmp;
	
  for (i=0; i<n-1; i++)
    if (s[i]->lseek > s[i+1]->lseek) goto l2;
  return;

l2:
  for (gap=n/2; gap>0; gap/=2)
    for (i=gap; i<n; i++)
      for (j=i-gap; j>=0; j -= gap) {
	if (s[j]->lseek <= s[j+gap]->lseek) break;
	btmp = s[j];
	s[j]=s[j+gap];
	s[j+gap]=btmp;
      }
}  

#ifdef MPI_SRC
void
aancpy(char *to, char *from, int count, struct pstruct pst)
{
  char *tp, *sq;
  int nsq;

  if (pst.ext_sq_set) {
    nsq = pst.nsqx;
    sq = pst.sqx;
  }
  else {
    nsq = pst.nsq;
    sq = pst.sq;
  }

  tp=to;
  while (count-- && *from) {
    if (*from <= nsq) *tp++ = sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}
#endif
