#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "param.h"
#include "tatstats.h"

#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

/* calc_priors() - calculate frequencies of amino-acids, possibly with counts */
/* generate_tatprobs() - build the table of score probabilities if the
   sequences are not too long */

double
det(double a11, double a12, double a13,
    double a21, double a22, double a23,
    double a31, double a32, double a33);

double power(double r, int p)
{
  double tr;
  int neg;

  if (r==0.0) return(p==0?1.0:0.0);
  if (neg = p<0) p = -p;
  tr = 1.0;
  while (p>0) {
    if (p & 1) tr *= r;
    p >>= 1;
    if (p) r *= r;
  }
  return((neg? 1.0/tr: tr));
}

double
factorial (int a, int b) {

  double res = 1.0;

  if(a == 0) { return 1.0; }

  while(a > b) {
    res *= (double) a;
    a--;
  }

  return res;
}

void
calc_priors(double *priors,
	    struct pstruct *ppst,
	    struct f_struct *f_str,
	    const unsigned char *aa1, int n1,
	    int pseudocts)
{
  long counts[25], sum;
  int i;

  if(n1 == 0 && f_str->priors[1] > 0.0) {
    for(i = 1 ; i <= ppst->nsq ; i++) {
      priors[i] = f_str->priors[i];
    }
    return;
  }

  if(n1 == 0) {
    if (ppst->dnaseq==SEQT_PROT ) {

      /* Robinson & Robinson residue counts from Stephen Altschul */
      counts[ 1] = 35155;    /*   A  */
      counts[ 2] = 23105;    /*   R  */  
      counts[ 3] = 20212;    /*   N  */  
      counts[ 4] = 24161;    /*   D  */  
      counts[ 5] =  8669;    /*   C  */  
      counts[ 6] = 19208;    /*   Q  */  
      counts[ 7] = 28354;    /*   E  */  
      counts[ 8] = 33229;    /*   G  */  
      counts[ 9] =  9906;    /*   H  */  
      counts[10] = 23161;    /*   I  */  
      counts[11] = 40625;    /*   L  */  
      counts[12] = 25872;    /*   K  */  
      counts[13] = 10101;    /*   M  */  
      counts[14] = 17367;    /*   F  */  
      counts[15] = 23435;    /*   P  */  
      counts[16] = 32070;    /*   S  */  
      counts[17] = 26311;    /*   T  */  
      counts[18] =  5990;    /*   W  */  
      counts[19] = 14488;    /*   Y  */  
      counts[20] = 29012;    /*   V  */  
      counts[21] =     0;    /*   B  */  
      counts[22] =     0;    /*   Z  */  
      counts[23] =     0;    /*   X   */ 
      counts[24] =     0;    /*   *   */
    }
    else { /* SEQT_DNA */
      counts[1] = 250;
      counts[2] = 250;
      counts[3] = 250;
      counts[4] = 250;
      for (i=5; i<=ppst->nsq; i++) counts[i]=0;
    }
  } else {
    memset(&counts[0], 0, sizeof(counts));

    for(i = 0 ; i < n1 ; i++) {
      if(aa1[i] > ppst->nsq || aa1[i] < 1) continue;
      counts[aa1[i]]++;
    }
  }

  sum = 0;
  for(i = 1 ; i <= ppst->nsq ; i++) sum += counts[i];
  
  for(i = 1 ; i <= ppst->nsq ; i++) {
    if(n1 == 0) {
      priors[i] = (double) counts[i] / (double) sum;
    } else {
      priors[i] = ( ((double) pseudocts * f_str->priors[i]) + (double) counts[i] ) / ( (double) sum + (double) pseudocts );
    }
  }

  return;
} 

int
max_score(int *scores, int nsq) {

  int max, i;

  max = -BIGNUM;
  for ( i = 1 ; i <= nsq ; i++ ) {
    if (scores[i] > max) max = scores[i];
  }
 
  return max;
}

int
min_score(int *scores, int nsq) {

  int min, i;

  min = BIGNUM;
  for (i = 1 ; i <= nsq ; i++ ) {
    if (scores[i] < min) min = scores[i];
  }

  return min;
}

double
calc_tatusov ( struct slink *last,
	       struct slink *this,
	       const unsigned char *aa0, int n0,
	       const unsigned char *aa1, int n1,
	       int **pam2, int nsq,
	       struct f_struct *f_str,
	       int pseudocts,
	       int do_opt,
	       int zsflag
	       )
{
  int i, is, j, k;

  double *priors, my_priors[MAXSQ], tatprob, left_tatprob, right_tatprob;
  unsigned char *query = NULL;
  int length, maxlength, sumlength, sumscore, tmp, seg;
  int start, stop;
  struct slink *sl;
  int N;
  double *tatprobsptr;

#if defined(FASTS) || defined(FASTM)
  int index = 0;
  int notokay = 0;
#endif

  struct tat_str *oldtat = NULL, *newtat = NULL;

#if defined(FASTS) || defined(FASTM)
  start = this->vp->start - this->vp->dp + f_str->noff;
  stop = this->vp->stop - this->vp->dp + f_str->noff;
  tmp = stop - start + 1;
#else
  /*
    FASTF alignments can also hang off the end of library sequences,
    but no query residues are used up in the process, but we have to
    keep track of which are
  */
  tmp = 0;
  for(i = 0, j = 0 ; i < n0 ; i++) {
    if (this->vp->used[i] == 1) {tmp++; }
  }
#endif

  sumlength = maxlength = length = tmp;
  seg = 1;
  sumscore = this->vp->score;

#if defined(FASTS) || defined(FASTM)
  if(f_str->aa0b[start] == start && f_str->aa0e[stop] == stop) {
    index |= (1 << f_str->aa0i[start]);
  } else {
    notokay |= (1 << f_str->aa0i[start]);
  }
#endif

  for(sl = last; sl != NULL ; sl = sl->prev) {

#if defined(FASTS) || defined(FASTM)
    start = sl->vp->start - sl->vp->dp + f_str->noff;
    stop = sl->vp->stop - sl->vp->dp + f_str->noff;
    tmp = stop - start + 1;
#else
    tmp = 0;
    for(i = 0, j = 0 ; i < n0 ; i++) {
      if(sl->vp->used[i] == 1) {
	tmp++;
      }
    }
#endif
    sumlength += tmp;
    maxlength = tmp > maxlength ? tmp : maxlength;
    seg++;
    sumscore += sl->vp->score;

#if defined(FASTS) || defined(FASTM)
    if(f_str->aa0b[start] == start && f_str->aa0e[stop] == stop) {
      index |= (1 << f_str->aa0i[start]);
    } else {
      notokay |= (1 << f_str->aa0i[start]);
    }
#endif

  }

  tatprob = -1.0; 
    
#if defined(FASTS) || defined(FASTM)

  /* for T?FASTS, we try to use what we've precalculated: */

  /* with z = 3, do_opt is true, but we can use precalculated - with
     all other z's we can use precalculated only if !do_opt */
  if(!notokay && f_str->tatprobs != NULL) {
    /* create our own newtat and copy f_str's tat into it */
    index--;

    newtat = (struct tat_str *) malloc(sizeof(struct tat_str));
    if(newtat == NULL) {
      fprintf(stderr, "Couldn't calloc memory for newtat.\n");
      exit(1);
    }
    
    memcpy(newtat, f_str->tatprobs[index], sizeof(struct tat_str));

    newtat->probs = (double *) calloc(f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore + 1, sizeof(double));
    if(newtat->probs == NULL) {
      fprintf(stderr, "Coudln't calloc memory for newtat->probs.\n");
      exit(1);
    }

    memcpy(newtat->probs, f_str->tatprobs[index]->probs,
	   (f_str->tatprobs[index]->highscore - f_str->tatprobs[index]->lowscore + 1) * sizeof(double)); 


    tatprob = f_str->intprobs[index][sumscore - f_str->tatprobs[index]->lowscore];

  } else { /* we need to recalculate from scratch */
#endif

    /* for T?FASTF, we're always recalculating from scratch: */

    query = (unsigned char *) calloc(length, sizeof(unsigned char));
    if(query == NULL) {
      fprintf(stderr, "Couldn't calloc memory for query.\n");
      exit(1);
    }
    
#if defined(FASTS) || defined(FASTM)
    start = this->vp->start - this->vp->dp + f_str->noff;
    for(i = 0, j = 0 ; i < length ; i++) {
      query[j++] = aa0[start + i];
    }
#else
    for(i = 0, j = 0 ; i < n0 ; i++) {
      if (this->vp->used[i] == 1) {query[j++] = aa0[i];}
    }
#endif

    /*  calc_priors - not currently implemented for aa1 dependent */
    /* 
    if( (do_opt && zsflag == 2) || zsflag == 4 ) {    
      priors = &my_priors[0];
      calc_priors(priors, f_str, aa1, n1, pseudocts);
    } else {
      priors = f_str->priors;
    }
    */

    priors = f_str->priors;
    oldtat = (last != NULL ? last->tat : NULL);

    generate_tatprobs(query, 0, length - 1, priors, pam2, nsq, &newtat, oldtat);

    free(query);
#if defined(FASTS) || defined(FASTM)
  } /* close the FASTS-specific if-else from above */
#endif

  this->newtat = newtat;
  
  if(tatprob < 0.0) { /* hasn't been set by precalculated FASTS intprobs */

    /* integrate probabilities >= sumscore */
    tatprobsptr = newtat->probs;

    is = i = newtat->highscore - newtat->lowscore;
    N = sumscore - newtat->lowscore;

    right_tatprob = 0;
    for ( ;  i >= N; i--) {
      right_tatprob += tatprobsptr[i];
    }

    left_tatprob = tatprobsptr[0];
    for (i = 1 ; i < N ; i++ ) {
      left_tatprob += tatprobsptr[i];
    }

    if (right_tatprob < left_tatprob) {tatprob = right_tatprob;}
    else {tatprob = 1.0 - left_tatprob;}

    tatprob /= (right_tatprob+left_tatprob);
  }

  if (maxlength > 0) {
    n1 += 2 * (maxlength - 1);
  }

#ifndef FASTM
  tatprob *= factorial(n1 - sumlength + seg, n1 - sumlength);
#else
  tatprob *= power(n1 - sumlength,seg)/(1<<seg);
#endif

  if(tatprob > 0.01)
    tatprob = 1.0 - exp(-tatprob);
 
  return tatprob;
}

void
generate_tatprobs(const unsigned char *query,
		  int begin,
		  int end,
		  double *priors,
		  int **pam2,
		  int nsq,
		  struct tat_str **tatarg,
		  struct tat_str *oldtat)
{

  int i, j, k, l, m, n, N, highscore, lowscore;
  int *lowrange = NULL, *highrange = NULL;
  double *probs = NULL, *newprobs = NULL, *priorptr, tmp;
  struct tat_str *tatprobs = NULL;
  int *pamptr, *pamptrsave;

  if((tatprobs = (struct tat_str *) calloc(1, sizeof(struct tat_str)))==NULL) {
    fprintf(stderr, "Couldn't allocate individual tatprob struct.\n");
    exit(1);
  }

  n = end - begin + 1;

  if ( (lowrange = (int *) calloc(n, sizeof(int))) == NULL ) {
    fprintf(stderr, "Couldn't allocate memory for lowrange.\n");
    exit(1);
  }
  
  if ( (highrange = (int *) calloc(n, sizeof(int))) == NULL ) {
    fprintf(stderr, "Couldn't allocate memory for highrange.\n");
    exit(1);
  }

  /* calculate the absolute highest and lowest score possible for this */
  /* segment.  Also, set the range we need to iterate over at each position */
  /* in the query: */
  if(oldtat == NULL) {
    highscore = lowscore = 0;
  } else {
    highscore = oldtat->highscore;
    lowscore = oldtat->lowscore;
  }

  for ( i = 0 ; i < n ; i++ ) {

    if (query[begin+i] == 0) break;

    highscore =
      (highrange[i] = highscore + max_score(pam2[query[begin + i]], nsq));

    lowscore =
      (lowrange[i] = lowscore + min_score(pam2[query[begin + i]], nsq));

    /*
    fprintf(stderr, "i: %d, max: %d, min: %d, high[i]: %d, low[i]: %d, high: %d, low: %d, char: %d\n",
	    i,
	    max_score(pam2[query[begin + i]], nsq),
	    min_score(pam2[query[begin + i]], nsq),
	    highrange[i], lowrange[i],
	    highscore, lowscore, query[begin + i]); 
    */
  }

  /* allocate an array of probabilities for all possible scores */
  /* i.e. if highest score possible is 50 and lowest score possible */
  /* is -20, then there are 50 - (-20) + 1 = 71 possible different */
  /* scores (including 0): */
  N = highscore - lowscore;
  if ( (probs = (double *) calloc(N + 1, sizeof(double))) == NULL ) {
    fprintf(stderr, "Couldn't allocate probability matrix : %d.\n", N + 1);
    exit(1);
  }

  if(oldtat == NULL) {
    /* for the first position, iterate over the only possible scores, */
    /* summing the priors for the amino acids that can yield each score. */
    pamptr = pam2[query[begin]];
    for ( i = 1 ; i <= nsq ; i++ ) {
      if(priors[i] > 0.0) {
	probs[(pamptr[i] - lowscore)] += priors[i];
      }
    }
  } else {
    /* Need to copy the data out of oldtat->probs into probs */
    memcpy( &probs[oldtat->lowscore - lowscore],
	    oldtat->probs,
	    (oldtat->highscore - oldtat->lowscore + 1) * sizeof(double));
  }

  if ( (newprobs = (double *) calloc(N + 1, sizeof(double))) == NULL ) {
    fprintf(stderr, "Couldn't allocate newprobs matrix.\n");
    exit(1);
  }

  /* now for each remaining residue in the segment ... */
  for ( i = (oldtat == NULL ? 1 : 0) ; i < n ; i++ ) {

    pamptrsave = pam2[query[begin + i]];

    /* ... calculate new probability distribution .... */

    /* ... for each possible score (limited to current range) ... */
    for ( j = lowrange[i] - lowscore,
	    k = highrange[i] - lowscore ;
	  j <= k ;
	  j++ ) {
      
      tmp = 0.0;
      pamptr = &pamptrsave[1];
      priorptr = &priors[1];
      /* ... for each of the possible alignment scores at this position ... */
      for ( l = 1 ;
	    l <= nsq ;
	    l++) {

	/* make sure we don't go past highest possible score, or past
           the lowest possible score; not sure why this can happen */
	m = j - *pamptr++;
	if ( m <= N && m >= 0 ) {
	  /* update the probability of getting score j: */
	  tmp += probs[m] * *priorptr++;
	}
      }
      newprobs[j] += tmp;
    }

    /* save the new set of probabilities, get rid of old; we don't
       necessarily have to copy/clear all N+1 slots, we could use
       high/low score boundaries -- not sure that's worth the
       effort. */
    memcpy(probs, newprobs, (N + 1) * sizeof(double));
    memset(newprobs, 0, (N + 1) * sizeof(double));
  }

  free(newprobs);
  free(highrange);
  free(lowrange);

  tatprobs->probs = probs;
  /*  tatprobs->intprobs = intprobs; */
  tatprobs->lowscore = lowscore;
  tatprobs->highscore = highscore;

  *tatarg = tatprobs;
}

