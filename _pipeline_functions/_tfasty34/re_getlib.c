/* re_getlib.c - re-acquire a sequence given lseek, lcont */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "mm_file.h"

#define GETLIB (m_fptr->getlib)

int
re_getlib(unsigned char *aa1,
	  int maxn,	/* longest aa1 */
	  int maxt3,	/* alternate maxn */
	  int loff,	/* overlap */
	  int lcont,
	  int term_code,
	  long *loffset,	/* offset from real start of sequence */
	  long *l_off_p,	/* coordinate of sequence start */
	  struct lmf_str *m_fptr) {

  unsigned char *aa1ptr;
  int icont, maxt, ccont, n1;
  char libstr[20];
  fseek_t lmark; 
  
  aa1ptr = aa1;
  icont=0;

  *loffset = 0l;
  maxt = maxn;
  n1 = -1;
  for (ccont=0; ccont<=lcont-1; ccont++) {

    n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);

    if (term_code && m_fptr->lib_aa && aa1ptr[n1-1]!=term_code) {
      aa1ptr[n1++]=term_code;
      aa1ptr[n1]=0;
    }

    if (aa1ptr!=aa1) n1 += loff;

    if (icont>lcont-1) break;

    if (icont) {
      maxt = maxt3;
      memcpy(aa1,&aa1[n1-loff],loff);
      aa1ptr= &aa1[loff];
      *loffset += n1 - loff;
    }
    else {
      maxt = maxn;
      aa1ptr=aa1;
    }
  }
  return n1;
}
