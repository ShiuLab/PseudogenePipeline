
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: nrand48.c,v 1.4 2006/04/12 18:00:02 wrp Exp $ */

#include <stdlib.h>
#include <time.h>

void 
irand(int n)	/* initialize random number generator */
{
  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;
  }
  srand48(n);
}

int
nrand(int n)	/* returns a random number between 0 and n-1 where n < 64K) */
{
  int rn;

  rn = lrand48();
  rn = rn >> 16;
  rn = (rn % n);
  return rn;
}

