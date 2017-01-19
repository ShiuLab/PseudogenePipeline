
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: nrand.c,v 1.2 2005/09/23 16:27:25 wrp Exp $ */

#include <stdlib.h>
#include <time.h>

int
irand(int n)	/* initialize random number generator */
{

  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;

  }
  srand(n);
}

int
nrand(int n)	/* returns a random number between 1 and n where n < 64K) */
{
  int rand();
  long rn;

  rn = rand();
#ifdef RAND32
  rn = rn >> 16;
#endif
  rn = rn % n;
  return (int)rn;
}




