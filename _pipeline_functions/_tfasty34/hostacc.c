
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: hostacc.c,v 1.7 2006/04/12 18:00:02 wrp Exp $ */

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef PVM_SRC
#include "pvm3.h"
#endif
#ifdef MPI_SRC
#include "mpi.h"
#endif

#include "msg.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

extern char prog_name[];

extern int nnodes;
#ifdef PVM_SRC
extern int pinums[];
#endif


#ifdef PVM_SRC
int tidtonode(tid)
     int tid;
{
  int i;
  for (i=FIRSTNODE; i< nnodes; i++) if (tid==pinums[i]) return i;
  return -1;
}
#endif

/* rand_nodes selects nnodes at random from max_nodes */

void
rand_nodes(int *node_map, int nnodes, int max_nodes)
{
  int node_used[MAXNOD];
  int i, j;
  struct timeval tv;

  gettimeofday(&tv,NULL);
  SRAND(tv.tv_usec);

  for (i=0; i<max_nodes; i++) node_used[i]=0;

  if (nnodes < (max_nodes+1)/2) {
    for (i=0; i<nnodes; ) {
      j = RAND()%max_nodes;
      if (node_used[j]) continue;
      else {
	node_map[i++]=j;
	node_used[j]=1;
      }
    }
  }
  else {
    for (i=0; i<(max_nodes-nnodes); ) {
      j = RAND()%max_nodes;
      if (node_used[j]) continue;
      else {
	node_used[j]=1;
	i++;
      }
    }
    for (i=j=0; i<nnodes; j++)
      if (node_used[j]) continue;
      else node_map[i++]=j;
  }
/*  for (i=0; i<nnodes; i++) fprintf(stderr,"%2d %2d\n",i,node_map[i]); */
}
