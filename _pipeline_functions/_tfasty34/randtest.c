#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

main(argc, argv)
     int argc; char **argv;
{
  int i, n, s;
  struct timeval t;

  if (argc < 2) n = 10;
  else n = atoi(argv[1]);

  gettimeofday(&t,NULL);
  printf(" seed: %d\n",t.tv_usec);
  srandom(t.tv_usec);

  for (i=0; i< n; i++)
    printf("%3d\n",random()%100);

}
