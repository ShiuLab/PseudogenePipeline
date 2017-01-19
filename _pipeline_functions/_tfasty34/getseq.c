/*	May, June 1987	- modified for rapid read of database

	copyright (c) 1987,1988,1989,1992,1995,2000 William R. Pearson

	This is one of three alternative files that can be used to
	read a database.  The three files are nxgetaa.c, nmgetaa.c, and
	mmgetaa.c.

	nxgetaa.c contains the original code for reading databases, and
	is still used for Mac and PC versions of fasta33 (which do not
	use mmap).

	nmgetaa.c and mmgetaa.c are used together.  nmgetaa.c provides
	the same functions as nxgetaa.c if memory mapping is not used,
	mmgetaa.c provides the database reading functions if memory
	mapping is used. The decision to use memory mapping is made on
	a file-by-file basis.

	June 2, 1987 - added TFASTA
	March 30, 1988 - combined ffgetaa, fgetgb;
	April 8, 1988 - added PIRLIB format for unix
	Feb 4, 1989 - added universal subroutines for libraries
	December, 1995 - added range option file.name:1-1000
	Feb 22, 2002 - fix to allow "plain" text file queries

	getnt.c	associated subroutines for matching sequences */

/* $Name: fa_34_26_5 $ - $Id: getseq.c,v 1.13 2006/10/05 18:22:07 wrp Exp $ */

/*
	8-April-88
	The compile time #define PIRLIB allows this routine to be used
	to read protein and DNA sequence libraries in the NBRF/PIR
	VAX/VMS library format.  That is:

	>P1;LCBO
	This is a line of description
	GTYH ... the sequence starts on this line

	This may ease conversion from UWGCG format libraries. It
	has not been extensively tested.

	In addition, sequence libraries with a '>' in the 4th position
	are recognized as NBRF format libraries for consistency with
	UWGCG
*/

/* 	Nov 12, 1987	- this version checks to see if the sequence
	is DNA or protein by asking whether > 85% is A, C, G, T

	May 5, 1988 - modify the DNA/PROTEIN checker by re-reading
	DNA sequences in order to check for 'U'.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"

#ifndef SFCHAR
#define SFCHAR ':'
#endif

#ifdef VMS
#define PIRLIB
#endif

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

#define YES 1
#define NO 0
#define MAXLINE 512

#ifndef min
#define min(x,y) ((x) > (y) ? (y) : (x))
#endif

#ifdef SUPERFAMNUM
extern int nsfnum;	/* number of superfamily numbers */
extern int sfnum[];	/* superfamily number from types 0 and 5 */
extern int nsfnum_n;
extern int sfnum_n[];
#endif

#define NO_FORMAT 0
#define FASTA_FORMAT 1
#define GCG_FORMAT 2

static int seq_format=NO_FORMAT;
static char seq_title[200];

int scanseq(unsigned char *, int, char *);
void sf_sort(int *, int);
extern void init_ascii(int is_ext, int *sascii, int is_dna);

/* getseq	- get a query sequence, possibly re-reading to set type
   returns	- length of query sequence or error = 0

   char *filen	- name of file to be opened
   char *seq	- destination for query sequence
   int maxs	- maximum length of query
   char libstr[20] - short description (locus or acc)
   int *dnaseq  - -1 => use scanseq to determine sequence type
   		   0 => must be protein
		   1 => must be DNA
   long *sq0off	- offset into query specified by query_file:1001-2000
*/   

int
getseq(char *filen, int *qascii, unsigned char *seq, int maxs, char *libstr, long *sq0off)
{
  FILE *fptr;
  char line[512],*bp, *bp1, *bpn, *tp;
  int i, rn, n;
  int ic;
  int sstart, sstop, sset=0;
  int llen, l_offset;
#ifdef SUPERFAMNUM
  static char tline[MAXLINE];
#endif

  seq_title[0]='\0';
  libstr[0]='\0';

  sstart = sstop = -1;
#ifndef DOS
  if ((bp=strchr(filen,':'))!=NULL && *(bp+1)!='\0') {
#else
  if ((bp=strchr(filen+3,':'))!=NULL && *(bp+1)!='\0') {
#endif
    *bp='\0';
    if (*(bp+1)=='-') {
      sstart = 0;
      sscanf(bp+2,"%d",&sstop);
    }
    else {
      sscanf(bp+1,"%d-%d",&sstart,&sstop);
      sstart--;
      if (sstop <= 0 ) sstop = BIGNUM;
    }
    sset=1;
  }
  else {
    sstart = 0;
    sstop = BIGNUM;
  }

  /* check for input from stdin */
  if (strcmp(filen,"-") && strcmp(filen,"@")) {
    if ((fptr=fopen(filen,"r"))==NULL) {
      fprintf(stderr," could not open %s\n",filen);
      return 0;
    }
  }
  else {
    fptr = stdin;
  }
  rn = n=0;

  while(fgets(line,sizeof(line),fptr)!=NULL) {
#ifdef PIRLIB
    if (line[0]=='>'&& (line[3]==';'||line[3]=='>'))
      fgets(line,sizeof(line),fptr);
    else
#endif
    l_offset = 0;
    if (line[0]=='>') {
      seq_format = FASTA_FORMAT;
#ifdef SUPERFAMNUM
      sfnum[nsfnum=0]= sfnum_n[nsfnum_n=0]=0;
      strncpy(tline,line+1,sizeof(tline));
      tline[sizeof(tline)-1]='\0';

      if ((bp=strchr(tline,' ')) && (bp=strchr(bp+1,SFCHAR))) {
	if ((bp1=strchr(bp+1,SFCHAR))==NULL) {
	  fprintf(stderr," second %c missing: %s\n",SFCHAR,tline);
	}
	else {
	  if ((bpn=strchr(bp+1,NSFCHAR))!=NULL) *bpn = '\0';
	  *bp1 = '\0';
	  i = 0;
	  if ((tp = strtok(bp+1," \t"))!=NULL)  {
	    sfnum[i++] = atoi(tp);
	    while ((tp = strtok((char *)NULL," \t")) != (char *)NULL) {
	      if (isdigit(*tp)) sfnum[i++] = atoi(tp);
	      if (i>=9) break;
	    }
	  }
	  sfnum[nsfnum=i]= 0;
	  if (nsfnum>1) sf_sort(sfnum,nsfnum);
	  else {
	    if (nsfnum < 1) fprintf(stderr," found | but no sfnum: %s\n",libstr);
	  }
	  if (bpn != NULL) {
	    tp = strtok(bpn+1," \t");
	    sfnum_n[0]=atoi(tp);
	    i = 1;
	    while ((tp=strtok(NULL," \t"))!=NULL) {
	      sfnum_n[i++] = atoi(tp);
	      if (i >= 10) {
		fprintf(stderr,
			" error - too many negative superfamilies: %d\n %s\n",
			i,tline);
		break;
	      }
	    }
	    sfnum_n[nsfnum_n=i]=0;
	    sf_sort(sfnum_n,nsfnum_n);
	  }
	}
      }
      else {
	sfnum[nsfnum = 0] = 0;
	sfnum_n[nsfnum_n = 0] = 0;
      }
#endif
      if ((bp=(char *)strchr(line,'\n'))!=NULL) *bp='\0';
      strncpy(seq_title,line+1,sizeof(seq_title));
      seq_title[sizeof(seq_title)-1]='\0';
      if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
      strncpy(libstr,line+1,12);
      libstr[12]='\0';
    }
    else if (seq_format==NO_FORMAT && strcmp(line,"..")==0) {
      seq_format = GCG_FORMAT;
/*
      if (*dnaseq != 1) qascii['*'] = qascii['X'];
*/
      l_offset = 10;
      llen = strlen(line);
      while (strncmp(&line[llen-3],"..\n",(size_t)3) != 0) {
	if (fgets(line,sizeof(line),fptr)==NULL) return 0;
	llen = strlen(line);
      }
      bp = strtok(line," \t");
/*
      if ((bp=(char *)strchr(line,' '))!=NULL) *bp='\0';
      else if ((bp=(char *)strchr(line,'\n'))!=NULL) *bp='\0';
*/
      if (bp!=NULL) strncpy(libstr,bp,12);
      else strncpy(libstr,filen,12);
      libstr[12]='\0';
      if (fgets(line,sizeof(line),fptr)==NULL) return 0;
    }
    else {
      if (libstr[0]=='\0') strncpy(libstr,filen,12);
      libstr[12]='\0';
    }

    if (seq_format==GCG_FORMAT && strlen(line)<l_offset) continue;

    if (line[0]!='>'&& line[0]!=';') {
      for (i=l_offset; (n<maxs && rn < sstop)&&
	     ((ic=qascii[line[i]&AAMASK])<EL); i++)
	if (ic<NA && ++rn > sstart) seq[n++]= ic;
      if (ic == ES || rn > sstop) break;
    }
  }

  if (n==maxs) {
    fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
    fflush(stderr);
  }
  if ((bp=strchr(libstr,'\n'))!=NULL) *bp = '\0';
  if ((bp=strchr(libstr,'\r'))!=NULL) *bp = '\0';
  seq[n]= EOSEQ;


  if (seq_format !=GCG_FORMAT) 
    while(fgets(line,sizeof(line),fptr)!=NULL) {
#ifdef PIRLIB
      if (line[0]=='>'&& (line[3]==';'||line[3]=='>'))
	fgets(line,sizeof(line),fptr);
      else
#endif
	if (line[0]!='>'&& line[0]!=';') {
	  for (i=0; (n<maxs && rn < sstop)&&
		 ((ic=qascii[line[i]&AAMASK])<EL); i++)
	    if (ic<NA && ++rn > sstart ) seq[n++]= ic;
	  if (ic == ES || rn > sstop) break;
	}
    }
  else {
    llen = strlen(line);
    while (strncmp(&line[llen-3],"..\n",(size_t)3) != 0) {
      if (fgets(line,sizeof(line),fptr)==NULL) return 0;
      llen = strlen(line);
    }
    while (fgets(line,sizeof(line),fptr)!=NULL) {
      if (strlen(line)<l_offset) continue;
      for (i=l_offset; (n<maxs && rn < sstop) &&
	     ((ic=qascii[line[i]&AAMASK])<EL); i++)
	if (ic<NA && ++rn > sstart ) seq[n++]= ic;
      if (ic == ES || rn > sstop ) break;
    }
  }

  if (n==maxs) {
    fprintf(stderr," sequence may be truncated %d %d\n",n,maxs);
    fflush(stderr);
  }
  seq[n]= EOSEQ;

  if (fptr!=stdin) fclose(fptr);

  if (sset==1) {
    sstart++;
    filen[strlen(filen)]=':';
    if (*sq0off==1 || sstart>=1) *sq0off = sstart;
  }

  return n;
}

int
gettitle(char *filen, char *title, int len) {
  FILE *fptr;
  char line[512];
  char *bp;
  int sset;
#ifdef WIN32
  char *strpbrk();
#endif

  sset = 0;

  if (strncmp(filen,"-",1)==0 || strncmp(filen,"@",1)==0) {
    strncpy(title,seq_title,len);
    title[len-1]='\0';
    return (int)strlen(title);
  }

  if ((bp=strchr(filen,':'))!=NULL) { *bp='\0'; sset=1;}
	  

  if ((fptr=fopen(filen,"r"))==NULL) {
    fprintf(stderr," file %s was not found\n",filen);
    fflush(stderr);
    return 0;
  }

  if (sset==1) filen[strlen(filen)]=':';

  while(fgets(line,sizeof(line),fptr)!=NULL) {
    if (line[0]=='>'|| line[0]==';') goto found;
  }
  fclose(fptr);
  title[0]='\0';
  return 0;

 found:

#ifdef PIRLIB
  if (line[0]=='>'&&(line[3]==';'||line[3]=='>')) {
    if ((bp = strchr(line,'\n'))!=NULL) *bp='\0';
    ll=strlen(line); line[ll++]=' '; line[ll]='\0';
    fgets(&line[ll],sizeof(line)-ll,fptr);
  }
#endif

#ifdef WIN32
  bp = strpbrk(line,"\n\r");
#else
  bp = strchr(line,'\n');
#endif
  if (bp!=NULL) *bp = 0;
  strncpy(title,line,len);
  title[len-1]='\0';
  fclose(fptr);
  return strlen(title);
}

