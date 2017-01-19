/* map_db.c - read a FASTA or GCG format database and generate a list
   of indices for rapid memory mapping */

/* copyright (c) 1999 William R. Pearson */

/* $Name: fa_34_26_5 $ - $Id: map_db.c,v 1.9 2005/09/27 15:32:58 wrp Exp $ */

/* input is a libtype 1,5, or 6 sequence database */
/* output is a BLAST2 formatdb type index file */

/* format of the index file:

1)  map_db version number ["MP"+2 bytes]
2)  number of sequences in database [4 bytes]
3)  total length of database        [8 bytes]  (MP1, 4 bytes for MP0)
4)  longest sequence in database    [8 bytes]  (MP1, 4 bytes for MP0)
5) list of offsets to definitions  [num_seq+1] int*8 (MP1, 4 bytes for MP0)
6) list of offsets to sequences    [num_seq+1] int*8 (MP1, 4 bytes for MP1)
7) list of flag characters for sequences [num_seq+1]bytes
    (used for GCG binary to encode 2bit or 4 bit representation)

    sequence files will be as defined by their format
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "uascii.h"
#include "ncbl2_head.h"

#define GCGBIN 6
#define LASTLIB 6

int (*get_entry) ();

int a_get_ent(long *, long *);
int v_get_ent(long *, long *);
int gcg_get_ent(long *, long *);
int gbf_get_ent(long *, long *);

void src_int4_write(FILE *, int);
void src_int4_read(FILE *, int *);
void src_long4_write(FILE *, long);
void src_long4_read(FILE *, long *);
void src_long8_write(FILE *, long);
void src_long8_read(FILE *, long *);

void newname(char *nname, char *oname, char *suff, int maxn);

int (*get_ent_arr[LASTLIB+1])()={a_get_ent, gbf_get_ent, NULL, NULL, NULL,
				 v_get_ent, gcg_get_ent};

long openlib(char *, int);

static int *sascii;

main(int argc, char **argv)
{
  FILE *libi;
  char lname[256];
  char iname[256];
  char format[4];
  char *bp;

  int i;
  int nlib;	/* number of entries */

  long max_len;	/* longest sequence */
  long tot_len;	/* total sequence length */

  int n1;
  
  long f_size;	/* file size from fstat() */
  int lib_size;	/* current space available - may be realloc'ed */
  int lib_inc;
  int lib_type; /* 1 for protein, 0 for DNA */
  int lib_aa;	/* dna=1; prot=0; */

  /* file offsets */
  long d_pos;	/* start of description */
  long s_pos;	/* start of sequence */
  long *d_pos_arr;	/* array of description pointers */
  long *s_pos_arr;	/* array of description pointers */

  lib_type = 0;
  lib_size = 200000;
  lib_inc  = 100000;

  lib_aa = 1;

  while (argc > 1 && *argv[1]=='-') {
    if (strcmp(argv[1],"-n")==0) lib_aa = 0;
    argv++;
    argc--;
  }

  /* open the database */
  if (argc > 1) strncpy(lname, argv[1],sizeof(lname));
  else {
    fprintf(stderr," Entry library name: ");
    fgets(lname,sizeof(lname),stdin);
    if ((bp=strchr(lname,'\n'))!=NULL) *bp='\0';
  }
    
  if ((bp=strchr(lname,' '))!=NULL) {
    lib_type = atoi(bp+1);
    *bp='\0';
  }
  else lib_type = 0;

  if (get_ent_arr[lib_type] == NULL) {
    fprintf(stderr," cannot index file %s type %d\n",lname,lib_type);
    exit(1);
  }
  
  if (lib_type == 6) lib_aa = 0;
  if (lib_type == 1) lib_aa = 0;
  
  if (lib_aa == 1) sascii = aascii;
  else sascii = nascii;

  if ((f_size=openlib(lname,lib_type))==0) {
    fprintf(stderr," cannot open %s (type: %d)\n",lname,lib_type);
    exit(1);
  }

  /* allocate array of description pointers */
  if ((d_pos_arr=(long *)calloc(lib_size, sizeof(long)))==NULL) {
    fprintf(stderr," cannot allocate %d for desc. array\n",lib_size);
    exit(1);
  }
  /* allocate array of sequence pointers */
  if ((s_pos_arr=(long *)calloc(lib_size, sizeof(long)))==NULL) {
    fprintf(stderr," cannot allocate %d for seq. array\n",lib_size);
    exit(1);
  }

  /* allocate array of sequence flags */

  nlib = 0; tot_len=0; max_len=-1;
  while ((n1=get_entry(&d_pos, &s_pos)) > 0) {
    d_pos_arr[nlib] = d_pos;
    s_pos_arr[nlib] = s_pos;
    nlib++;
    tot_len += n1;
    if (n1 > max_len) max_len = n1;
    if (nlib >= lib_size) { /* too many entries */
      lib_size += lib_inc;
      if ((d_pos_arr=(long *)realloc(d_pos_arr,lib_size*sizeof(long)))==NULL) {
	fprintf(stderr," cannot realloc allocate %d for desc.. array\n",
		lib_size);
	exit(1);
      }
      if ((s_pos_arr=(long *)realloc(s_pos_arr,lib_size*sizeof(long)))==NULL) {
	fprintf(stderr," cannot realloc allocate %d for seq. array\n",
		lib_size);
	exit(1);
      }
    }
  }

  d_pos_arr[nlib]= d_pos;	/* put in the end of the file */
  s_pos_arr[nlib]=0;

  /* all the information is in, write it out */
  
  newname(iname,lname,"xin",sizeof(iname));

  if ((libi=fopen(iname,"w"))==NULL) {
    fprintf(stderr," cannot open %s for writing\n",iname);
    exit(1);
  }

  /* write out format version */
  format[0]='M';
  format[1]='P';
#ifdef BIG_LIB64
  format[2]= 1;		/* format 1 for 8-byte offsets */
#else
  format[2]='\0';	/* format '\0' for original 4-byte */
#endif

  format[3]=lib_type;
  fwrite(format,4,sizeof(char),libi);

  /* write out sequence type */
  src_int4_write(libi, lib_aa);

  /* write out file fstat as integrity check */
#ifdef BIG_LIB64
  src_long8_write(libi, f_size);
#else
  src_int4_write(libi, f_size);
#endif

  /* write out num_seq */
  src_int4_write(libi, nlib);

#ifdef BIG_LIB64
  /* write out tot_len, max_len */
  src_long8_write(libi, tot_len);
#else
  src_int4_write(libi, tot_len);
#endif
  src_int4_write(libi, max_len);

#ifdef BIG_LIB64
  for (i=0; i<=nlib; i++) src_long8_write(libi,d_pos_arr[i]);
  for (i=0; i<=nlib; i++) src_long8_write(libi,s_pos_arr[i]);
#else
  for (i=0; i<=nlib; i++) src_int4_write(libi,d_pos_arr[i]);
  for (i=0; i<=nlib; i++) src_int4_write(libi,s_pos_arr[i]);
#endif

  fclose(libi);

#ifdef BIG_LIB64
  fprintf(stderr," wrote %d sequences (tot=%ld, max=%ld) to %s\n",
	  nlib,tot_len,max_len,iname);
#else
  fprintf(stderr," wrote %d sequences (tot=%ld, max=%ld) to %s\n",
	  nlib,tot_len,max_len,iname);
#endif
}


FILE *libf=NULL;
long lpos;

#define MAXLINE 4096
char lline[MAXLINE+1];

long
openlib(char *lname, int lib_type)
{
  long f_size;
  struct stat stat_buf;

  if (stat(lname,&stat_buf)<0) {
    fprintf(stderr," cannot stat library: %s\n",lname);
    return 0;
  }

  if ((libf=fopen(lname,"r"))==NULL) {
    fprintf(stderr," cannot open library: %s (type: %d)\n",
	    lname, lib_type);
    return 0;
  }
  
  f_size = stat_buf.st_size;

  get_entry = get_ent_arr[lib_type];

  lpos = ftell(libf);
  if (fgets(lline,MAXLINE,libf)==NULL) return 0;
  return f_size;
}

int
a_get_ent(long *d_pos, long *s_pos)
{
  register char *cp;
  register int *ap, n1;

  ap = sascii;

  while (lline[0]!='>' && lline[0]!=';') {
    lpos = ftell(libf);
    if (fgets(lline,sizeof(lline),libf)==NULL) {
      *d_pos = lpos;
      return 0;
    }
  }

  *d_pos = lpos;

  /* make certain we have the end of the line */
  while (strchr((char *)lline,'\n')==NULL) {
    if (fgets(lline,sizeof(lline),libf)==NULL) break;
  }

  *s_pos = ftell(libf);
  lline[0]='\0';
  n1 = 0;
  while (fgets(lline,sizeof(lline),libf)!=NULL) {
    if (lline[0]=='>') break;
    if (lline[0]==';') {
      if (strchr(lline,'\n')==NULL) {
	fprintf(stderr," excessive continuation\n%s",lline);
	return -1;
      }
    }

    for (cp=lline; *cp; ) if (ap[*cp++]<NA) n1++;
    lpos = ftell(libf);
  }
  return n1;
}

int
v_get_ent(long *d_pos, long *s_pos)
{
  register char *cp;
  register int *ap;
  int n1;

  ap = sascii;

  /* check for seq_id line */
  while (lline[0]!='>' && lline[0]!=';') {
    lpos = ftell(libf);
    if (fgets(lline,sizeof(lline),libf)==NULL) {
      *d_pos = lpos;
      return 0;
    }
  }
  *d_pos = lpos;

  /* get the description line */
  if (fgets(lline,sizeof(lline),libf)==NULL) return 0;
  /* make certain we have the end of the line */
  while (strchr((char *)lline,'\n')==NULL) {
    if (fgets(lline,sizeof(lline),libf)==NULL) break;
  }

  *s_pos = ftell(libf);
  lline[0]='\0';
  n1 = 0;
  while (fgets(lline,sizeof(lline),libf)!=NULL) {
    if (lline[0]=='>') break;

    for (cp=lline; *cp; ) if (ap[*cp++]<NA) n1++;
    lpos = ftell(libf);
  }
  return n1;
}

static char gcg_type[10];
static long gcg_len;
static int gcg_bton[4]={2,4,1,3};

int
gcg_get_ent(long *d_pos, long *s_pos)
{
  register char *cp;
  register int *ap;
  char libstr[20], dummy[20];
  char gcg_date[6];
  int r_block;
  int n1;

  /* check for seq_id line */
  while (lline[0]!='>') {
    lpos = ftell(libf);
    if (fgets(lline,sizeof(lline),libf)==NULL) {
      *d_pos = lpos;
      return 0;
    }
  }
  *d_pos = lpos;

  /* get the encoding/sequence length info */

  sscanf(&lline[4],"%s %s %s %s %ld",
	 libstr,gcg_date,gcg_type,dummy,&gcg_len);

  /* get the description line */
  if (fgets(lline,MAXLINE,libf)==NULL) return;

  *s_pos = ftell(libf);
  /* seek to the end of the sequence; +1 to jump over newline */
  if (gcg_type[0]=='2') {
    r_block = (gcg_len+3)/4;
    fseek(libf,r_block+1,SEEK_CUR);
  }
  else fseek(libf,gcg_len+1,SEEK_CUR);

  lpos = ftell(libf);
  fgets(lline,MAXLINE,libf);

  return gcg_len;
}

int
gbf_get_ent(long *d_pos, long *s_pos)
{
  int n1;
  char *cp;
  register int *ap;

#if !defined(TFAST)
  ap = sascii;
#else
  ap = nascii;
#endif

  while (lline[0]!='L' || lline[1]!='O' || 
	 strncmp(lline,"LOCUS",5)) { /* find LOCUS */
    lpos = ftell(libf);
    if (fgets(lline,MAXLINE,libf)==NULL) return (-1);
  }
  *d_pos=lpos;

  while (lline[0]!='O' || lline[1]!='R' ||
	 strncmp(lline,"ORIGIN",6)) { /* find ORIGIN */
    if (fgets(lline,MAXLINE,libf)==NULL) return (-1);
  }
  *s_pos = ftell(libf);

  lline[0]='\0';
  n1=0;
  while (fgets(lline,MAXLINE,libf)!=NULL) {
    if (lline[0]=='/') break;
    for (cp=lline; *cp; ) if (ap[*cp++]<NA) n1++;
  }
  lpos = ftell(libf);
  fgets(lline,MAXLINE,libf);

  return n1;
}

void src_int4_read(FILE *fd,  int *val)
{
#ifdef IS_BIG_ENDIAN
  fread((char *)val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
#endif
}

void src_int4_write(FILE *fd,  int val)
{
#ifdef IS_BIG_ENDIAN
  fwrite(&val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  b[3] = val & 255;
  b[2] = (val=val>>8)&255;
  b[1] = (val=val>>8)&255;
  b[0] = (val=val>>8)&255;

  fwrite(b,(size_t)1,(size_t)4,fd);
#endif
}

void src_long8_write(FILE *fd,  long val)
{
#ifdef IS_BIG_ENDIAN
  fwrite(&val,(size_t)8,(size_t)1,fd);
#else
  unsigned char b[8];

  b[7] = val & 255;
  b[6] = (val=val>>8)&255;
  b[5] = (val=val>>8)&255;
  b[4] = (val=val>>8)&255;
  b[3] = (val=val>>8)&255;
  b[2] = (val=val>>8)&255;
  b[1] = (val=val>>8)&255;
  b[0] = (val=val>>8)&255;

  fwrite(b,(size_t)1,(size_t)8,fd);
#endif
}

void
newname(char *nname, char *oname, char *suff, int maxn)
{
  strncpy(nname,oname,maxn-1);
  strncat(nname,".",1);
  strncat(nname,suff,maxn-strlen(nname));
}
