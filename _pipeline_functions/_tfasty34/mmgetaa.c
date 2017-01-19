/* mmgetaa.c - functions for mmap()ed access to libraries */

/* copyright (c) 1999,2000 William R. Pearson */

/* version 0 September, 1999 */

/*
  This is one of two alternative files that can be used to
  read a database.  The two files are nmgetaa.c, and mmgetaa.c
  (nxgetaa.c has been retired).

  nmgetlib.c and mmgetaa.c are used together. nmgetlib.c provides
  the same functions as nxgetaa.c if memory mapping is not used,
  mmgetaa.c provides the database reading functions if memory
  mapping is used. The decision to use memory mapping is made on
  a file-by-file basis.
*/

/* $Name: fa_34_26_5 $ - $Id: mmgetaa.c,v 1.41 2006/04/12 18:00:02 wrp Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#define MAXLINE 512
#define EOSEQ 0

#define XTERNAL
#include "uascii.h"
/* #include "upam.h" */
#undef XTERNAL

#ifdef SUPERFAMNUM
extern int nsfnum;	/* number of superfamily numbers */
extern int sfnum[10];	/* superfamily number from types 0 and 5 */
extern int nsfnum_n;
extern int sfnum_n[10];
static char tline[MAXLINE];
#endif

#define GCGBIN 6

#ifndef MAP_FILE
#define MAP_FILE 0
#endif

#include "defs.h"
#include "mm_file.h"

extern MM_OFF bl2_long8_cvt(int64_t);
extern int bl2_uint4_cvt(int);


long crck(char *, int);
extern void src_int4_read(FILE *fd,  int *val);
extern void src_long4_read(FILE *fd,  long  *valp);
extern void src_long8_read(FILE *fd,  int64_t *val);

/* load_mmap() loads the d_pos[] and s_pos[] arrays for rapid access */

struct lmf_str *
load_mmap(FILE *libi,	/* fd for already open ".xin" file */
	  char *sname,	/* name of sequence database file */
	  int lib_type,	/* 0-Fasta, 5-vms_pir, 6-gcg_binary */
	  int ldnaseq,	/* 1 for DNA, 0 for protein */
	  struct lmf_str *m_fd)
{
  char format[4];
  int i, lib_aa;
  MM_OFF f_size;
  long lf_size;
  struct stat statbuf;
  int max_cnt;
  MM_OFF *d_pos_arr, *s_pos_arr;
  int mm_flag, mm64_flag;
  int *tmp_pos_arr;

  /* first check that the necessary indices are up-to-date */
  /* read the offsets in ".xin" file */
  if (fread(format,1,4,libi)==0) {
    fprintf(stderr," cannot read .xin format\n");
    return NULL;
  }
    
  mm64_flag = (format[2]==1);	/* 4 bytes or 8 bytes for long? */

#ifndef BIG_LIB64
  if (mm64_flag) {return NULL;}
#endif

  if (format[3]!=lib_type) {
    fprintf(stderr," cannot read format %d != lib_type %d\n",
	    format[3],lib_type);
    return NULL;
  }

  src_int4_read(libi,&lib_aa);
  if (lib_aa == ldnaseq) { /* database residue mismatch */
    fprintf(stderr," residue type mismatch %s != %s (.xin) in %s\n",
	    (lib_aa ? "DNA" : "prot."),(ldnaseq ? "prot." : "DNA"),
	    sname);
    return NULL;
  }
    
  /* everything looks good, allocate an lmf_str */

  m_fd->lib_aa = lib_aa;

  /* get get file size from index */
  if (mm64_flag) src_long8_read(libi,&f_size);
  else {
    src_long4_read(libi,&lf_size);
    f_size = lf_size;
  }

  /* now, start to open mmap()ed file */
  mm_flag=((m_fd->mmap_fd=open(sname,O_RDONLY))>=0);
  if (!mm_flag) {
    fprintf(stderr," cannot open %s for mmap()", sname);
    perror("...");
    return NULL;	/* file did not open */
  }

  /* fstat the library file and get size */
  if(fstat(m_fd->mmap_fd, &statbuf) < 0) {
    fprintf(stderr," cannot stat %s for mmap()", sname);
    perror("...");
    m_fd->mm_flg = 0;
    goto finish;
  }

  /* check for identical sizes - if different, do not mmap */
  if (f_size != statbuf.st_size) {
    fprintf(stderr," %s file size (%lld) and expected size (%ld) don't match\n",
	    sname,statbuf.st_size,f_size);
    mm_flag = 0;
    goto finish;    
  }

  /* the index file and library file are open and the sizes match */
  /* allocate the m_file struct and  map the file */

  m_fd->st_size = statbuf.st_size;
  if((m_fd->mmap_base = 
      mmap(NULL, m_fd->st_size, PROT_READ,
	   MAP_FILE | MAP_SHARED, m_fd->mmap_fd, 0)) == (char *) -1) {
    mm_flag = 0;
#ifdef DEBUG
    fprintf(stderr," cannot mmap %s", sname);
    perror("...");
#endif
  }  
 finish:
  close(m_fd->mmap_fd);
  if (!mm_flag) { return NULL; }

  /* now finish reading the index file */
  src_int4_read(libi,&max_cnt);

  if (mm64_flag) {
    src_long8_read(libi,&m_fd->tot_len);
  }
  else {
    src_long4_read(libi,&lf_size);
    m_fd->tot_len = lf_size;
  }
  src_long4_read(libi,&lf_size);
  m_fd->max_len = lf_size;

#ifdef DEBUG
  fprintf(stderr,
	  "%s\tformat: %c%c%d %d; max_cnt: %d; tot_len: %lld max_len: %ld\n",
	  sname,format[0],format[1],format[2],format[3],
	  max_cnt,m_fd->tot_len,m_fd->max_len);
#endif

  /* allocate array of description pointers */
  if (!mm64_flag) {
    if ((tmp_pos_arr=(int *)calloc(max_cnt+1,sizeof(int)))==NULL) {
      fprintf(stderr," cannot allocate %d for tmp_pos array\n",
	      max_cnt+1);
    }
  }

  if ((d_pos_arr=(MM_OFF *)calloc(max_cnt+1, sizeof(MM_OFF)))==NULL) {
    fprintf(stderr," cannot allocate %d for desc. array\n",max_cnt+1);
    exit(1);
  }

  /* read m_fd->d_pos[max_cnt+1] */
  if (mm64_flag) {
    if (fread(d_pos_arr,sizeof(MM_OFF),max_cnt+1,libi)!=
	max_cnt+1) {
      fprintf(stderr," error reading desc. offsets: %s\n",sname);
      return NULL;
    }
  }
  else {
    if (fread(tmp_pos_arr,sizeof(int),max_cnt+1,libi)!=
	max_cnt+1) {
      fprintf(stderr," error reading desc. offsets: %s\n",sname);
      return NULL;
    }
#ifdef DEBUG
    fprintf(stderr,"d_pos_crc: %ld\n",
	    crck((char *)tmp_pos_arr,sizeof(int)*(max_cnt+1)));
#endif
  }


#ifndef IS_BIG_ENDIAN
  if (mm64_flag)
    for (i=0; i<=max_cnt; i++) {
      d_pos_arr[i] = bl2_long8_cvt(d_pos_arr[i]);
    }
  else
    for (i=0; i<=max_cnt; i++) {
      d_pos_arr[i] = bl2_uint4_cvt(tmp_pos_arr[i]);
    }
#else
  if (!mm64_flag) {
    for (i=0; i<=max_cnt; i++) {
      d_pos_arr[i] = tmp_pos_arr[i];
    }
  }
#endif

#ifdef DEBUG
  for (i=0; i<max_cnt-1; i++) {
    if (d_pos_arr[i+1] <= d_pos_arr[i] )
      fprintf(stderr," ** dpos_error [%d]\t%ld\t%ld\n",
	      i,d_pos_arr[i],d_pos_arr[i+1]);
  }
#endif

  /* allocate array of sequence pointers */
  if ((s_pos_arr=(MM_OFF *)calloc(max_cnt+1,sizeof(MM_OFF)))==NULL) {
    fprintf(stderr," cannot allocate %d for seq. array\n",max_cnt+1);
    exit(1);
  }

  /* read m_fd->s_pos[max_cnt+1] */
  if (mm64_flag) {
    if (fread(s_pos_arr,sizeof(long),max_cnt+1,libi)!=
	max_cnt+1) {
      fprintf(stderr," error reading seq offsets: %s\n",sname);
      return NULL;
    }
  }
  else {
    if (fread(tmp_pos_arr,sizeof(int),max_cnt+1,libi)!=
	max_cnt+1) {
      fprintf(stderr," error reading seq offsets: %s\n",sname);
      return NULL;
    }
#ifdef DEBUG
    fprintf(stderr,"s_pos_crc: %ld\n",
	    crck((char *)tmp_pos_arr,sizeof(int)*(max_cnt+1)));
#endif
  }

#ifndef IS_BIG_ENDIAN
  if (mm64_flag)
    for (i=0; i<=max_cnt; i++)
      s_pos_arr[i] = bl2_long8_cvt(s_pos_arr[i]);
  else
    for (i=0; i<=max_cnt; i++)
      s_pos_arr[i] = (long)bl2_uint4_cvt(tmp_pos_arr[i]);
#else
  if (!mm64_flag) 
    for (i=0; i<=max_cnt; i++)
      s_pos_arr[i] = (long)tmp_pos_arr[i];
#endif

#ifdef DEBUG
  for (i=1; i<max_cnt-1; i++) {
    if (s_pos_arr[i+1]<s_pos_arr[i])
      fprintf(stderr," ** spos_error [%d]\t%ld\t%ld\n",
	      i,s_pos_arr[i],s_pos_arr[i]);
  }
#endif

  if (!mm64_flag) free(tmp_pos_arr);

  m_fd->max_cnt = max_cnt;
  m_fd->d_pos_arr = d_pos_arr;
  m_fd->s_pos_arr = s_pos_arr;
  m_fd->lpos = 0;

  /*   check_mmap(m_fd,-2); */

  return m_fd;
}  
 
char *mgets (char *s, int n, struct lmf_str *m_fd)
{
  char *cs, *mfp;

  mfp = m_fd->mmap_addr;
  cs = s;

  while (--n > 0 && (*mfp != (char)EOF))
    if ((*cs++ = *mfp++) == '\n') break;
  *cs = '\0';

  m_fd->mmap_addr = mfp;
  return (*mfp == (char)EOF && cs == s) ? NULL : s;
}

int
agetlibm(unsigned char *seq,
	 int maxs,
	 char *libstr,
	 int n_libstr,
	 fseek_t *libpos,
	 int *lcont,
	 struct lmf_str *m_fd,
	 long *l_off)
{
  register unsigned char *cp, *seqp;
  register int *ap;
  char *desc;
  int lpos;		/* entry number in library */
  long l;
  unsigned char *seqm, *seqm1;
  char *bp;
  static long seq_len;
  static unsigned char *cp_max;
#ifdef SUPERFAMNUM
  char *bp1, *bpa, *tp;
  int i;
#endif

  *l_off = 1;

  lpos = m_fd->lpos;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = m_fd->sascii;

  if (*lcont==0) {
    if (lpos >= m_fd->max_cnt) return (-1);
    seq_len = m_fd->d_pos_arr[lpos+1] - m_fd->s_pos_arr[lpos];
    if (seq_len < 0 || (seq_len > m_fd->max_len && seq_len > (m_fd->max_len*5)/4)) {
      fprintf(stderr," ** sequence over-run: %ld at %d\n",seq_len,lpos);
      return(-1);
    }
    *libpos = (fseek_t)lpos;

    desc = m_fd->mmap_base+m_fd->d_pos_arr[lpos]+1;
    strncpy(libstr,desc,n_libstr-1);
    libstr[n_libstr-1]='\0';
    if ((bp=strchr(libstr,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
    if (n_libstr > MAX_UID) {
      bp = libstr;
      while (*bp++) if ( *bp=='\001' || *bp=='\t') *bp=' ';
    }

    for (bp = desc; *bp && (*bp != '\n'); *bp++ )
      if (*bp == '@' && !strncmp(bp+1,"C:",2)) sscanf(bp+3,"%ld",l_off);

#ifdef SUPERFAMNUM
    sfnum[0]=nsfnum=0;
    strncpy(tline,desc,sizeof(tline));
    tline[MAXLINE-1]='\0';
    if ((bp=strchr(tline,'\n'))!=NULL) *bp='\0';
    if ((bp=strchr(tline,' ')) && (bp=strchr(bp+1,SFCHAR))) {
      if ((bpa = strchr(bp+1,'\001'))!=NULL) *bpa = '\0';
      if ((bp1=strchr(bp+1,SFCHAR))==NULL) {
	fprintf(stderr," second %c missing: %s\n",SFCHAR,tline);
      }
      else {
	*bp1 = '\0';
	i = 0;
	if ((tp = strtok(bp+1," \t"))!=NULL) {
	  sfnum[i++] = atoi(tp);
	  while ((tp = strtok((char *)NULL," \t")) != (char *)NULL) {
	    sfnum[i++] = atoi(tp);
	    if (i>=9) break;
	  }
	}
	sfnum[nsfnum=i]= 0;
	if (nsfnum>1) sf_sort(sfnum,nsfnum);
	else {
	  if (nsfnum<1) fprintf(stderr," found | but no sfnum: %s\n",libstr);
	}
      }
    }
    else {
      sfnum[0] = nsfnum = 0;
      }
#endif

    m_fd->mmap_addr = m_fd->mmap_base+m_fd->s_pos_arr[lpos];
    cp_max = (unsigned char *)(m_fd->mmap_addr+seq_len);
  }

  for (cp=(unsigned char *)m_fd->mmap_addr; seqp<seqm1; ) {
    if ((*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA) continue;
    --seqp;
    if (cp >= cp_max) break;
  }
  m_fd->mmap_addr = (char *)cp;

  if (seqp>=seqm1) (*lcont)++;
  else {
    *lcont=0;
    lpos++;
    m_fd->lpos = lpos;
  }
  *seqp = EOSEQ;
  /*   if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
aranlibm(char *str,
	 int cnt,
	 fseek_t libpos,
	 char *libstr,
	 struct lmf_str *m_fd)
{
  char *bp;
  int llen;
  int lpos;

  lpos = (int) libpos;

  llen = m_fd->s_pos_arr[lpos]-m_fd->d_pos_arr[lpos];
  if (llen >= cnt) llen = cnt-1;

  strncpy(str,m_fd->mmap_base+m_fd->d_pos_arr[lpos]+1,llen);
  str[llen]='\0';
  if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
  bp = str;
  while (*bp++) if ( *bp=='\001' || *bp=='\t') *bp=' ';
  m_fd->lpos = lpos;
}

/* there is no vgetlibm() because vgetlibm() and agetlibm() are
   identical - the difference in the two file formats relates to the
   location of the sequence, which is already available in spos_arr[].

   however vranlibm must accomodate both type 5 and 6 files;
   type 6 has extra stuff after the seq_id.
*/

void
vranlibm(char *str,
	 int cnt,
	 fseek_t libpos,
	 char *libstr,
	 struct lmf_str *m_fd)
{
  char *bp, *mp;
  int llen;
  int lpos;

  lpos = (int)libpos;

  llen = m_fd->s_pos_arr[lpos]-m_fd->d_pos_arr[lpos];

  mp = m_fd->mmap_base+m_fd->d_pos_arr[lpos];
  
  strncpy(str,mp+4,20);
  str[20]='\0';
  if ((bp=strchr(str,' '))!=NULL) *(bp+1) = '\0';
  else if ((bp=strchr(str,'\n'))!=NULL) *bp = ' ';
  bp = strchr(mp,'\n');

  llen -= (bp-mp)-5;
  if (llen >  cnt-strlen(str)) llen = cnt-strlen(str)-1;

  strncat(str,bp+1,llen);
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
  str[cnt-1]='\0';
  m_fd->lpos = lpos;
}

void
close_mmap(struct lmf_str *m_fd) {
  free(m_fd->s_pos_arr);
  free(m_fd->d_pos_arr);
  if (m_fd->mm_flg) {
    munmap(m_fd->mmap_base,m_fd->st_size);
    free(m_fd);
  }
  m_fd->mm_flg=0;
}  

#ifndef min
#define min(x,y) ((x) > (y) ? (y) : (x))
#endif

static int gcg_bton[4]={2,4,1,3};

int
gcg_getlibm(unsigned char *seq,
	    int maxs,
	    char *libstr,
	    int n_libstr,
	    fseek_t *libpos,
	    int *lcont,
	    struct lmf_str *m_fd,
	    long *l_off)
{
  char dummy[20];
  char gcg_date[6];
  char gcg_type[10];
  register unsigned char *cp, *seqp, stmp;
  register int *ap, lpos;
  unsigned char *seqm, *seqm1;
  long r_block, b_block, r_fact, r16_block;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = m_fd->sascii;
  lpos = m_fd->lpos; 

  if (*lcont==0) {
    if (lpos >= m_fd->max_cnt) return (-1);
    sscanf(m_fd->mmap_base+m_fd->d_pos_arr[lpos]+4,"%s %s %s %s %ld\n",
	   libstr,gcg_date,gcg_type,dummy,&(m_fd->gcg_len));

    m_fd->gcg_binary = (gcg_type[0]=='2');

    libstr[12]='\0';
    *libpos = lpos;
    m_fd->mmap_addr = m_fd->mmap_base+m_fd->s_pos_arr[lpos];
  }

  r_block = b_block = min((size_t)(seqm-seqp),m_fd->gcg_len);
  if (m_fd->gcg_binary) {
    r_block = (r_block+3)/4;
  }

  cp=(unsigned char *)m_fd->mmap_addr; 
  if (!m_fd->gcg_binary) {
    r_fact = 1;
    r16_block = r_block/16;
    while (r16_block-- > 0) {
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
      *seqp++ = ap[*cp++];
    }
    while (seqp<seq+r_block) *seqp++ = ap[*cp++];
  }
  else if (m_fd->gcg_binary) {
    r_fact = 4;
    r16_block = r_block/8;
    while(r16_block-- > 0) {
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
    }

    while (seqp < seq+4*r_block) {
      stmp = *cp++;
      *seqp++ = gcg_bton[(stmp>>6) &3];
      *seqp++ = gcg_bton[(stmp>>4) &3];
      *seqp++ = gcg_bton[(stmp>>2) &3];
      *seqp++ = gcg_bton[(stmp) &3];
    }
  }
  if (r_fact * r_block >= m_fd->gcg_len) {
    *lcont = 0;
    m_fd->lpos++;
  }
  else {
    if (m_fd->gcg_binary) b_block = 4*r_block;
    m_fd->gcg_len -= b_block;
    (*lcont)++;
  }

  seq[b_block] = EOSEQ;
  /*   if (b_block==0) return 1; else */
  return b_block;
}

void lget_ann_m(struct lmf_str *lm_fd, char *libstr, int n_libstr);

int
lgetlibm(unsigned char *seq,
	 int maxs,
	 char *libstr,
	 int n_libstr,
	 fseek_t *libpos,
	 int *lcont,
	 struct lmf_str *m_fd,
	 long *l_off)
{
  register unsigned char *cp, *seqp;
  register int *ap, lpos;
  unsigned char *seqm, *seqm1;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  lpos = m_fd->lpos;
  ap = m_fd->sascii;

  if (*lcont==0) {
    if (lpos >= m_fd->max_cnt) return (-1);

    if (n_libstr <= 21) {
      strncpy(libstr,m_fd->mmap_base+m_fd->d_pos_arr[lpos]+12,12);
      libstr[12]='\0';
    }
    else {
      lget_ann_m(m_fd,libstr,n_libstr);
    }
    *libpos = lpos;

    m_fd->mmap_addr = m_fd->mmap_base+m_fd->s_pos_arr[lpos];
    cp = (unsigned char *)m_fd->mmap_addr;
  }
  else cp = (unsigned char *)m_fd->mmap_addr;

  while (seqp<seqm1) {
    if (*cp=='/' && *(cp-1)=='\n') break;
    if ((*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA &&
	(*seqp++=ap[*cp++])<NA) continue;
    --seqp;
    if (*cp=='\n' && *(cp+1)==' ') cp += 11;
  }

  if (seqp>=seqm1) {
    (*lcont)++;
    m_fd->mmap_addr = (char *)cp;
  }
  else {
    *lcont=0;
    m_fd->lpos++;
  }

  *seqp = EOSEQ;
  return (int)(seqp-seq);
}

void
lget_ann_m(struct lmf_str *lm_fd, char *libstr, int n_libstr) {
  char *bp, *bp_gid, locus[120], desc[120], acc[120], ver[120];

  /* copy in locus from lm_fd->lline */
  strncpy(locus,&lm_fd->mmap_addr[12],sizeof(locus));
  if ((bp=strchr(locus,' '))!=NULL) *(bp+1) = '\0';

  /* get description */
  mgets(desc,sizeof(desc),lm_fd);
  while (desc[0]!='D' || desc[1]!='E' || strncmp(desc,"DEFINITION",10))
    mgets(desc,sizeof(desc),lm_fd);
  if ((bp = strchr(&desc[12],'\n'))!=NULL) *bp='\0';

  /* get accession */
  mgets(acc,sizeof(acc),lm_fd);
  while (acc[0]!='A' || acc[1]!='C' || strncmp(acc,"ACCESSION",9)) {
    mgets(acc,sizeof(acc),lm_fd);
    if (acc[0]=='O' && acc[1]=='R' && strncmp(acc,"ORIGIN",6)==0)
      break;
  }
  if ((bp = strchr(&acc[12],'\n'))!=NULL) *bp='\0';
  if ((bp = strchr(&acc[12],' '))!=NULL) *bp='\0';

  /* get version */
  mgets(ver,sizeof(ver),lm_fd);
  while (ver[0]!='V' || ver[1]!='E' || strncmp(ver,"VERSION",7)) {
    mgets(ver,sizeof(ver),lm_fd);
    if (ver[0]=='O' && ver[1]=='R' && strncmp(ver,"ORIGIN",6)==0)
      break;
  }
  if ((bp = strchr(&ver[12],'\n'))!=NULL) *bp='\0';

      /* extract gi:123456 from version line */
  bp_gid = strchr(&ver[12],':');
  if (bp_gid != NULL) {
    if ((bp=strchr(bp_gid+1,' '))!=NULL) *bp='\0';
    bp_gid++;
  }
  if ((bp = strchr(&ver[12],' '))!=NULL) *bp='\0';

      /* build up FASTA header line */
  if (bp_gid != NULL) {
    strncpy(libstr,"gi|",n_libstr-1);
    strncat(libstr,bp_gid,n_libstr-4);
    strncat(libstr,"|gb|",n_libstr-20);
  }
  else {libstr[0]='\0';}

  /* if we have a version number, use it, otherwise accession, 
	 otherwise locus/description */

  if (ver[0]=='V') {
    strncat(libstr,&ver[12],n_libstr-1-strlen(libstr));
    strncat(libstr,"|",n_libstr-1-strlen(libstr));
  }
  else if (acc[0]=='A') {
    strncat(libstr,&acc[12],n_libstr-1-strlen(libstr));
    strncat(libstr," ",n_libstr-1-strlen(libstr));
  }

  strncat(libstr,locus,n_libstr-1-strlen(libstr));
  strncat(libstr,&desc[11],n_libstr-1-strlen(libstr));
  libstr[n_libstr-1]='\0';
}

void
lranlibm(char *str,
	 int cnt,
	 fseek_t seek,
	 char *libstr,
	 struct lmf_str *m_fd)
{
  char *bp, *llp;
  char acc[MAXLINE], desc[MAXLINE];

  llp = m_fd->mmap_addr = m_fd->mmap_base + m_fd->d_pos_arr[seek];

  lget_ann_m(m_fd,str,cnt);

  str[cnt-1]='\0';

  m_fd->lpos = seek;
}

static int check_status=0;

void
check_mmap(struct lmf_str *m_fd,long ntt) {

  int i, seq_len, ok_stat;
  
  ok_stat = 1;
  if ( ++check_status > 5) return;

  fprintf(stderr," ** checking %s %ld**\n", m_fd->lb_name,ntt);
  for (i=0; i<m_fd->max_cnt; i++) {
    seq_len = m_fd->d_pos_arr[i+1] - m_fd->s_pos_arr[i];
    if (seq_len < 0 || (seq_len > m_fd->max_len && seq_len > (m_fd->max_len*5)/4)) {
      fprintf(stderr,"%d:\t%ld\t%ld\t%ld\n",
	      i,m_fd->d_pos_arr[i],m_fd->s_pos_arr[i],
	      m_fd->d_pos_arr[i+1]-m_fd->s_pos_arr[i]);
      ok_stat=0;
    }
  }
  if (ok_stat) {
    if (check_status) fprintf(stderr," ** check_mmap OK %s %ld**\n",
			      m_fd->lb_name,ntt);
  }
}

#ifdef DEBUG
/*  C H K 3  --  Compute a type-3 Kermit block check.  */
/*
 Calculate the 16-bit CRC of a null-terminated string using a byte-oriented
 tableless algorithm invented by Andy Lowry (Columbia University).  The
 magic number 010201 is derived from the CRC-CCITT polynomial x^16+x^12+x^5+1.
 Note - this function could be adapted for strings containing imbedded 0's
 by including a length argument.
*/
long
crck(s,n)
    char *s; int n;
{
    unsigned int c, q;
    long crc = 0;

    while (n-->0) {
	c = *s++;
	/* if (parity)*/
	c &= 0177;
	q = (crc ^ c) & 017;		/* Low-order nibble */
	crc = (crc >> 4) ^ (q * 010201);
	q = (crc ^ (c >> 4)) & 017;	/* High order nibble */
	crc = (crc >> 4) ^ (q * 010201);
    }
    return(crc);
}
#endif
