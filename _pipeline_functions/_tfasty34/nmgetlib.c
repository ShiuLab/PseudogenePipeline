
/* $Name: fa_34_26_5 $ - $Id: nmgetlib.c,v 1.35 2007/01/08 15:38:46 wrp Exp $ */

/*	May, June 1987	- modified for rapid read of database

	copyright (c) 1987,1988,1989,1992,1995,2000 William R. Pearson

	revised (split) version of nmgetaa.c -> renamed nmgetlib.c

	This version seeks to be a thread safe, no global, library
	reading program.  While adjusting the routines in this file
	should be relatively easy, ncbl2_mlib.c and mysql_lib.c may be
	more difficult.

	nmgetlib.c and mmgetaa.c are used together.  nmgetlib.c provides
	the same functions as nxgetaa.c if memory mapping is not used,
	mmgetaa.c provides the database reading functions if memory
	mapping is used. The decision to use memory mapping is made on
	a file-by-file basis.

	June 2, 1987 - added TFASTA
	March 30, 1988 - combined ffgetaa, fgetgb;
	April 8, 1988 - added PIRLIB format for unix
	Feb 4, 1989 - added universal subroutines for libraries
	December, 1995 - added range option file.name:1-1000
	September, 1999 - added option for mmap()ed files using ".xin" */


/*
	February 4, 1988 - this starts a major revision of the getaa
	routines.  The goal is to be able to seach the following format
	libraries:

	0 - normal FASTA format
	1 - full Genbank tape format
	2 - NBRF/PIR CODATA format
	3 - EMBL/Swiss-prot format
	4 - Intelligentics format
	5 - NBRF/PIR VMS format
	6 - GCG 2bit format

	11 - NCBI setdb/blastp (1.3.2) AA/NT
	12 - NCBI setdb/blastp (2.0) AA/NT
	16 - mySQL queries

	see file altlib.h to confirm numbers

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "structs.h"

#ifndef SFCHAR
#define SFCHAR ':'
#endif

#define EOSEQ 0

#include "uascii.h"
/* #include "upam.h" */

#define LFCHAR '\015'  /* for MWC 5.5 */

#include "altlib.h"


#include <fcntl.h>
#ifndef O_RAW
#ifdef O_BINARY
#define O_RAW O_BINARY
#else
#define O_RAW 0
#endif		/* O_BINARY */
#endif		/* O_RAW */

#ifdef WIN32
#define RBSTR "rb"	/* read file in binary mode */
#else
#define RBSTR "r"
#endif

#include "mm_file.h"
struct lmf_str *load_mmap(FILE *, char *, int, int, struct lmf_str *);
struct lmf_str *ncbl2_reopen(struct lmf_str *);
struct lmf_str *ncbl2_openlib(char *, int);

static struct lmf_str *last_m_fptr=NULL;

#ifdef MYSQL_DB
struct lmf_str *mysql_openlib(char *, int, int *);
struct lmf_str *mysql_reopen(struct lmf_str *);
#endif

#ifdef PGSQL_DB
struct lmf_str *pgsql_openlib(char *, int, int *);
struct lmf_str *pgsql_reopen(struct lmf_str *);
#endif

void closelib(struct lmf_str *m_fptr);
extern void newname(char *nname, char *oname, char *suff, int maxn);

/* a file name for openlib may include a library type suffix */

struct lmf_str *
openlib(char *lname, int ldnaseq, int *sascii,
	int outtty, struct lmf_str *om_fptr)
{
  char rline[10],sname[MAX_FN], iname[MAX_FN], *bp;
  char opt_text[MAX_FN];	/* save text after ':' */
  int wcnt, opnflg;
  int libtype;
  FILE *libi=NULL;
  FILE *libf;
  int use_stdin;
  struct lmf_str *m_fptr=NULL;

  /* this is currently unavailable - later it can return a value somewhere */
  /*
  if (lname[0]=='#') {return -9;}
  */

  if (om_fptr != NULL && om_fptr->mm_flg) {
    om_fptr->lpos = 0;
    return om_fptr;
  }

  wcnt = 0;	/* number of times to ask for file name */

  /* check to see if there is a file option ":1-100" */
#ifndef WIN32
  if ((bp=strchr(lname,':'))!=NULL && *(bp+1)!='\0') {
#else
  if ((bp=strchr(lname+3,':'))!=NULL && *(bp+1)!='\0') {
#endif
    strncpy(opt_text,bp+1,sizeof(opt_text));
    opt_text[sizeof(opt_text)-1]='\0';
    *bp = '\0';
  }
  else opt_text[0]='\0';

  if (lname[0] == '-' || lname[0] == '@') {
    use_stdin = 1;
  }
  else use_stdin=0;

  strncpy(sname,lname,sizeof(sname));
  sname[sizeof(sname)-1]='\0';
    /* check for library type */
  if ((bp=strchr(sname,' '))!=NULL) {
    *bp='\0';
    sscanf(bp+1,"%d",&libtype);
    if (libtype<0 || libtype >= LASTLIB) {
      fprintf(stderr," invalid library type: %d (>%d)- resetting\n%s\n",
	      libtype,LASTLIB,lname);
      libtype=0;
    }
  }
  else libtype=0;

  if (use_stdin && libtype !=0) {
    fprintf(stderr," @/- STDIN libraries must be in FASTA format\n");
    return NULL;
  }

  /* check to see if file can be open()ed? */

 l1:
  if (libtype<=LASTTXT) {
    if (!use_stdin) {
      opnflg=((libf=fopen(sname,RBSTR))!=NULL);
    }
    else {
      libf=stdin;
      strncpy(sname,"STDIN",sizeof(sname));
      sname[sizeof(sname)-1]='\0';
      opnflg=1;
    }
  } 
#ifdef NCBIBL13
  else if (libtype==NCBIBL13) opnflg=(ncbl_openlib(sname,ldnaseq)!= -1);
#endif
#ifdef NCBIBL20
  else if (libtype==NCBIBL20) {
    opnflg=((m_fptr=ncbl2_openlib(sname,ldnaseq))!=NULL);
  }
#endif

#ifdef MYSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (libtype==MYSQL_LIB) {
    opnflg=((m_fptr=mysql_openlib(sname,ldnaseq,sascii))!=NULL);
  }
#endif
#ifdef PGSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (libtype==PGSQL_LIB) {
    opnflg=((m_fptr=pgsql_openlib(sname,ldnaseq,sascii))!=NULL);
  }
#endif

  if (!opnflg) {	/* here if open failed */
    if (outtty) {
      fprintf(stderr," cannot open %s library\n",sname);
      fprintf(stderr," enter new file name or <RET> to quit ");
      fflush(stderr);
      if (fgets(sname,sizeof(sname),stdin)==NULL) return NULL;
      if ((bp=strchr(sname,'\n'))!=0) *bp='\0';
      if (strlen(sname)==0) return NULL;
      if (++wcnt > 10) return NULL;
      strncpy(lname,sname,sizeof(lname)-1);
      lname[sizeof(lname)-1]='\0';
      goto l1;
    }
    else return NULL;
  }	/* !openflg */

  if (libtype <= LASTTXT) {
    /* now allocate a buffer for the opened text file */
    if ((m_fptr = calloc(1,sizeof(struct lmf_str)))==NULL) {
      fprintf(stderr," cannot allocate lmf_str (%ld) for %s\n",
	      sizeof(struct lmf_str),sname);
      return NULL;
    }
    if ((m_fptr->lline = calloc(MAX_STR,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate lline (%d) for %s\n",
	      MAX_STR,sname);
      return NULL;
    }

    strncpy(m_fptr->lb_name,sname,MAX_FN);
    m_fptr->lb_name[MAX_FN-1]='\0';
    strncpy(m_fptr->opt_text,opt_text,MAX_FN);
    m_fptr->opt_text[MAX_FN-1]='\0';
    m_fptr->sascii = sascii;

    m_fptr->libf = libf;
    m_fptr->lb_type = libtype;
    m_fptr->getlib = getliba[libtype];
    m_fptr->ranlib = ranliba[libtype];
    m_fptr->mm_flg = 0;
    m_fptr->tot_len = 0;
    m_fptr->max_len = 0;
    m_fptr->lib_aa = (ldnaseq==0);
  }
  last_m_fptr = m_fptr;

#ifdef USE_MMAP
  /* check for possible mmap()ed files */
  if (!use_stdin && (libtype <= LASTTXT) && (getlibam[libtype]!=NULL)) {
    /* this is a file we can mmap() */
    /* look for .xin file */
    newname(iname,sname,"xin",sizeof(iname));
    if ((libi=fopen(iname,"r"))!=NULL) { /* have a *.xin file, use mmap */
      if (load_mmap(libi,sname,libtype,ldnaseq,m_fptr)!=NULL) {
	fclose(libi);	/* close index file */
	m_fptr->lb_type = libtype;
	m_fptr->getlib = getlibam[libtype];
	m_fptr->ranlib = ranlibam[libtype];
	m_fptr->mm_flg = 1;
	return m_fptr;
      }
    fclose(libi);	/* memory mapping failed, but still must close file */
    }
  }
#endif

  if (libtype <= LASTTXT) {
    m_fptr->lpos = 0;
    if (fgets(m_fptr->lline,MAX_STR,libf)==NULL) return NULL;
  }
  return m_fptr;
}

void
closelib(struct lmf_str *m_fptr) {


#ifdef MMAP
  if (m_fptr->mm_flag) {
/* don't close memory mapped files
    close_mmap(m_fptr);
*/
    return;
  }
#endif

  if (m_fptr->libf!=NULL && m_fptr->libf != stdin) {
    fclose(m_fptr->libf);
    m_fptr->libf = NULL;
  }

#ifdef NCBIBL13
  if (m_fptr->lb_type == NCBIBL13) ncbl_closelib(m_fptr);
#endif
#ifdef NCBIBL20
  if (m_fptr->lb_type == NCBIBL20) ncbl2_closelib(m_fptr);
#endif
#ifdef MYSQL_DB
  if (m_fptr->lb_type == MYSQL_LIB) mysql_closelib(m_fptr);
#endif
}

struct lmf_str *
re_openlib(struct lmf_str *om_fptr, int outtty)
{
  int opnflg;

  /* if the file mmap()ed and has been opened - use it and return */
  if (om_fptr->mm_flg) {
    return om_fptr;
  }
#ifdef MYSQL_DB
  /* if this is a mysql database - use it and return */
  else if (om_fptr->lb_type == MYSQL_LIB) {
    return om_fptr;
  }
#endif

  /* data is available, but file is closed or not memory mapped, open it */
  /* no longer check to memory map - because we could not do it before */

  opnflg = 1;
  if (om_fptr->lb_type<=LASTTXT && om_fptr->libf==NULL)
    opnflg=((om_fptr->libf=fopen(om_fptr->lb_name,RBSTR))!=NULL);
#ifdef NCBIBL13
  else if (om_fptr->lb_type==NCBIBL13)
    opnflg=(ncbl_openlib(om_fptr->lb_name,!om_fptr->lib_aa)!= -1);
#endif
#ifdef NCBIBL20
  else if (om_fptr->lb_type==NCBIBL20) {
    opnflg=((om_fptr=ncbl2_openlib(om_fptr->lb_name,!om_fptr->lib_aa))!=NULL);
  }
#endif
#ifdef MYSQL_DB
  /* a mySQL filename contains mySQL commands, not sequences */
  else if (om_fptr->lb_type==MYSQL_LIB) 
    opnflg=(mysql_reopen(om_fptr)!=NULL);
#endif

  if (!opnflg) {
    fprintf(stderr,"*** could not re_open %s\n",om_fptr->lb_name);
    return NULL;
  }

  /* use the old buffer for the opened text file */
  om_fptr->mm_flg = 0;
  last_m_fptr =  om_fptr;

  return om_fptr;
}

#ifdef SUPERFAMNUM
static char tline[512];
extern int nsfnum;	/* number of superfamily numbers */
extern int sfnum[10];	/* superfamily number from types 0 and 5 */
extern int nsfnum_n;
extern int sfnum_n[10];
#endif

void sf_sort(int *, int);

int
agetlib(unsigned char *seq, int maxs,
	char *libstr, int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int i;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  /* int ic, l_start, l_stop, l_limit, rn; */
  char *bp, *bp1, *bpa, *tp;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    *l_off = 1;
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';') {
      if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }
#ifdef SUPERFAMNUM
    strncpy(tline,lm_fd->lline+1,sizeof(tline));
    tline[sizeof(tline)-1]='\0';
    sfnum[0]=nsfnum=0;
    if ((bp=strchr(tline,' ')) && (bp=strchr(bp+1,SFCHAR))) {
      if ((bpa = strchr(bp+1,'\001'))!=NULL) *bpa = '\0';
      if ((bp1=strchr(bp+1,SFCHAR))==NULL) {
/*	fprintf(stderr," second %c missing: %s\n",SFCHAR,libstr); */
      }
      else {
	*bp1 = '\0';
	i = 0;
	if ((tp = strtok(bp+1," \t"))!=NULL) {
	  sfnum[i++] = atoi(tp);
	  while ((tp = strtok((char *)NULL," \t")) != (char *)NULL) {
	    if (isdigit(*tp)) sfnum[i++] = atoi(tp);
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

    if ((bp=strchr(lm_fd->lline,'@'))!=NULL && !strncmp(bp+1,"C:",2)) {
      sscanf(bp+3,"%ld",l_off);
    }

    strncpy(libstr,lm_fd->lline+1,n_libstr-1);
    libstr[n_libstr-1]='\0';
    if ((bp=strchr(libstr,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
    if (n_libstr > MAX_UID) {
      tp = libstr;
      while (*tp++) if (*tp == '\001' || *tp== '\t') *tp = ' ';
    }

    *libpos = lm_fd->lpos;

    /* make certain we have the end of the line */
    while (strchr((char *)lm_fd->lline,'\n')==NULL) {
      if (strlen(lm_fd->lline)<MAX_STR/2) 
	fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR/2,lm_fd->libf);
      else 
	fgets(&lm_fd->lline[MAX_STR/2],MAX_STR/2,lm_fd->libf);
    }
    lm_fd->lline[MAX_STR-1]='\0';
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),lm_fd->libf)!=NULL) {
    if (*seqp=='>') goto new;
    if (*seqp==';') {
      if (strchr((char *)seqp,'\n')==NULL) goto cont;
      continue;
    }

    /* removed - used for @P:1-n 
       if (l_limit) {
       for (cp=seqp; seqp<seqm1 && rn < l_stop && (ic=ap[*cp++])<EL; )
       if (ic < NA && ++rn > l_start) *seqp++ = (unsigned char)ic;
       if (rn > l_stop) goto finish;
       }
       else {
    */
    for (cp=seqp; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
    if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
  }
  goto done;
 new:
  strncpy(lm_fd->lline,(char *)seqp,MAX_STR);
  lm_fd->lline[MAX_STR-1]='\0';
  /* be certain to get complete line, if possible */
  if (strchr(lm_fd->lline,'\n')==NULL)
    fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
  lm_fd->lline[MAX_STR-1]='\0';
  if (strchr(lm_fd->lline,'\n')==NULL && strchr((char *)seqp,'\n')!=NULL)
    lm_fd->lline[strlen(lm_fd->lline)-1]='\n';
  goto done;

  /* removed - used for @P:1-n
finish: 
   while (lm_fd->lline[0]!='>' && 
	  fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
     if (lm_fd->libf != stdin) lm_fd->lpos = FTELL(lm_fd->libf);
   }
   goto done;
*/
 cont:
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  seqm1 = seqp;
 done:
  if (seqp>=seqm1) (*lcont)++;
  else {
    *lcont=0;
  }

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
aranlib(char *str, int cnt, fseek_t seek, char *libstr, struct lmf_str *lm_fd)
{
  char *bp;

  if (lm_fd->libf != stdin) {
    FSEEK(lm_fd->libf, seek, 0);
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

    if (lm_fd->lline[0]=='>' || lm_fd->lline[0]==';') {
      strncpy(str,lm_fd->lline+1,cnt);
      str[cnt-1]='\0';
      if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
      if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
      /*
	if ((bp = strchr(str,SFCHAR))!=NULL) *bp='\0';
	else if ((bp = strchr(str,'\001'))!=NULL) *bp='\0';
	else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
	else str[cnt-1]='\0';
      */
      bp = str;
      while (*bp++) if (*bp=='\001' || *bp=='\t') *bp=' ';
    }
    else {
      str[0]='\0';
    }
  }
  else str[0]='\0';
}

void lget_ann(struct lmf_str *, char *, int);

int
lgetlib(unsigned char *seq,
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *bp, *bp_gid;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='L' || lm_fd->lline[1]!='O' || 
	   strncmp(lm_fd->lline,"LOCUS",5)) { /* find LOCUS */
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    *libpos= lm_fd->lpos;

    if (n_libstr <= 21) {
      strncpy(libstr,&lm_fd->lline[12],12);
      libstr[12]='\0';
    }
    else {
      lget_ann(lm_fd,libstr,n_libstr);
      fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    }

    while (lm_fd->lline[0]!='O' || lm_fd->lline[1]!='R' ||
	   strncmp(lm_fd->lline,"ORIGIN",6)) { /* find ORIGIN */
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag) getc(lm_fd->libf);
    if (lm_fd->lline[0]=='/') goto new;
    for (cp= (unsigned char *)&lm_fd->lline[10]; seqp<seqm1; ) {
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
      if (*(--seqp)>NA) break;
    }
  }
  goto done;
new:
  lm_fd->lpos = FTELL(lm_fd->libf);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

done:
  if (seqp>=seqm1) {
    lm_fd->cpsave = cp;
    (*lcont)++;
  }
  else *lcont=0;

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
lget_ann(struct lmf_str *lm_fd, char *libstr, int n_libstr) {
  char *bp, *bp_gid, locus[120], desc[120], acc[120], ver[120];

  /* copy in locus from lm_fd->lline */
  strncpy(locus,&lm_fd->lline[12],sizeof(locus));
  if ((bp=strchr(locus,' '))!=NULL) *(bp+1) = '\0';

  /* get description */
  fgets(desc,sizeof(desc),lm_fd->libf);
  while (desc[0]!='D' || desc[1]!='E' || strncmp(desc,"DEFINITION",10))
    fgets(desc,sizeof(desc),lm_fd->libf);
  if ((bp = strchr(&desc[12],'\n'))!=NULL) *bp='\0';

  /* get accession */
  fgets(acc,sizeof(acc),lm_fd->libf);
  while (acc[0]!='A' || acc[1]!='C' || strncmp(acc,"ACCESSION",9)) {
    fgets(acc,sizeof(acc),lm_fd->libf);
    if (acc[0]=='O' && acc[1]=='R' && strncmp(acc,"ORIGIN",6)==0)
      break;
  }
  if ((bp = strchr(&acc[12],'\n'))!=NULL) *bp='\0';
  if ((bp = strchr(&acc[12],' '))!=NULL) *bp='\0';

  /* get version */
  fgets(ver,sizeof(ver),lm_fd->libf);
  while (ver[0]!='V' || ver[1]!='E' || strncmp(ver,"VERSION",7)) {
    fgets(ver,sizeof(ver),lm_fd->libf);
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


/* this code seeks to provide both the various accession numbers
   necessary to identify the sequence, and also some description.

   Unfortunately, the various contributors to Genbank use three
   slightly different formats for including the accession number.

(1)LOCUS       HSJ214M20  107422 bp    DNA             HTG       16-JUN-2000
   DEFINITION  Homo sapiens chromosome 6 clone RP1-214M20 map p12.1-12.3, ***
               SEQUENCING IN PROGRESS ***, in unordered pieces.
   ACCESSION   AL121969

(2)LOCUS       AL359201   117444 bp    DNA             HTG       15-JUN-2000
   DEFINITION  Homo sapiens chromosome 1 clone RP4-671C13 map p13.2-21.1, ***
               SEQUENCING IN PROGRESS ***, in unordered pieces.
   ACCESSION   AL359201

(3)LOCUS       BB067000      280 bp    mRNA            EST       19-JUN-2000
   DEFINITION  BB067000 RIKEN full-length enriched, 15 days embryo male testis Mus
               musculus cDNA clone 8030456L01 3', mRNA sequence.
   ACCESSION   BB067000

This makes it more difficult to both provide the accession number in a
standard location and to conserve definition space
*/

void
lranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp, acc[MAX_STR], desc[MAX_STR];

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  lget_ann(lm_fd, str, cnt);
  str[cnt-1]='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

int
pgetlib(unsigned char *seq,
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int ic;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='E' || lm_fd->lline[1]!='N' || strncmp(lm_fd->lline,"ENTRY",5))
      { /* find ENTRY */
	lm_fd->lpos = FTELL(lm_fd->libf);
	if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      }
    strncpy(libstr,&lm_fd->lline[16],8);
    libstr[8]='\0';
    *libpos = lm_fd->lpos;
    while (lm_fd->lline[2]!='Q' || lm_fd->lline[0]!='S' || strncmp(lm_fd->lline,"SEQUENCE",8))
      { /* find SEQUENCE */
	if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      }
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf); /* get the extra line */
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lline[0]=='/') goto new;
    for (cp= (unsigned char *)&lm_fd->lline[8]; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    };
    if (*seqp==ES) goto done;
  }
  goto done;
new:
  lm_fd->lpos = FTELL(lm_fd->libf);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

done:
  if (seqp>=seqm1) {
    lm_fd->cpsave = cp;
    (*lcont)++;
  }
  else *lcont=0;

  *seqp = EOSEQ;
  /*  if ((int)(seqp-seq)==0) return 1; */
  return (int)(seqp-seq);
}

void
pranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  strncpy(str,&lm_fd->lline[16],8);
  str[8]='\0';
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  while (lm_fd->lline[0]!='T' || lm_fd->lline[1]!='I' || strncmp(lm_fd->lline,"TITLE",5))
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  strncpy(&str[8],&lm_fd->lline[16],cnt-9);
  str[cnt-9]='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
}

int
egetlib(unsigned char *seq,
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int ll;
  int ic;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char id[11];  /* Holds Identifier */

  *l_off=1;

  seqp = seq;
  seqm = &seq[maxs-11];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='I' || lm_fd->lline[1]!='D') { /* find ID */
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    sscanf(&lm_fd->lline[5],"%s",id);
    sprintf(libstr,"%-12.12s",id);
    libstr[12]='\0';
    *libpos = lm_fd->lpos;
    while (lm_fd->lline[0]!='S' || lm_fd->lline[1]!='Q') { /* find ORIGIN */
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }
    sscanf(&lm_fd->lline[14],"%ld",&lm_fd->gcg_len);
  }
  else {
    for (cp= lm_fd->cpsave; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets(lm_fd->lline,MAX_STR,lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag) getc(lm_fd->libf);
    if (lm_fd->lline[0]=='/') goto new;
    lm_fd->lline[70]='\0';
    for (cp= (unsigned char *)&lm_fd->lline[5]; seqp<seqm1; ) {
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
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
  }
  goto done;
new:	lm_fd->lpos = FTELL(lm_fd->libf);
fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
if (lm_fd->lfflag) getc(lm_fd->libf);
goto done;

done:	if (seqp>=seqm1) {
  lm_fd->cpsave = cp;
  (*lcont)++;
  lm_fd->gcg_len -= (long)(seqp-seq);
}
else *lcont=0;

*seqp = EOSEQ;
/* if ((int)(seqp-seq)==0) return 1; */
/*	if (*lcont==0 && (long)(seqp-seq)!=lm_fd->gcg_len)
	printf("%s read %d of %d\n",libstr,(int)(seqp-seq),lm_fd->gcg_len);
	*/
return (int)(seqp-seq);
}

void
eranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp;
  char id[11];  /* Holds Identifier */

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  sscanf(&lm_fd->lline[5],"%s",id);
  sprintf(str,"%-10.10s ",id);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
  while (lm_fd->lline[0]!='D' || lm_fd->lline[1]!='E') fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  strncpy(&str[11],&lm_fd->lline[5],cnt-11);
  str[cnt-11]='\0';
  if ((bp = strchr(str,'\r'))!=NULL) *bp='\0';
  if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

int
igetlib(unsigned char *seq,
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
	register unsigned char *cp, *seqp;
	register int *ap;
	unsigned char *seqm, *seqm1;
	char *bp;

	*l_off = 1;

	seqp = seq;
	seqm = &seq[maxs-9];
	seqm1 = seqm-1;

	ap = lm_fd->sascii;

	if (*lcont==0) {
		while (lm_fd->lline[0]!=';') {
			lm_fd->lpos = FTELL(lm_fd->libf);
			if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
			}
		*libpos = lm_fd->lpos;
		while (lm_fd->lline[0]==';') fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
		strncpy(libstr,lm_fd->lline+1,12);
		libstr[12]='\0';
		if((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
		}

	lm_fd->lline[0]='\0';
	while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),lm_fd->libf)!=NULL) {
		if (*seqp=='>') goto new;
		if (*seqp==';') {
			if (strchr((char *)seqp,'\n')==NULL) goto cont;
			continue;
			}
		for (cp=seqp; seqp<seqm1; ) {
			if ((*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
		            (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA &&
			    (*seqp++=ap[*cp++])<NA) continue;
			    if (*(--seqp)>NA) break;
			    }
		if (*seqp==ES) goto done;
		lm_fd->lpos = FTELL(lm_fd->libf);
		}
	goto done;
new:	strncpy(lm_fd->lline,(char *)seqp,MAX_STR);
	lm_fd->lline[MAX_STR-1]='\0';
	if (strchr((char *)seqp,'\n')==NULL)
	    fgets(lm_fd->lline,MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
	goto done;

cont:
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	seqm1 = seqp;

done:	if (seqp>=seqm1) {
		(*lcont)++;
		}
	else {
	*lcont=0;
		}


	*seqp = EOSEQ;
	/*	if ((int)(seqp-seq)==0) return 1; */
	return (int)(seqp-seq);
	}

void
iranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
	char *bp;
	char tline[MAX_FN];

	FSEEK(lm_fd->libf, seek, 0);
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

	if (lm_fd->lline[0]=='>' || lm_fd->lline[0]==';') {
		strncpy(tline,lm_fd->lline+1,sizeof(tline));
		tline[sizeof(tline)-1]='\0';
		if ((bp = strchr(tline,'\n'))!=NULL) *bp='\0';
		}
	else {
		tline[0]='\0';
		}

	while (lm_fd->lline[0]==';') fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	if ((bp=strchr(lm_fd->lline,'\n'))!=NULL) *bp=0;
	if ((bp=strchr(lm_fd->lline,' '))!=NULL) *bp=0;
	strncpy(str,lm_fd->lline,cnt);
	str[cnt-1]='\0';
	strncat(str,"  ",cnt-strlen(str)-1);
	strncat(str,tline,cnt-strlen(str)-1);
	
	FSEEK(lm_fd->libf,seek,0);
	fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
	}

int
vgetlib(unsigned char *seq,
	int maxs,
	char *libstr,
	int n_libstr,
	fseek_t *libpos,
	int *lcont,
	struct lmf_str *lm_fd,
	long *l_off)
{
  int i, ich;
  register unsigned char *cp, *seqp;
  register int *ap;
  unsigned char *seqm, *seqm1;
  char *bp, *tp;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';') {
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
      if (lm_fd->lfflag) getc(lm_fd->libf);
    }

#ifdef SUPERFAMNUM
    if ((bp=strchr(&lm_fd->lline[1],' ')) &&
	 (bp=strchr(bp+1,SFCHAR))) {
      i=0;
      if ((tp = strtok(bp+1," \t\n"))!=NULL) sfnum[i++] = atoi(tp);
      while ((tp = strtok(NULL," \t")) != NULL) {
	sfnum[i++] = atoi(tp);
	if (i>=10) break;
      }
      sfnum[nsfnum=i]= 0;
      if (nsfnum>1) sf_sort(sfnum,nsfnum);
      else {
	if (nsfnum < 1) fprintf(stderr," found | but no sfnum: %s\n",libstr);
      }
    }
    else sfnum[0]=nsfnum=0;
#endif

    if ((bp=strchr(lm_fd->lline,'\n'))!=NULL) *bp='\0';
    strncpy(libstr,&lm_fd->lline[4],12);
    libstr[12]='\0';
    if ((bp=strchr(libstr,' '))!=NULL) *bp='\0';
    if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
    
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    if (lm_fd->lfflag) getc(lm_fd->libf);

    if (n_libstr > 21) {
      strcat(libstr," ");
      strncat(libstr,lm_fd->lline,n_libstr-1-strlen(libstr));
      if ((bp=strchr(libstr,'\n'))!=NULL) *bp='\0';
      libstr[n_libstr-1]='\0';
    }
    *libpos = lm_fd->lpos;
  }

  lm_fd->lline[0]='\0';
  while (seqp<seqm1 && fgets((char *)seqp,(size_t)(seqm-seqp),lm_fd->libf)!=NULL) {
    if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
    if (*seqp=='>') goto new;
    if (*seqp==';') {
      if (strchr((char *)seqp,'\n')==NULL) goto cont;
      continue;
    }
    for (cp=seqp; seqp<seqm1; ) {
      if ((*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA &&
	  (*seqp++=ap[*cp++])<NA) continue;
      if (*(--seqp)>NA) break;
    }
    if (*seqp==ES) goto done;
    lm_fd->lpos = FTELL(lm_fd->libf);
  }
  goto done;
new:
  strncpy(lm_fd->lline,(char *)seqp,MAX_STR);
  lm_fd->lline[MAX_STR-1]='\0';
  if (strchr((char *)seqp,'\n')==NULL) {
    fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR-strlen(lm_fd->lline),lm_fd->libf);
    if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
  }
  goto done;

cont:
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag && (ich=getc(lm_fd->libf))!=LFCHAR) ungetc(ich,lm_fd->libf);
  seqm1 = seqp;

done:
  if (seqp>=seqm1) {
    (*lcont)++;
  }
  else {
    *lcont=0;
  }

  *seqp = EOSEQ;
  /*   if ((int)(seqp-seq)==0) return 1;*/
  return (int)(seqp-seq);
}

void
vranlib(char *str,
	int cnt,
	fseek_t seek,
	char *libstr,
	struct lmf_str *lm_fd)
{
  char *bp, *llp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);

  if (lm_fd->lline[0]=='>'&&(lm_fd->lline[3]==';'||lm_fd->lline[3]=='>')) {
    strncpy(str,&lm_fd->lline[4],cnt-1);
    str[cnt-1]='\0';

    if ((bp = strchr(str,':'))!=NULL) *bp='\0';
    if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    if (lm_fd->lfflag) getc(lm_fd->libf);

    /* skip over redundant stuff */
    for (llp=lm_fd->lline,bp=str; *llp==*bp; llp++,bp++);
    if ((int)(llp-lm_fd->lline)<5) llp = lm_fd->lline;

    if ((bp=strchr(llp,'\r'))!=NULL) *bp=' ';
    if ((bp=strchr(llp,'\n'))!=NULL) *bp='\0';
    strncat(str," ",(size_t)1);
    strncat(str,llp,(size_t)cnt-strlen(str)-1);
  }
  else {
    str[0]='\0';
  }

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
  if (lm_fd->lfflag) getc(lm_fd->libf);
}

static int gcg_bton[4]={2,4,1,3};

int
gcg_getlib(unsigned char *seq,
	   int maxs,
	   char *libstr,
	   int n_libstr,
	   fseek_t *libpos,
	   int *lcont,
	   struct lmf_str *lm_fd,
	   long *l_off)
{
  char dummy[20];
  char gcg_date[10];
  register unsigned char *cp, *seqp, stmp;
  register int *ap;
  char gcg_type[10];
  unsigned char *seqm, *seqm1;
  long r_block, b_block;
  char *bp;

  *l_off = 1;

  seqp = seq;
  seqm = &seq[maxs-9];
  seqm1 = seqm-1;

  ap = lm_fd->sascii;

  if (*lcont==0) {
    while (lm_fd->lline[0]!='>' && lm_fd->lline[0]!=';') {
      lm_fd->lpos = FTELL(lm_fd->libf);
      if (fgets(lm_fd->lline,MAX_STR,lm_fd->libf)==NULL) return (-1);
    }
    sscanf(&lm_fd->lline[4],"%s %s %s %s %ld",
	   libstr,gcg_date,gcg_type,dummy,&(lm_fd->gcg_len));

    lm_fd->gcg_binary = (gcg_type[0]=='2');

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    while (strchr((char *)lm_fd->lline,'\n')==NULL) {
      if (strlen(lm_fd->lline)<MAX_STR/2) 
	fgets(&lm_fd->lline[strlen(lm_fd->lline)],MAX_STR/2,lm_fd->libf);
      else 
      fgets(&lm_fd->lline[strlen(lm_fd->lline)-MAX_STR/2],MAX_STR/2,lm_fd->libf);
    }
    lm_fd->lline[MAX_STR-1]='\0';
    if (n_libstr <= 21) {
      libstr[12]='\0';
    }
    else {
      strncat(libstr," ",1);
      strncat(libstr,lm_fd->lline,n_libstr-1-strlen(libstr));
      if ((bp = strchr(libstr,'\n'))!=NULL) *bp='\0';
      libstr[n_libstr-1]='\0';
    }
    *libpos = lm_fd->lpos;
  }

  lm_fd->lline[0]='\0';

  r_block = b_block = min((size_t)(seqm-seqp),lm_fd->gcg_len);
  if (lm_fd->gcg_binary) { r_block = (r_block+3)/4; }

  fread((char *)seqp,(size_t)r_block,(size_t)1,lm_fd->libf);
  if (!lm_fd->gcg_binary) 
    for (cp=seqp; seqp<seq+r_block; ) *seqp++ = ap[*cp++];
  else if (lm_fd->gcg_binary) {
    seqp = seq + r_block;
    cp = seq + 4*r_block;
    while (seqp > seq) {
      stmp = *--seqp;
      *--cp = gcg_bton[stmp&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
      *--cp = gcg_bton[(stmp >>= 2)&3];
    }
  }
  if (4 * r_block >= lm_fd->gcg_len) {
    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
    *lcont = 0;
  }
  else {
    if (lm_fd->gcg_binary) b_block = 4*r_block;
    lm_fd->gcg_len -= b_block;
    (*lcont)++;
  }

  seq[b_block] = EOSEQ;
  /*   if (b_block==0) return 1; else */
  return b_block;
}

void
gcg_ranlib(char *str,
	   int cnt,
	   fseek_t seek,
	   char *libstr,
	   struct lmf_str *lm_fd)
{
  char *bp, *bp1, *llp;

  FSEEK(lm_fd->libf, seek, 0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

  if (lm_fd->lline[0]=='>'&&(lm_fd->lline[3]==';'||lm_fd->lline[3]=='>')) {
    strncpy(str,&lm_fd->lline[4],cnt-1);
    str[cnt-1]='\0';
    if ((bp = strchr(str,' '))!=NULL) *bp='\0';
    else if ((bp=strchr(str,'\r'))!=NULL) *bp='\0';
    else if ((bp = strchr(str,'\n'))!=NULL) *bp='\0';
    else str[cnt-1]='\0';

    fgets(lm_fd->lline,MAX_STR,lm_fd->libf);

    /* check beginning of line it is a duplicate */
    for (llp=lm_fd->lline,bp=str; *llp == *bp; llp++,bp++);
    if ((int)(llp-lm_fd->lline)<5) llp = lm_fd->lline;

    /* here we would like to skip over some species stuff */
	/*
    if ((bp1 = strchr(llp,';'))!=NULL && (int)(bp1-llp)<50) {
      if ((bp2 = strchr(bp1+1,';'))!=NULL && (int)(bp2-bp1)<50) {
	*(bp2+1)='\0'; bp1 = bp2+2;
      }
      else {bp1=llp;}
    }
    else if ((bp1=strchr(llp,'.'))!=NULL && *(bp1+1)==' ') {
      *(bp1+1) = '\0'; bp1 += 2;}
    else bp1 = llp;
    */
    
    bp1 = llp;
    if ((bp=strchr(bp1,'\r'))!=NULL) *bp='\0';
    if ((bp=strchr(bp1,'\n'))!=NULL) *bp='\0';
    strncat(str," ",(size_t)1);
    strncat(str,bp1,(size_t)cnt-strlen(str));
    if (bp1!=llp) strncat(str,llp,(size_t)cnt-strlen(str));
  }
  else {
    str[0]='\0';
  }

  FSEEK(lm_fd->libf,seek,0);
  fgets(lm_fd->lline,MAX_STR,lm_fd->libf);
}

void
sf_sort(s,n)
     int *s, n;
{
  int gap, i, j;
  int itmp;

  if (n == 1) return;
	
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
