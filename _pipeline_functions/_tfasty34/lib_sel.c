/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: lib_sel.c,v 1.16 2006/12/06 17:30:52 wrp Exp $ */

/*	modified Dec 13, 1989 requires different FASTLIBS */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "defs.h"
#include "structs.h"

#ifdef NCBIBL13
#define LASTLIB NCBIBL13+1
#else
#define LASTLIB 11
#endif


int getlnames(char *tname, struct mngmsg *m_msg);
void addfile(char *, char *, struct mngmsg *);
void libchoice(char *lname, int nl, struct mngmsg *m_msg);
void libselect(char *lname, struct mngmsg *m_msg);
void subs_env(char *dest, char *src, int dest_size);
char *ulindex(char *str, char *chr);

static char ldname[MAX_FN];
static char *libenv;

int
getlnames(char *iname, struct mngmsg *m_msg)	/* read in the library names */
{
  char *bp, tsave[MAX_STR], *tname;
  char lline[MAX_FN], *llp;
  FILE *tptr;

  /* expand environment variables */
  
  tname = tsave;
  subs_env(tname, iname, sizeof(tsave));

  if (*tname != '@') {addfile(tname,"\0",m_msg); return 1;}
  else tname++;

  /* remove ' ' before deftype if present */
  if ((bp=strchr(tname,' '))!=NULL) *bp='\0';

  if ((tptr=fopen(tname,"r"))==NULL) {
    fprintf(stderr," could not open file of names: %s\n",tname);
    return 0;
  }

  while (fgets(lline,sizeof(lline),tptr)!=NULL) {
    if (lline[0]==';') continue;
    if ((bp=strchr(lline,'\n'))!=NULL) *bp='\0';
    subs_env(tsave, lline, sizeof(tsave));
    if (tsave[0]=='<') {
      strncpy(ldname,&tsave[1],sizeof(ldname));
      ldname[sizeof(ldname)-1]='\0';
      libenv=ldname;
    }
    else addfile(tsave,libenv,m_msg);
  }
  fclose(tptr);
  return 1;
}

/* libchoice displays a list of potential library files
   in the new &lib& version, only traditional 1-letter files will be
   shown initially
*/

void
libchoice(char *lname, int nl, struct mngmsg *m_msg)
{
  FILE *fch;
  char line[MAX_STR], *bp;
  char *chstr[MAX_CH],*chfile[MAX_CH];
  char *chtmp, *charr;
  int i,j,k,chlen;

  charr = NULL;
  if (strlen(m_msg->flstr)> (size_t)0) {
    chlen = MAX_CH*MAX_FN;
    if ((chtmp=charr=calloc((size_t)chlen,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate choice file array\n");
      goto l1;
    }
    chlen--;
    if ((fch=fopen(m_msg->flstr,"r"))==NULL) {
      fprintf(stderr," cannot open choice file: %s\n",m_msg->flstr);
      goto l1;
    }
    fprintf(stderr,"\n Choose sequence library:\n\n");

    for (i=j=0; j<MAX_CH; i++) {
      if (fgets(line,sizeof(line),fch)==NULL) break;/* check for comment */
      if (line[0]==';') continue;
      if ((bp=strchr(line,'\n'))!=NULL) *bp='\0'; /* remove \n */
      if ((bp=strchr(line,'$'))==NULL) continue;  /* if no '$', continue */
      *bp++='\0';	      /* replace $ with \0, bp points to libtype */

      /* if libtypes don't match, continue */
      if ((*bp++ -'0')!=m_msg->ldnaseq) continue;

      /* if the library file name is too long, quit */
      if ((k=strlen(line))>chlen) break;

      /* save the library file name */
      strncpy(chstr[j]=chtmp,line,chlen);
      chtmp += k+1; chlen -= k+1;

      if ((k=strlen(bp))>chlen) break;
      strncpy(chfile[j]=chtmp,bp,chlen);
      chtmp += k+1; chlen -= k+1;
      fprintf(stderr,"    %c: %s\n",*chfile[j++],line);
    }
  l2:  fprintf(stderr,"\n Enter library filename (e.g. %s), letter (e.g. P)\n",
	       (m_msg->ldnaseq==0)? "prot.lib" : "dna.lib");
    fprintf(stderr," or a %% followed by a list of letters (e.g. %%PN): ");
    fflush(stderr);
    if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
    if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
    if (strlen(line)==0) goto l2;
    strncpy(lname,line,nl);
  }
  else {
  l1: fprintf(stderr," library file name: ");
    fflush(stderr);
    if (fgets(line,sizeof(line),stdin)==NULL) exit(0);
    if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';
    if (strlen(line)> (size_t)0) strncpy(lname,line,nl);
    else goto l1;
  }
  if (charr!=NULL) {
    fclose(fch);
    free(charr);
  }
}

/* libselect parses the choices in char *lname and builds the list
   of library files
*/
void
libselect(char *lname, struct mngmsg *m_msg)
{
  char line[MAX_FN*2], *bp, *bp1;
  char *llnames[MAX_LF]; /* pointers into new list of names */
  int new_abbr,ich, nch;	  /* use new multi-letter abbr */
  FILE *fch;

  new_abbr = 0;
  m_msg->nln = 0;
  if (strlen(lname) > (size_t)1 && *lname != '%' && *lname != '+') {
    getlnames(lname,m_msg); /* file name */ 
    return;
  }
  else {
    if (*m_msg->flstr=='\0') {
      fprintf(stderr," abbrv. list request but FASTLIBS undefined, cannot use %s\n",lname);
      exit(1);
    }

    if (strchr(lname,'+')) {
      /* indicates list of database abbrevs (not files) */
      new_abbr=1;
      nch = 0;
      bp = lname+1; if (*bp == '+') bp++;
      for (bp1=bp; bp!=NULL && bp1!=NULL; bp=bp1+1) {
	if ((bp1=strchr(bp,'+'))!=NULL) *bp1='\0';
	llnames[nch++] = bp;
      }
    }
    else if (*lname=='%') {     /* list of single letter abbreviations */
      lname++;	/* bump over '%' to get letters */
    }

  /* else just use a single character abbreviation */

  if (strlen(m_msg->flstr) > (size_t)0) {
    if ((fch=fopen(m_msg->flstr,"r"))==NULL) {
      fprintf(stderr," cannot open choice file: %s\n",m_msg->flstr);
      return;
    }
  }
  else {
    fprintf(stderr," FASTLIBS undefined\n");
    addfile(lname,"\0",m_msg);
    return;
  }

  /* read each line of FASTLIBS */
    while (fgets(line,sizeof(line),fch)!=NULL) { 
      if (line[0]==';') continue;	/* skip comments */
      if ((bp=strchr(line,'\n'))!=NULL) *bp='\0';	/* remove '\n' */
      if ((bp=strchr(line,'$'))==NULL) continue; /* no delim, continue */
      *bp++='\0';	/* point to library type */
      if ((*bp++ -'0')!=m_msg->ldnaseq) continue; /* doesn't match, continue */

      /* if !new_abbr, match on one letter with ulindex() */
      if (!new_abbr) {
	if (*bp=='+') continue; /* not a &lib& */
	else if (ulindex(lname,bp)!=NULL) { 
	  strncpy(m_msg->ltitle,line,MAX_FN);
	  getlnames(bp+1,m_msg);
	}
      }
      else {
	if (*bp!='+') continue;
	else {
	  bp++;
	  if ((bp1 = strchr(bp,'+'))!=NULL) {
	    *bp1='\0';
	    for (ich = 0; ich<nch; ich++) {
	      if (strcmp(llnames[ich],bp)==0) {
		strncpy(m_msg->ltitle,line,MAX_FN);
		getlnames(bp1+1,m_msg);
		break;
	      }
	    }
	    *bp1='+';
	  }
	  else fprintf(stderr,"%s missing final '+'\n",bp);
	}
      }
    }
    fclose(fch);
  }
}

void
addfile(char *fname, char *env, struct mngmsg *m_msg)
{
  char tname[MAX_STR], *bp, *bp1;
  char *lbptr;
  int len, lenv, l_size;

  /*  check for default directory for files  */
  if (env != NULL && *env != '\0') lenv = strlen(env)+1;
  else lenv = 0;

  len=strlen(fname)+1+lenv;

  if (lenv > 1 && *fname != '#') {	/* add default directory to file name */
    strncpy(tname,env,sizeof(tname)-1);
#ifdef UNIX
    strcat(tname,"/");
#endif
    }
  else tname[0]='\0';

  /* add fname to tname, allocate space, and move to space */
  strncat(tname,fname,sizeof(tname)-strlen(tname)-1);
  len=strlen(tname)+1;
  if ((lbptr=calloc(len,sizeof(char)))==NULL) {
    fprintf(stderr,"no more space for filenames: %s ignored\n",fname);
    return;
  }
  else {
    strncpy(lbptr,tname,len);
    lbptr[len-1]='\0';
  }

  if (m_msg->nln< MAX_LF) {
    m_msg->lbnames[m_msg->nln++]=lbptr;
  }
  else fprintf(stderr," no more file name slots: %s ignored\n",lbptr);
}

char *
ulindex(char *str, char *chr)
{
  char c;
 
  c = tolower((int)(*chr));

  while (*str != '\0' && tolower(*str) !=c ) str++;
  if (*str=='\0') return NULL;
  else return str;
}
