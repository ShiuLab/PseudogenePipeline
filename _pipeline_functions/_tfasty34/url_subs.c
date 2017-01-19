
/* copyright (c) 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: url_subs.c,v 1.9 2006/08/20 18:18:33 wrp Exp $ */

/* 30 Dec 2004 - modify REF_URL to accomodate current Entrez */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

#ifndef DEF_PROT_LIB
#define DEF_PROT_LIB "q"
#endif

#ifndef FASTA_HOST
#define FASTA_HOST "your.fasta.host.here/fasta/cgi"
#endif

void do_url1(FILE *fp, struct mngmsg m_msg, struct pstruct pst,
	     char *l_name, int n1, struct a_struct aln, long loffset)
{
  char my_l_name[200];
  char *db;
  char pgm[10], lib[MAX_FN];
  char *ref_url, *lbp=NULL;
  char *srch_url, *srch_url1;

  if (m_msg.ldnaseq==SEQT_DNA) db="nucleotide";
  else db="Protein";

  if (strncmp(m_msg.f_id0,"rss",3)==0) {
    strncpy(pgm,"fa",sizeof(pgm));
  }
  else if (strncmp(m_msg.f_id0,"rfx",3)==0) {
    strncpy(pgm,"fx",sizeof(pgm));
  }
  else { strncpy(pgm,m_msg.f_id0,sizeof(pgm)); }

  if (m_msg.lname[0]!='%') {
    strncpy(lib,m_msg.lname,sizeof(lib));
  }
  else {
    strncpy(lib,"%25",sizeof(lib));
    strncat(lib,&m_msg.lname[1],sizeof(lib));
  }
  lib[sizeof(lib)-1]='\0';

  strncpy(my_l_name,l_name,sizeof(my_l_name));
  my_l_name[sizeof(my_l_name)-1]='\0';

  if (pgm[0]=='t' || strcmp(pgm,"fx") || strcmp(pgm,"fy")==0 ) {
    if ((lbp=strchr(my_l_name,':'))!=NULL) *lbp='\0';
    lbp = &my_l_name[strlen(my_l_name)-2];
    if ( *lbp == '_' ) *lbp = '\0';
  }

  /* change the program name for fastx, tfastx, tfasta */
  /* fastx returns proteins */
  if (strcmp(pgm,"fx")==0 || strcmp(pgm,"fy")==0) strncpy(pgm,"fa",sizeof(pgm));
  else if (strcmp(pgm,"ff")==0) strncpy(pgm,"fa",sizeof(pgm));
  else if (pgm[0]=='t') {
    strncpy(pgm,"fx",sizeof(pgm));
    strncpy(lib,DEF_PROT_LIB,sizeof(lib));
  }

  fflush(fp);
  if ((ref_url = getenv("REF_URL"))==NULL)
    fprintf(fp,"<A HREF=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=%s&fcmd=Search&doptcmd1=DocSum&term=%s\">Entrez lookup</A>&nbsp;&nbsp;",
	    db,my_l_name);
  else
    fprintf(fp,ref_url,db,my_l_name);

  if ((srch_url = getenv("SRCH_URL"))==NULL)
    fprintf(fp,"<A HREF=\"http://%s/searchfa.cgi?query=%s&db=%s&lib=%s&pgm=%s&start=%ld&stop=%ld&n1=%d\">Re-search database</A>&nbsp;&nbsp;",
	    FASTA_HOST,my_l_name,db,lib,pgm,
	    loffset+aln.amin1+1,loffset+aln.amax1,n1);
  else 
    fprintf(fp,srch_url,my_l_name,db,lib,pgm,
	    loffset+aln.amin1+1,loffset+aln.amax1,n1);

  if ((srch_url1 = getenv("SRCH_URL1"))==NULL)
    fprintf(fp,"<A HREF=\"http://%s/searchxf.cgi?query=%s&db=%s&lib=%s&pgm=%s&start=%ld&stop=%ld&n1=%d\">General re-search</A>\n<p>\n",
	    FASTA_HOST,my_l_name,db,lib,pgm,
	    loffset+aln.amin1+1,loffset+aln.amax1,n1);
  else 
    fprintf(fp,srch_url1,my_l_name,db,lib,pgm,
	    loffset+aln.amin1+1,loffset+aln.amax1,n1);

  /* put back "_r"  */
  if (lbp!=NULL) *lbp = '_';

  /*
  if ((srch_url2 = getenv("SRCH_URL2"))==NULL)
    fprintf(fp,"<A HREF=\"http://fasta.bioch.virginia.edu/fasta/cgi/lalignx.cgi?seq1=\"%s\"&in_seq1=\"FASTA\"&seq2=\"%s\"&in_seq2=\"Accession\"&ssr2=%ld:%ld\">lalign</A>\n<p>\n",my_l_name,db,lib,pgm,loffset+aln.amin1+1,loffset+aln.amax1,n1);
  else 
    fprintf(fp,srch_url1,my_l_name,db,lib,pgm,
	    loffset+aln.amin1+1,loffset+aln.amax1,n1);
  */
  fflush(fp);

}
