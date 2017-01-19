
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: mshowbest.c,v 1.44 2006/06/30 19:46:36 wrp Exp $ */

/*   29-Oct-2003 - changes so that bbp->cont < 0 => aa1 sequence is
     already in aa1, no re_openlib or re_getlib required
*/

/*   14-May-2003 Changes to use a more consistent coordinate numbering
     system for displays.  aln->d_start[01] is now consistently used
     to report the start of the alignment in all functions, and
     mshowbest.c has been modified to use d_start[01] instead of
     d_start[01]-1.  aln->min[01] now starts at 0 for all functions;
     instead of 1 for some functions (dropnfa.c, dropgsw.c, dropfs2.c
     earlier).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

#ifndef PCOMPLIB
#include "mm_file.h"
#include "mw.h"
#else
#include "p_mw.h"
#endif


#define MAX_BLINE 256

#ifndef PCOMPLIB
/* function calls necessary to re_getlib() the sequence and, do
   alignments, if necessary
*/

#define RANLIB (m_fptr->ranlib)

int
re_getlib(unsigned char *, int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

struct lmf_str *re_openlib(struct lmf_str *, int outtty);
#endif

extern void cal_coord(int n0, int n1, long sq0off, long loffset,
		      struct a_struct *aln);

void header_aux(FILE *);
void show_aux(FILE *, struct beststr *);
void w_abort (char *p, char *p1);

/* BBP_INFO get stuff directly from beststr or from beststr->desptr */
#ifdef PCOMPLIB
#define BBP_INFO(info) bbp->desptr->info
#else
#define BBP_INFO(info) bbp->info
#endif

extern double zs_to_bit(double, int, int);

/* showbest() shows a list of high scoring sequence descriptions, and
   their scores.  If -m 9, then an additional complete set of
   alignment information is provided.

   If PCOMPLIB or m_msg.quiet then the number of high scores to be
   shown is pre-determined by m_msg.mshow before showbest is called.

   The comp_lib.c version re_getlib()'s the sequence for its
   discription, and then does another alignment for -m 9 (Thus, it
   needs an f_str.  The PCOMPLIB version has everything available in
   beststr before showbest() is called.
*/

void showbest (FILE *fp, 
#ifndef PCOMPLIB
	       unsigned char **aa0, unsigned char *aa1, int maxn,
#endif
	       struct beststr **bptr,int nbest, int qlib, struct mngmsg *m_msg,
	       struct pstruct pst, struct db_str db,
	       char *gstring2
#ifndef PCOMPLIB
	       ,void **f_str
#endif
)
{
  int ntmp = 0;
  char bline[MAX_BLINE], fmt[40], pad[MAX_BLINE], rline[40];
  char l_name[128];
  int istart = 0, istop, ib;
  int nshow;
  int quiet;
  int r_margin;
  struct beststr *bbp;
  int n1tot;
  char *bp;
  char rel_label[12];
  char tmp_str[20], *seqc;
  int seqc_len;
  long loffset, l_off;
  int n0, n1;
  struct rstruct rst;
  int lc, maxc, nident, ngap;
  float percent, gpercent;
  struct a_struct *aln_p;
  int *tres;
  int gi_num;

#ifndef PCOMPLIB
  struct lmf_str *m_fptr;
#endif

  strncpy(rel_label,"\0",2);
#ifdef SHOWREL
  strncpy(rel_label," related",sizeof(rel_label));
#endif
#ifdef SHOWUN
  strncpy(rel_label," unrelated",sizeof(rel_label));
#endif
  rel_label[sizeof(rel_label)-1]='\0';

#ifdef PCOMPLIB
  quiet = 1;
#else
  quiet = m_msg->quiet;
#endif

  n0 = m_msg->n0;

  if (m_msg->aln.llen > MAX_BLINE) m_msg->aln.llen = MAX_BLINE;

  if (pst.zsflag < 0) r_margin = 10;
  else if (pst.zsflag>=0  && m_msg->srelv > 1 ) r_margin = 19;
  else r_margin = 10;

  if (m_msg->markx & MX_M9SUMM && m_msg->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
    r_margin += 15;
#else
    r_margin += 10;
#endif
  }

  if (m_msg->nframe < 0) 
#ifndef SUPERFAMNUM
    sprintf(fmt,"%%-%ds (%%4d)",m_msg->aln.llen-r_margin);
#else
    sprintf(fmt,"%%-%ds [%%4d](%%4d)",m_msg->aln.llen-(r_margin+4));
#endif
  else 
    sprintf(fmt,"%%-%ds (%%4d)",m_msg->aln.llen-(r_margin+4));

  memset(pad,' ',m_msg->aln.llen-(r_margin+6));
  pad[m_msg->aln.llen-(r_margin+12)]='\0';

  if (quiet != -1) {	/* quiet is set to -1 in comp_mlib.c to force
			   all significant hits to be shown */
    nshow = 20;
    if (m_msg->mshow == -1) nshow = nbest;		/* show all */
    /* show specified number */
    else if (m_msg->mshow_flg) {
      nshow = min (m_msg->mshow, nshow);
    }
  }
  else nshow = m_msg->nshow;

  if (quiet==0) istop = 20;
  else istop = nshow;

  if (quiet==0) {
    printf(" How many scores would you like to see? [%d] ",m_msg->nshow);
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    nshow = m_msg->nshow;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
    if (nshow<=0) nshow = min(20,nbest);
  }

  if ((bp = strchr (m_msg->qtitle, '\n')) != NULL) *bp = '\0';
/*   fprintf (fp, "%3d %s\n", qlib,m_msg->qtitle); */

  if (m_msg->markx & MX_HTML) fprintf(fp,"<p><tt><pre>\n");

  if (pst.zsflag >= 0) {
    if (bptr[0]->escore < m_msg->e_cut) {
      if (m_msg->z_bits==1) {/* show bit score */
	fprintf(fp,"The best%s scores are:%s%s bits E(%ld)",
	       rel_label,pad,m_msg->label,pst.zdb_size);
      }
      else {/* show z-score */
	fprintf(fp,"The best%s scores are:%s%s z-sc E(%ld)",
	      rel_label,pad,m_msg->label,pst.zdb_size);
      }
      header_aux(fp);
      if (m_msg->markx & MX_M9SUMM) {
	if (m_msg->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sim  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
	if (m_msg->markx & MX_HTML && m_msg->show_code !=1) { fprintf(fp,"<!-- ");}
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msg->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msg->f_id1);
#endif
	}
	if (m_msg->show_code == SHOW_CODE_ALIGN) {	fprintf(fp," aln_code"); }
	if (m_msg->markx & MX_HTML && m_msg->show_code!=1) { fprintf(fp," -->");}
      }
      fprintf(fp,"\n");
    }
    else {
      fprintf(fp,"!! No library sequences with E() < %.2g\n",m_msg->e_cut);
      m_msg->nshow = 0;
      if (m_msg->markx & MX_HTML) fprintf(fp,"<p></tt></pre>\n");
      return;
    }
  }
  else {
    fprintf(fp,"The best%s scores are:%s%s",rel_label,pad,m_msg->label);
    header_aux(fp);
    if (m_msg->markx & MX_M9SUMM) {
      if (m_msg->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	fprintf(fp," %%_id  %%_sm  alen");
#else
	fprintf(fp," %%_id  alen");
#endif
      }
      else {
#ifndef SHOWSIM
	fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msg->f_id1);
#else
	fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msg->f_id1);
#endif
      }
    }
    if (m_msg->show_code == SHOW_CODE_ALIGN) {	fprintf(fp," aln_code"); }
    fprintf(fp,"\n");
  }

  istart = 0;
l1:
  istop = min(nbest,nshow);
  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];
#ifdef SUPERFAMNUM
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msg->qsfnum_n,BBP_INFO(sfnum))) continue;
#ifdef SHOWUN
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msg->qsfnum,BBP_INFO(sfnum))) {
      istop = min(istop+1,nbest);
    /*
      fprintf(stderr,"skipping %d: %d==%d\n",ib,m_msg->qsfnum,BBP_INFO(sfnum));
      */
      continue;
    }
#endif
#ifdef SHOWREL
    if (BBP_INFO(nsfnum) == 0 || (BBP_INFO(nsfnum) > 0 && !sfn_cmp(m_msg->qsfnum,BBP_INFO(sfnum)))) {
      istop = min(istop+1,nbest);
    /*
      fprintf(stderr,"skipping %d: %d==%d\n",ib,m_msg->qsfnum,BBP_INFO(sfnum));
      */
      continue;
    }
#endif
#endif
    if (quiet==1 && pst.zsflag>=0) {
      if (bbp->escore > m_msg->e_cut) {
	nshow = ib;
	goto done;
      }
      else if (bbp->escore < m_msg->e_low) continue;
    }

#ifndef PCOMPLIB
    if ((m_fptr=re_openlib(bbp->m_file_p,!m_msg->quiet))==NULL) {
      fprintf(stderr,"*** cannot re-open %s\n",bbp->m_file_p->lb_name);
      exit(1);
    }
    RANLIB(bline,m_msg->aln.llen,bbp->lseek,bbp->libstr,m_fptr);
#else
  strncpy(bline,BBP_INFO(bline),m_msg->aln.llen-r_margin);
  bline[m_msg->aln.llen]='\0';
#endif

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  This requires that fast_pan use -m 6. */

  strncpy(l_name,bline,sizeof(l_name)); /* get rid of text after second "|" */
  l_name[sizeof(l_name)-1]='\0';
  if ((bp=strchr(l_name,' '))!=NULL) *bp=0;
  if ((bp=strchr(&l_name[3],'|'))!=NULL) *bp='\0';
  if (m_msg->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
  else if (m_msg->nframe > 0 && bbp->frame == 1)
    strncat(l_name,"_r",sizeof(l_name));
  if (bbp->cont-1 > 0) {
    sprintf(tmp_str,":%d",bbp->cont-1);
    strncat(l_name,tmp_str,sizeof(l_name)-strlen(l_name));
  }


#ifndef PCOMPLIB
  if (m_msg->stages>1 || m_msg->markx & MX_M9SUMM) {
    if (bbp->cont >= 0) {
      n1 = re_getlib(aa1,maxn,m_msg->maxt3,m_msg->loff,bbp->cont,m_msg->term_code,
		     &loffset,&l_off,bbp->m_file_p);
    }
    else { n1 = maxn;}
    if (! m_msg->markx & MX_M9SUMM) {
      do_opt (aa0[bbp->frame], m_msg->n0, aa1, n1, bbp->frame, &pst, f_str[bbp->frame], &rst);
      bbp->score[2]=rst.score[2];
    }
    else {
      bbp->sw_score = 
	do_walign(aa0[bbp->frame],m_msg->n0, aa1, n1, bbp->frame, 
		  &pst, f_str[bbp->frame], &bbp->a_res, &bbp->have_ares);

      
      /* save the alignment encoding for future use */
      if (bbp->have_ares && ((tres = calloc(bbp->a_res.nres+1,sizeof(int)))!=NULL)) {
	memcpy(tres,bbp->a_res.res,sizeof(int)*bbp->a_res.nres);
	bbp->a_res.res = tres;
      }

      aln_func_vals(bbp->frame, &m_msg->aln);

      maxc = bbp->a_res.nres + 4*m_msg->aln.llen+4;
      seqc = NULL;
      seqc_len = 0;
      if (m_msg->show_code == SHOW_CODE_ALIGN) {
	if ((seqc=(char *)calloc(maxc,sizeof(char)))!=NULL) {
	  lc=calc_code(aa0[bbp->frame],m_msg->n0,
		       aa1,n1, 
		       &m_msg->aln,bbp->a_res,
		       pst,seqc,maxc,f_str[bbp->frame]);
	  seqc_len = strlen(seqc);
	}
      }
      else {
	lc=calc_id(aa0[bbp->frame],m_msg->n0,aa1,n1,
		   &m_msg->aln, bbp->a_res,
		   pst,f_str[bbp->frame]);
      }
      m_msg->aln.a_len = lc;

      nident = m_msg->aln.nident;
      if (lc > 0) percent = (100.0*(float)nident)/(float)lc;
      else percent = -1.00;

      ngap = m_msg->aln.ngap_q + m_msg->aln.ngap_l;
#ifndef SHOWSIM
      if (lc-ngap > 0) gpercent = (100.0*(float)nident)/(float)(lc-ngap);
      else gpercent = -1.00;
#else
      if (lc-ngap > 0) gpercent = (100.0*(float)m_msg->aln.nsim)/(float)(lc);
      else gpercent = -1.00;
#endif

    }
  }
#endif

  n1tot = (BBP_INFO(n1tot_p)) ? *BBP_INFO(n1tot_p) : bbp->n1;

  bp = bline;
  if ((m_msg->markx & MX_HTML) && !strncmp(bline,"gi|",3)) {
    bp = strchr(bline+4,'|')+1;
    *(bp-1) = 0;
    gi_num = atoi(bline+3);
  }

#ifndef SUPERFAMNUM
  bp[m_msg->aln.llen-r_margin]='\0';
#else
  bp[m_msg->aln.llen-r_margin-5]='\0';
#endif

  if (m_msg->nframe == -1) bp[m_msg->aln.llen-r_margin]='\0';
  else bp[m_msg->aln.llen-(r_margin+4)]='\0';

#ifndef SUPERFAMNUM
  fprintf (fp, fmt,bp,n1tot);
#else
  if (m_msg->nframe == -1) {
    fprintf (fp, fmt,bp,BBP_INFO(sfnum[0]),n1tot);
  }
  else {fprintf (fp, fmt,bp,n1tot);}
#endif

  if (m_msg->nframe > 2) fprintf (fp, " [%d]", bbp->frame+1);
  else if (m_msg->nframe >= 0) fprintf(fp," [%c]",(bbp->frame > 0 ?'r':'f'));

  if (m_msg->srelv == 1) fprintf (fp, " %4d", bbp->score[pst.score_ix]);
  else {
    if (m_msg->srelv-1 > 0) fprintf (fp, " %4d", bbp->score[0]);
    if (m_msg->srelv-1 > 1 || m_msg->stages>1)
      fprintf (fp, " %4d", bbp->score[1]);
    fprintf (fp, " %4d", bbp->score[pst.score_ix]);
  }

  if (pst.zsflag>=0) { 
    if (m_msg->z_bits==1) {
      fprintf (fp, " %.1f %7.2g",zs_to_bit(bbp->zscore,m_msg->n0,bbp->n1),bbp->escore);
    }
    else fprintf (fp, " %.1f %7.2g",bbp->zscore,bbp->escore);
  }
  show_aux(fp,bbp);

#ifdef PCOMPLIB
  n1 = bbp->n1;
  percent = bbp->percent;
  gpercent = bbp->gpercent;
  aln_p = bbp->aln_d;
  seqc = bbp->aln_code;
  seqc_len = bbp->aln_code_n;
  loffset = bbp->desptr->loffset;
  l_off = 0;
#else
  aln_p = &(m_msg->aln);
#endif

  if (m_msg->markx & MX_M9SUMM) {
    if (m_msg->show_code != SHOW_CODE_ID) {
      if (m_msg->markx & MX_HTML) fprintf(fp,"<!-- ");
      cal_coord(m_msg->n0,bbp->n1,m_msg->sq0off,loffset+l_off-1,aln_p);

      /*            %_id  %_sim s-w alen an0  ax0  pn0  px0  an1  ax1  pn1  px1 gapq gapl fs  */
      /*                    alignment    min  max            min  max */
      /*                    sequence coordinate    min  max            min  max */
      fprintf(fp,"\t%5.3f %5.3f %4d %4d %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %3d %3d %3d",
	      percent/100.0,gpercent/100.0, bbp->sw_score,aln_p->a_len,
	      aln_p->d_start0,aln_p->d_stop0,
	      m_msg->sq0off, m_msg->sq0off+m_msg->n0-1,
	      aln_p->d_start1,aln_p->d_stop1,
	      loffset+l_off, loffset+l_off+bbp->n1-1,
	      aln_p->ngap_q,aln_p->ngap_l,aln_p->nfs);
      if (m_msg->show_code == SHOW_CODE_ALIGN
	  && seqc_len > 0 && seqc != NULL) {
	fprintf(fp,"\t%s",seqc);
      /*      fprintf(fp," [%2d:%d]",bbp->wrkr,bbp->seqnm); */
	free(seqc);
	seqc = NULL;
      }
      if (m_msg->markx & MX_HTML) fprintf(fp," -->");
    }
    else {
#ifdef SHOWSIM
      fprintf(fp," %5.3f %5.3f %4d", percent/100.0,(float)aln_p->nsim/(float)aln_p->a_len,aln_p->a_len);
#else
      fprintf(fp," %5.3f %4d", percent/100.0,aln_p->a_len);
#endif
    }
  }
  if (m_msg->markx & MX_HTML) fprintf(fp," <A HREF=\"#%s\">align</A>",l_name);
  fprintf (fp, "\n");
  fflush(fp);
  }

  if (quiet==0) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow += ntmp;
      goto l1;
    }
  }
  else if (quiet == 1)
    if (ib < nbest && (pst.zsflag>=0 && bbp->escore < m_msg->e_cut)) {
      if (m_msg->mshow_flg && istop >= m_msg->mshow) goto done;
      istart=istop;
      nshow += 10;
      goto l1;
    }

 done:
  m_msg->nshow = nshow;
  if (m_msg->markx & MX_HTML) fprintf(fp,"</pre></tt><p><hr><p>\n");
  if (fp!=stdout) fprintf(fp,"\n");
}

/*
  q[] has one set of sfnums, 0 terminated
  s[] has second
  return first match or 0
*/
