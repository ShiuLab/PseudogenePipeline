
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: mshowalign.c,v 1.43 2007/01/08 15:38:46 wrp Exp $ */

/* mshowalign.c - show sequence alignments in pvcomplib */

/* 
   this is a merged version of showalign.c that works properly with
   both the comp_lib (serial, threaded) and PCOMPLIB parallel versions
   of the programs.

   In the serial and current threaded versions of the programs,
   showalign gets a list of high scoring sequences and must
   re_getlib() the sequence, do_walign(), and then calculate the
   alignment.

   In the PCOMPLIB parallel versions, the worker programs do the
   aligning, so showalign() must send them the appropriate messages to
   have the alignment done, and then collect the alignment results

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "msg.h"
#include "structs.h"
#include "param.h"

#ifdef PCOMPLIB
#ifdef PVM_SRC
#include "pvm3.h"
extern int pinums[];
#endif
#ifdef MPI_SRC
#include "mpi.h"
#endif
#include "p_mw.h"
#else
#include "mm_file.h"
#include "mw.h"
#endif

#ifndef PCOMPLIB

/* used to position the library sequence for re_getlib - also gets
   description */
#define RANLIB (m_fptr->ranlib)

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

int
re_getlib(unsigned char *aa1, int maxn, int maxt,
	  int loff, int cont, int term_code,
	  long *loffset, long *l_off, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

#endif


extern void cal_coord(int n0, int n1, long sq0off, long loffset,
		      struct a_struct *aln);

void initseq(char **, char **, char **, char **, int);
void freeseq(char **, char **, char **, char **);

void do_show(FILE *fp, int n0, int n1, int score,
	     char *name0, char *name1, int nml,
	     struct mngmsg m_msg, struct pstruct pst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqca, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln, long loffset);

extern void discons(FILE *fd, struct mngmsg m_msg, struct pstruct pst,
		    char *seqc0, char *seqc0a, char *seqc1, char *seqca,
		    int nc, 
		    int n0, int n1, char *name0, char *name1, int nml,
		    struct a_struct *aln,
		    long loffset);

extern void disgraph(FILE *fd, int n0, int n1,
		     float percent, int score,
		     int min0, int min1, int max0, int max1, long sq0off,
		     char *name0, char *name1, int nml, int llen, int markx);

extern double zs_to_bit(double, int, int);

extern void
do_url1(FILE *, struct mngmsg, struct pstruct,	char *, int,
	struct a_struct , long);

#ifndef A_MARK
#define A_MARK ">>"
#endif

static char l_name[200];	/* link name */

#ifdef PCOMPLIB
#define BBP_INFO(info) bbp->desptr->info
#else
#define BBP_INFO(info) bbp->info
#endif

/* this version does not check for m_msg.e_cut because nshow/nbest has
   already been set to limit on e_cut */

void showalign (FILE *fp,
#ifndef PCOMPLIB
		unsigned char **aa0, unsigned char *aa1, int maxn,
#endif
		struct beststr **bptr, int nbest, int qlib, 
		struct mngmsg m_msg, struct pstruct pst, char *gstring2
#ifndef PCOMPLIB
		, void **f_str
#endif
)
{
  char tmp_str[20];
  char info_str[200];
  char bline[2048], *bl_ptr, *bp, fmt[40];
  int tmp_len, l_llen;
  int t_have_ares;
  char name0[80], name0s[80], name1[200];
  int istart = 0, istop, i = 0, ib, nml;
  int  n1tot;
  struct beststr *bbp;
  int nc, lc, maxc;
  float percent, gpercent;
  char *seqc0, *seqc0a, *seqc1, *seqca;
  long loffset, l_off;
#ifdef PCOMPLIB
  struct stage2_str liblist;
  struct qmng_str qm_msg;
#ifdef MPI_SRC
  int int_msg_b[10];
  MPI_Status mpi_status;
#endif
#else
  int n1;
  struct lmf_str *m_fptr;
  int ngap;
#endif

#ifdef PCOMPLIB
  /* this function has its own copy of qm_msg, so we must fill it
     appropriately */
  qm_msg.n0 = m_msg.n0;
  strncpy(qm_msg.libstr,m_msg.qtitle,sizeof(qm_msg.libstr));
#endif

  /* set the name0,1 label length */
  if (m_msg.markx & MX_M10FORM) nml = 12;
  else nml = m_msg.nmlen;

  if (strlen(m_msg.qtitle) > 0) {
    if (m_msg.qtitle[0]=='>') strncpy(name0s,&m_msg.qtitle[1],sizeof(name0s));
    else strncpy(name0s,m_msg.qtitle,sizeof(name0s));
  }
  else {
    strncpy(name0s,m_msg.tname,sizeof(name0s));
  }
  name0s[sizeof(name0s)-1]='\0';

  if ((bp=strchr(name0s,' '))!=NULL) *bp='\0';

  if (m_msg.revcomp) name0[nml-1]='-';

  l_llen = m_msg.aln.llen;
  if ((m_msg.markx & MX_M9SUMM) && m_msg.show_code != SHOW_CODE_ID) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  sprintf(fmt,"%s%%-%ds (%%d %s)\n",A_MARK,l_llen-5,m_msg.sqnam);

  if (!(m_msg.markx&MX_M10FORM)) fprintf(fp,"\n");

  if (m_msg.ashow < 0) m_msg.ashow = m_msg.nshow;
  istart = 0; istop = min(min(nbest,m_msg.ashow),m_msg.nshow);

  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];

#ifdef SHOWUN
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msg.qsfnum,BBP_INFO(sfnum))) {
      istop = min(istop+1,nbest);
      continue;
    }
#endif
    if (bbp->score[0] <= 0) break;

    if (m_msg.quiet==1 && pst.zsflag>=0 
	&& bbp->escore < m_msg.e_low) continue;

#ifndef PCOMPLIB
    /* get the alignment and score by re-aligning */

    if ((m_fptr=re_openlib(bbp->m_file_p,!m_msg.quiet))==NULL)
      exit(1);

    /* get the description - do not "edit" it yet */

    if (!(m_msg.markx & MX_M10FORM)){
      if (m_msg.long_info) {tmp_len = sizeof(bline)-1;}
      else {tmp_len = l_llen-5;}
      RANLIB(bline,tmp_len,bbp->lseek,bbp->libstr,bbp->m_file_p);
      bline[tmp_len]='\0';
    }
    else {
      RANLIB(bline,sizeof(bline),bbp->lseek,bbp->libstr,bbp->m_file_p);
      bline[sizeof(bline)-1]='\0';
    }

    n1 = re_getlib(aa1,maxn,m_msg.maxt3,m_msg.loff,bbp->cont,m_msg.term_code,
		   &loffset,&l_off,bbp->m_file_p);
#ifdef DEBUG
    if (n1 != bbp->n1) {
      fprintf(stderr," library sequence: %s lengths differ: %d != %d\n",
	      bline,bbp->n1, n1);
      fprintf(stderr, "offset is: %lld\n",bbp->lseek);
    }
#endif

    if (!bbp->have_ares) {
      bbp->sw_score = 
	do_walign(aa0[bbp->frame],m_msg.n0, aa1, n1, bbp->frame, &pst,
		  f_str[bbp->frame], &bbp->a_res, &t_have_ares);
    }
    else {
      pre_cons(aa1,n1,bbp->frame,f_str[bbp->frame]);
    }

    aln_func_vals(bbp->frame, &m_msg.aln);

#else	/* PCOMPLIB - get the alignment information from a worker */

    /* we have a sequence that we need an alignment for -
       send a message to the appropriate worker to produce an alignment 
       qm_msg.slist == 1  -> one alignment
       qm_msg.s_func == DO_ALIGN_FLG -> use the alignment function
       send mngmsg (MSEQTYPE)
       then send number of sequence to be aligned
    */

    qm_msg.slist = 1;
    qm_msg.s_func = DO_ALIGN_FLG;

    liblist.seqnm = bbp->seqnm;
    liblist.frame = bbp->frame;
#ifdef PVM_SRC
    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)&qm_msg,sizeof(struct qmng_str),1);
    pvm_send(pinums[bbp->wrkr],MSEQTYPE);

    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)&liblist,sizeof(struct stage2_str),1);
    pvm_send(pinums[bbp->wrkr],LISTTYPE);
#endif
#ifdef MPI_SRC
    MPI_Send(&qm_msg,sizeof(struct qmng_str),MPI_BYTE,bbp->wrkr,
	     MSEQTYPE,MPI_COMM_WORLD);
    MPI_Send(&liblist,sizeof(struct stage2_str),MPI_BYTE,bbp->wrkr,
	     LISTTYPE,MPI_COMM_WORLD);
#endif
    /* information should be sent */
    /* pick up description */
    strncpy(bline,bbp->desptr->bline,l_llen-5);
    bline[l_llen-5]='\0';
#endif	/* PCOMPLIB */

    if (strlen(bline)==0) {
      bline[0]='>';
      strncpy(&bline[1],m_msg.lname,l_llen-5);
      bline[l_llen-5]='\0';
    }
    /* re-format bline */
    while ((bp=strchr(bline,'\n'))!=NULL) *bp=' ';
    if (m_msg.long_info) {
      tmp_len = strlen(bline);
      bl_ptr = bline;
      if (!(m_msg.markx & MX_M10FORM)) while (tmp_len > l_llen) {
	for (i=l_llen; i>10; i--)
	  if (bl_ptr[i]==' ') {
	    bl_ptr[i]='\n';
	    break;
	  }
	if (i <= 10) break;
	tmp_len -= i;
	bl_ptr += i;
      }
      bline[sizeof(bline)-1]='\0';
    }

    n1tot = (BBP_INFO(n1tot_p)) ? *BBP_INFO(n1tot_p) : bbp->n1;

    strncpy(name1,bline,sizeof(name1));

    if ((!m_msg.markx & MX_M10FORM)) name1[nml]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  Text must match that in showbest.c */

    strncpy(name1,bline,sizeof(name1));
    name1[sizeof(name1)-1]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';
    strncpy(l_name,name1,sizeof(l_name));
    l_name[sizeof(l_name)-1]='\0';
    if ((bp=strchr(&l_name[3],'|'))!=NULL) *bp='\0';
    if (m_msg.nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
    else if (m_msg.qframe >= 0 && bbp->frame == 1)
      strncat(l_name,"_r",sizeof(l_name));
    if (bbp->cont-1 > 0) {
      sprintf(tmp_str,":%d",bbp->cont-1);
      strncat(l_name,tmp_str,sizeof(l_name)-strlen(l_name));
    }

    if (!(m_msg.markx & MX_M10FORM)) name1[nml]='\0';

    /* print out score information; */

    if (m_msg.markx & MX_HTML ) {
      fprintf (fp,"<A name=%s>\n<tt><pre>\n",l_name);
    }
    strncpy(name0,name0s,nml);
    name0[nml]='\0';

    if (pst.zsflag%10 == 6) {
      sprintf(info_str," comp: %.5f H: %.5f",bbp->comp,bbp->H);
    }
    else info_str[0]='\0';

    if ((m_msg.markx & MX_ATYPE)!=7 && !(m_msg.markx & MX_M10FORM)) {
      fprintf (fp, fmt,bp=bline,n1tot);
      if (m_msg.nframe > 2) 
	fprintf (fp, "Frame: %d",bbp->frame+1);
      else if (m_msg.nframe > 1) 
	fprintf (fp, "Frame: %c",(bbp->frame? 'r': 'f'));
      else if (m_msg.qframe >= 0 && bbp->frame > 0 ) {
	  fputs("rev-comp",fp);
	  name0[nml-1]='\0';
	  strcat(name0,"-");
      }

      if (m_msg.arelv > 0)
	fprintf (fp, " %s: %3d", m_msg.alab[0],bbp->score[0]);
      if (m_msg.arelv > 1)
	fprintf (fp, " %s: %3d", m_msg.alab[1],bbp->score[1]);
      if (m_msg.arelv > 2)
	fprintf (fp, " %s: %3d", m_msg.alab[2],bbp->score[2]);
      fprintf(fp,"%s",info_str);
      if (pst.zsflag>=0) 
	fprintf (fp, "  Z-score: %4.1f  bits: %3.1f E(): %4.2g", 
		 bbp->zscore,zs_to_bit(bbp->zscore,m_msg.n0,bbp->n1),bbp->escore);
      fprintf (fp, "\n");
    }
    else if (m_msg.markx & MX_M10FORM) {
      fprintf(fp,">>%s\n",bline);
      if (m_msg.qframe > -1) {
	if (m_msg.nframe > 2) {
	  fprintf(fp,"; %s_frame: %d\n",m_msg.f_id0,bbp->frame+1);
	}
	else {
	  fprintf(fp,"; %s_frame: %c\n",m_msg.f_id0,(bbp->frame > 0? 'r':'f'));
	}
      }
      fprintf (fp, "; %s_%s: %3d\n", m_msg.f_id0,m_msg.alab[0],bbp->score[0]);
      if (m_msg.arelv > 1)
	fprintf (fp,"; %s_%s: %3d\n", m_msg.f_id0,m_msg.alab[1],bbp->score[1]);
      if (m_msg.arelv > 2)
	fprintf (fp,"; %s_%s: %3d\n", m_msg.f_id0,m_msg.alab[2],bbp->score[2]);
      if (info_str[0]) fprintf(fp,"; %s_info: %s\n",m_msg.f_id0,info_str);
      if (pst.zsflag>=0) 
     fprintf (fp,"; %s_z-score: %4.1f\n; %s_bits: %3.1f\n; %s_expect: %6.2g\n",
	      m_msg.f_id0,bbp->zscore,
	      m_msg.f_id0,zs_to_bit(bbp->zscore,m_msg.n0,bbp->n1),
	      m_msg.f_id0,bbp->escore);
    }
      

#ifdef PCOMPLIB
    /*  get the sw_score, alignment  information,  get seqc0, seqc1 */

#ifdef PVM_SRC
    /* get alignment lengths, percents */
    pvm_recv(pinums[bbp->wrkr],ALN1TYPE);
    pvm_upkint(&nc,1,1);
    pvm_upkint(&lc,1,1);
    pvm_upkint(&maxc,1,1);

    pvm_upkfloat(&percent,1,1);
    pvm_upkfloat(&gpercent,1,1);

    pvm_upkint(&bbp->sw_score,1,1);
    pvm_upkbyte((char *)&m_msg.aln,sizeof(struct a_struct),1);

    initseq(&seqc0, &seqc0a, &seqc1, &seqca, maxc);

    pvm_recv(pinums[bbp->wrkr],ALN2TYPE);
    pvm_upkbyte(seqc0,maxc,1);
    if (m_msg.ann_flg) pvm_upkbyte(seqc0a,maxc,1);
    pvm_upkbyte(seqc1,maxc,1);
    pvm_upkbyte(seqca,maxc,1);
#endif
#ifdef MPI_SRC
    MPI_Recv(int_msg_b,4,MPI_INT,bbp->wrkr,ALN1TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    nc = int_msg_b[0];
    lc = int_msg_b[1];
    maxc = int_msg_b[2];
    bbp->sw_score = int_msg_b[3];
    MPI_Recv(&percent,1,MPI_FLOAT,bbp->wrkr,ALN2TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(&gpercent,1,MPI_FLOAT,bbp->wrkr,ALN2TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(&m_msg.aln,sizeof(struct a_struct),MPI_BYTE,
	     bbp->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);

    initseq(&seqc0, &seqc0a, &seqc1, &seqca, maxc);
    MPI_Recv(seqc0,maxc,MPI_BYTE,bbp->wrkr,ALN2TYPE,MPI_COMM_WORLD,&mpi_status);
    if (m_msg.ann_flg)
      MPI_Recv(seqc0a,maxc,MPI_BYTE,bbp->wrkr,ALN2TYPE,MPI_COMM_WORLD,&mpi_status);
    MPI_Recv(seqc1,maxc,MPI_BYTE,bbp->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);
    MPI_Recv(seqca,maxc,MPI_BYTE,bbp->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);
#endif

    /* l_off is the coordinate of the first residue */
    l_off = 1;
    /* loffset is the offset of the aa1 in the full sequence */
    loffset = bbp->desptr->loffset-l_off;

#else	/* not PCOMPLIB */

    /* estimate space for alignment consensus */
    if (m_msg.aln.showall==1) {
      maxc = bbp->a_res.nres + max(bbp->a_res.min0,bbp->a_res.min1)+
	max((m_msg.n0-bbp->a_res.max0),(n1-bbp->a_res.max1))+4;
    }
    else {
      maxc = bbp->a_res.nres + 4*m_msg.aln.llen+4;
    }

    /* get space to put the sequence alignment consensus */
    initseq(&seqc0, &seqc0a, &seqc1, &seqca, maxc);

    /* build consensus from res, nres (done by workers if PCOMPLIB) */
    if (!m_msg.ann_flg) {
      nc=calcons(aa0[bbp->frame],m_msg.n0,aa1,n1,
		 &lc,&m_msg.aln, bbp->a_res, pst, seqc0, seqc1, seqca,
		 f_str[bbp->frame]);
      memset(seqc0a,' ',nc);
      seqc0a[nc]='\0';
    }
    else {
      nc=calcons_a(aa0[bbp->frame],m_msg.aa0a,m_msg.n0,aa1,n1,
		   &lc,&m_msg.aln,bbp->a_res,pst, seqc0, seqc0a, 
		   seqc1, seqca, m_msg.ann_arr,f_str[bbp->frame]);
    }

    /* PCOMPLIB workers return percent, gpercent, so calculate it here */
    if (lc > 0) percent = (100.0*(float)m_msg.aln.nident)/(float)lc;
    else percent = -1.00;
    ngap = m_msg.aln.ngap_q + m_msg.aln.ngap_l;
#ifndef SHOWSIM
    if (lc-ngap> 0) gpercent =(100.0*(float)m_msg.aln.nident)/(float)(lc-ngap);
#else
    if (lc > 0) gpercent =(100.0*(float)m_msg.aln.nsim)/(float)lc;
#endif
    else gpercent = -1.00;
#endif

    if (max(strlen(seqc0),strlen(seqc1)) > nc) {
      fprintf(stderr," mshowalign: nc/maxc: %d/%d seqc0/1: %u/%u\n",
	      nc,maxc,strlen(seqc0),strlen(seqc1));
    }

    /* here PCOMPLIB/comp_lib logic is the same */

#ifdef DEBUG
    if (bbp->sw_score < bbp->score[pst.score_ix]) {
      fprintf(stderr," *** warning - SW score=%d < opt score=%d ***\n",
	      bbp->sw_score, bbp->score[pst.score_ix]);
    }
#endif

    cal_coord(m_msg.n0,bbp->n1,m_msg.sq0off,loffset+l_off-1,&m_msg.aln);

#ifndef PCOMPLIB
    if (bbp->a_res.nres > 0)
#endif
      do_show(fp, m_msg.n0, bbp->n1, bbp->sw_score, name0, name1, nml,
	      m_msg, pst, seqc0, seqc0a, seqc1, seqca,
	      nc, percent, gpercent, lc, &m_msg.aln,
	      loffset+l_off-1);

    if (m_msg.markx & MX_HTML) fprintf(fp,"</pre></tt>\n<hr>\n");
    fflush(fp);

    freeseq(&seqc0,&seqc0a,&seqc1, &seqca);
  }
  if (fp!=stdout) fprintf(fp,"\n");
}

void do_show(FILE *fp, int n0,int n1, int score,
	     char *name0, char *name1, int nml,
	     struct mngmsg m_msg, struct pstruct pst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqca, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln, long loffset)
{
  int tmp;

  if (m_msg.markx & MX_AMAP && (m_msg.markx & MX_ATYPE)==7)
    disgraph(fp, n0, n1, percent, score,
	     aln->amin0, aln->amin1, aln->amax0, aln->amax1, m_msg.sq0off,
	     name0, name1, nml, aln->llen, m_msg.markx);
  else if (m_msg.markx & MX_M10FORM) {
    if (pst.sw_flag && m_msg.arelv>0)
      fprintf(fp,"; %s_score: %d\n",m_msg.f_id1,score);
    fprintf(fp,"; %s_ident: %5.3f\n",m_msg.f_id1,percent/100.0);
#ifndef SHOWSIM
    fprintf(fp,"; %s_gident: %5.3f\n",m_msg.f_id1,gpercent/100.0);
#else
    fprintf(fp,"; %s_sim: %5.3f\n",m_msg.f_id1,gpercent/100.0);
#endif

    fprintf(fp,"; %s_overlap: %d\n",m_msg.f_id1,lc);
    discons(fp, m_msg, pst, seqc0, seqc0a, seqc1, seqca, nc,
	    n0, n1, name0, name1, nml, aln, loffset);
  }
  else {
    if (pst.sw_flag) fprintf(fp,"Smith-Waterman score: %d; ",score);
    else fprintf(fp,"banded Smith-Waterman score: %d; ",score);
#ifndef SHOWSIM
    fprintf(fp," %6.3f%% identity (%6.3f%% ungapped) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	    percent,gpercent,lc,m_msg.sqnam,aln->d_start0,aln->d_stop0,
	    aln->d_start1,aln->d_stop1);
#else
    fprintf(fp," %6.3f%% identity (%6.3f%% similar) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	    percent,gpercent,lc,m_msg.sqnam,aln->d_start0,aln->d_stop0,
	    aln->d_start1,aln->d_stop1);
#endif

    if (m_msg.markx & MX_HTML) {
      do_url1(fp, m_msg, pst, l_name,n1,*aln,loffset);
    }

    if (m_msg.markx & MX_AMAP && (m_msg.markx & MX_ATYPE)!=7) {
      fputc('\n',fp);
      tmp = n0;

      if (m_msg.qdnaseq == SEQT_DNA && m_msg.ldnaseq== SEQT_PROT)
	tmp /= 3;

      disgraph(fp, tmp, n1, percent, score,
	       aln->amin0, aln->amin1,
	       aln->amax0, aln->amax1,
	       m_msg.sq0off,
	       name0, name1, nml, aln->llen,m_msg.markx);
    }

    discons(fp, m_msg, pst, seqc0, seqc0a, seqc1, seqca, nc,
	    n0, n1, name0, name1, nml, aln, loffset);

    fputc('\n',fp);

  }
}


#ifndef MPI_SRC
void	/* initialize consensus arrays */
initseq(char **seqc0, char **seqc0a, char **seqc1, char **seqca, int seqsiz)
{
  *seqc0=(char *)calloc((size_t)seqsiz*4,sizeof(char));
  if (*seqc0==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  *seqc0a=*seqc0 + seqsiz;
  *seqc1=*seqc0a + seqsiz;
  *seqca=*seqc1 + seqsiz;
}

void freeseq(char **seqc0, char **seqc0a, char **seqc1, char **seqca)
{
  free(*seqc0);
}
#endif
