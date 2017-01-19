
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa_34_26_5 $ - $Id: p2_workcomp.c,v 1.49 2007/01/02 17:24:36 wrp Exp $ */

/* This version is modifed to read all files, query and database,
   through the manager process. Workers will now receive their
   database from the manager, rather than reading it themselves.  This
   cuts down considerably on NFS traffic, simplifies searches of
   multiple files, and allows use of clusters of slave nodes that do
   not have NFS access */

/* September, 1994 - this version has been modified to do two kinds of
   searches, a general library search, or list of library sequences search.
   The latter would be used to generate optimized scores for fasta and
   to produce alignments */

/* modified July, 2002, to provide query shuffle */

/* modified October, 2005, to support struct a_res_str a_res -
   coordinates of alignment in aa0[], aa1[].  Future modifications
   will cause do_walign to be run only once - subsequent calls for
   seqc[0,1] can be filled using a_res, by adding a_res to the
   struct sqs2 array.

   19-March-2006 - modifications to call do_walign() only once, and
   use the resulting a_res structure for subsequent calls to calc_id,
   calcons, calcons_a, have been implemented.  Also, the -V option is
   now valid with the parallel programs.

   31-May-2006 - some functions (e.g. dropfs and dropff do not store
   complete information in a_res - thus they cannot use this shortcut
   (yet).

*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef PVM_SRC
#include "pvm3.h"
#endif

#ifdef MPI_SRC
#include "mpi.h"
#endif

/*
#define PvmDataDefault 0
#define PvmTaskDefault 0
*/
#include "msg.h"
#include "defs.h"
#include "param.h"
#include "w_mw.h"
#include "structs.h"

#ifdef MPI_SRC
#define XTERNAL
#endif
#include "upam.h"
#include "uascii.h"

#ifdef PVM_SRC
int worker, mytid;
int nnodes, pinums[MAXNOD];
#endif

#include "drop_func.h"

extern void alloc_pam (int d1, int d2, struct pstruct *ppst); /* allocate ppst->pam12,pam12x */
extern int **alloc_pam2p (int len, int nsq); 
extern void w_init ();
extern void irand(int);
extern void revcomp(unsigned char *, int, int *);



extern void initseq(char **seqc0, char **seqc0a, char **seqc1, char **seqca, int seqsiz);
extern void freeseq(char **seqc0, char **seqc0a, char **seqc1, char **seqca);

void send_bestr(int, int, struct comstr *, int, int);
void send_bestr2(int, struct comstr2 *, int);
void send_code(int, char *, int);

extern void get_param (struct pstruct *ppst, char *pstring2, char *pstring3);
extern void update_param(struct qmng_str *qm_msg, struct mngmsg *m_msg,
			 struct pstruct *ppst);
extern int shuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int, int *);

extern char err_str[];

/* local function declarations */
void free_ares(struct sqs2 *, int itt, int *, int walign_cnt, int worker);



void w_abort (p, p1)
char *p, *p1;
{
    fprintf (stderr, " %s %s\n", p, p1);
#ifdef PVM_SRC
    pvm_exit();
    exit (1);
#endif
#ifdef MPI_SRC
  MPI_Abort(MPI_COMM_WORLD,1);
#endif
}

#ifdef PVM_SRC
main ()
#endif
#ifdef MPI_SRC
void
workcomp(int worker)
#endif
{
  unsigned char *aa0[6], *aa1s, *aa0s;	/* Query and library sequences */
  struct mngmsg m_msg;	/* start message from manager to worker 1 */
  struct qmng_str qm_msg; /* updated for each query */
  int last_msg_b[10];	/* last set of numbers */
  struct sqs2 *seqpt; 		/* sequence pointers for chunk */
  int seqbuf_n,seqbuf_s;	/* number of sequences, length of sequences */
  int max_sql;			/* maximum number of sequences/node */
  int *n1_arr;			/* array of sequence lengths in buffer */
  int *m_seqnm_arr;		/* array of sequence numbers in buffer */
  int *aa1i_arr;		/* array of offsets into the buffer */
  unsigned char *seq_buf;	/* space for sequence data */
  int ntx;
  int nsq;			/* effective alphabet size */
  long curtype = ONETYPE;	/* current send message type */
  int ieven=0;			/* flag for window shuffle */
  int cur_n0;
  int n1, n1over;		/* length of query, library sequences */
  struct comstr bestr[BFR+1];	/* best structures */
  struct comstr2 bestr2[BFR2+1];	/* best structures */
  struct a_struct aln, *aln_dp;
  int qres_bufsize;		/* results buffer size */
  int bestcnt = 0;		/* how many best structures are full */
  char gstring2[MAX_STR];		/* parameter string for manager */
  char gstring3[MAX_STR];		/* parameter string for manager */
  struct pstruct pst;		/* parameter structure */
  struct rstruct rst, qrst, rrst;	/* results structure */
  void *f_str[6], *qf_str;
  int sw_score;
  int lcnt, count, seqnm;	/* counters */
  int *walign_done[2], walign_cnt[2];	/* index of current valid a_res in seqpt */
  int have_walign;
  int *tres;			/* allocated storage for seqpt[].a_res[].res */
  int lend;			/*global library sequence number information */
  int lsn;			/* library sequence number */
  struct stage2_str *liblist=NULL;	/* list of sequences to search */
  int i, j;		/* my turn to send sequence descriptions */
  char libstr[21];
  char errstr[128];
  int itt=0;
  int bufid;
  char *seqc0, *seqc0a, *seqc1, *seqca;
  char *seqc, *seqc_buff;
  int seqc_buff_cnt, seqc_buff_len, seqc_flag;
  int maxc, lc, nc, nident, ngap, aln_code_n;
  float percent, gpercent;
  int old_shuffle=0;	/* did a qshuffle last time */
  int hosttid=0;
  char worker_str[5];

#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif
  
#ifdef PVM_SRC
  mytid = pvm_mytid();
  hosttid = pvm_parent();
#endif

  w_init();	/* sets up default sascii, hsq, sq */

  /* Allocate space for the query sequence */
  if ((aa0[0] = (unsigned char *) malloc ((MAXTST+2+SEQ_PAD)*sizeof (char))) == NULL) {
    w_abort ("Unable to allocate sequence array[0] - exiting!","");
  }
  *aa0[0]='\0';
  aa0[0]++;

  /* initial messages set up various parameter structures:

     STARTTYPE0: &nnodes
     		 pinums
		 &m_msg

     STARTTYPE1  &pst

     STARTTYPE2	 pam12
     STARTTYPE3	 pam12x
  */

#ifdef PVM_SRC
#ifdef ROUTE_DIRECT
  pvm_setopt(PvmRoute,PvmRouteDirect);
#endif
  /* get number of nodes, pinums */
  bufid = pvm_recv(hosttid,STARTTYPE0);
  pvm_upkint(&nnodes,1,1);
  pvm_upkint(pinums,nnodes,1);
  pvm_upkbyte((char *)&m_msg,(int)sizeof(m_msg),1);
  worker = tidtonode(mytid);
  pvm_freebuf(bufid);
#endif

  sprintf(worker_str,"@%d",worker);

#ifdef MPI_SRC
  MPI_Recv(&m_msg,sizeof(m_msg),MPI_BYTE,hosttid,STARTTYPE0,MPI_COMM_WORLD,
	   &mpi_status);
#endif

  /* the aln structure needs some information from m_msg0.aln */
  memcpy(&aln,&m_msg.aln,sizeof(struct a_struct));

  /*
  fprintf(stderr,"d1: %d d2: %d\n",m_msg.pamd1,m_msg.pamd2);
  */

  /* get pst params */
#ifdef PVM_SRC
  bufid = pvm_recv(hosttid,STARTTYPE1);
  pvm_upkbyte((char *)&pst,(int)sizeof(pst),1);
  /* 31t nsq = pst.nsq; */
  pvm_upkbyte((char *)pascii,(int)sizeof(aascii),1);
  pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
  MPI_Recv(&pst,(int)sizeof(pst),MPI_BYTE,hosttid,STARTTYPE1,MPI_COMM_WORLD,
	   &mpi_status);

  MPI_Recv(pascii,(int)sizeof(aascii)/sizeof(int),MPI_INT,hosttid,STARTTYPE1,MPI_COMM_WORLD,
	   &mpi_status);
#endif

  if (pst.ext_sq_set) { nsq = pst.nsqx;}
  else { nsq = pst.nsq;}

  aa0[5] = aa0[4] = aa0[3] = aa0[2] = aa0[1] = aa0[0];
  if (m_msg.qframe == 2) {
    if ((aa0[1]=(unsigned char *)malloc((MAXTST+2)*sizeof (char)))==NULL)
      w_abort ("Unable to allocate sequence[1] array - exiting!","");
    *aa0[1]='\0';
    aa0[1]++;
  }

  if ((aa1s=(unsigned char *)malloc((m_msg.max_tot+1)*sizeof (char)))==NULL)
      w_abort ("Unable to allocate shuffled library sequence", "");
  *aa1s=0;
  aa1s++;

  irand(0);	/* necessary for shuffled sequences */

  /* this function allocates pam12, pam12x
     assigns pst.pam[0][0]=pam12, pst.pam[1][0] = pam12x
     and sets up the correct pst.pam[0][0][0] pointers */

  alloc_pam(m_msg.pamd1,m_msg.pamd2,&pst);

#ifdef PVM_SRC
  bufid = pvm_recv(hosttid,STARTTYPE2);
  pvm_upkint(pam12,m_msg.pamd1*m_msg.pamd2,1);
  pvm_freebuf(bufid);

  bufid = pvm_recv(hosttid,STARTTYPE3);
  pvm_upkint(pam12x,m_msg.pamd1*m_msg.pamd2,1);
  pvm_freebuf(bufid);
#endif  

#ifdef DEBUG
  if (worker==FIRSTNODE) {
    fprintf(stderr,"ext?: %d\tnsq: %d\tnsqx: %d\n",pst.ext_sq_set,pst.nsq, pst.nsqx);
    for (i=1; i<5; i++) {
      for (j=1; j <= i; j++) fprintf(stderr," %c,%c:%2d",pst.sq[i],pst.sq[j],pst.pam2[0][i][j]);
      fprintf(stderr,"\n");
    }
    for (i=pst.nsq+1; i<pst.nsq+5; i++) {
      for (j=pst.nsq+1; j <= i; j++) fprintf(stderr," %c,%c:%2d",pst.sqx[i],pst.sqx[j],pst.pam2[0][i][j]);
      fprintf(stderr,"\n");
    }

    for (i=1; i<5; i++) {
      for (j=1; j <= i; j++) fprintf(stderr," %c,%c:%2d",pst.sqx[i],pst.sqx[j],pst.pam2[1][i][j]);
      fprintf(stderr,"\n");
    }
    for (i=pst.nsq+1; i<pst.nsq+5; i++) {
      for (j=pst.nsq+1; j <= i; j++) fprintf(stderr," %c,%c:%2d",pst.sqx[i],pst.sqx[j],pst.pam2[1][i][j]);
      fprintf(stderr,"\n");
    }
  }
#endif

#ifdef MPI_SRC
  MPI_Recv(pam12,m_msg.pamd1*m_msg.pamd2,MPI_INT,hosttid,STARTTYPE2,
	   MPI_COMM_WORLD,&mpi_status);

  MPI_Recv(pam12x,m_msg.pamd1*m_msg.pamd2,MPI_INT,hosttid,STARTTYPE3,
	   MPI_COMM_WORLD,&mpi_status);
#endif

/*
  We have the PAM matrices - get the library sequences
*/

  /* Allocate space for the sequences */
  max_sql = MAXSQL/2;

  if ((seqpt=(struct sqs2 *)calloc(max_sql,sizeof(struct sqs2)))==NULL)
    w_abort("cannot allocate seqpt(sqs2)","");

  if ((n1_arr=(int *)calloc(m_msg.pbuf_siz+1,sizeof(int)))==NULL)
    w_abort("cannot allocate n1_arr","");

  if ((aa1i_arr=(int *)calloc(m_msg.pbuf_siz+1,sizeof(int)))==NULL)
    w_abort("cannot allocate n1_arr","");

  if ((m_seqnm_arr=(int *)calloc(m_msg.pbuf_siz+1,sizeof(int)))==NULL)
    w_abort("cannot allocate m_seqnm_arr","");

/*****************************************************************/
/* This section gets all the database sequences from the manager */
/*****************************************************************/

  lcnt = 0;
  while (1) {
#ifdef PVM_SRC
    /* get the number of sequences, sequence lengths */
    bufid = pvm_recv(hosttid,STARTTYPE4);
    pvm_upkint(&seqbuf_n,1,1);	/* number of sequences */
    pvm_upkint(&seqbuf_s,1,1);	/* size of sequence buffer */
    pvm_upkint(n1_arr,seqbuf_n,1);	/* length of each sequence in buffer */
    pvm_upkint(aa1i_arr,seqbuf_n,1);    /* indexes for each sequence */
    pvm_upkint(m_seqnm_arr,seqbuf_n,1);	/* number of each library sequence */
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv(&seqbuf_n,1,MPI_INT,hosttid,STARTTYPE4,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(&seqbuf_s,1,MPI_INT,hosttid,STARTTYPE4,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(n1_arr,seqbuf_n,MPI_INT,hosttid,STARTTYPE4,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(aa1i_arr,seqbuf_n,MPI_INT,hosttid,STARTTYPE4,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(m_seqnm_arr,seqbuf_n,MPI_INT,hosttid,STARTTYPE4,MPI_COMM_WORLD,
	     &mpi_status);
#endif

    if (seqbuf_n <= 0) break;
#ifdef DEBUG    
    /*
    fprintf(stderr,"[%d] seqbuf_n: %d seqbuf_s: %d\n",
	    worker,seqbuf_n,seqbuf_s);
    fprintf(stderr,"[%d] lcnt: %d n1: %d seqnm %d\n",
	    worker,0,n1_arr[0],m_seqnm_arr[0]);
    fprintf(stderr,"[%d] lcnt: %d n1: %d seqnm %d\n",
	    worker,1,n1_arr[1],m_seqnm_arr[1]);
    */
#endif

    /* allocate space for sequences */
    if ((seq_buf = (unsigned char *)calloc((size_t)seqbuf_s+1,sizeof(char)))
	==NULL) {
      w_abort("cannot allocate tmp_seq","");
    }
    seq_buf++; /* leave a '\0' at the start */

    /* get the sequence buffer */
#ifdef PVM_SRC
    bufid = pvm_recv(hosttid,STARTTYPE5);
    pvm_upkbyte((char *)seq_buf,seqbuf_s,1);
    pvm_freebuf(bufid);
#endif  
#ifdef MPI_SRC
    MPI_Recv(seq_buf,seqbuf_s,MPI_BYTE,hosttid,STARTTYPE5,MPI_COMM_WORLD,
	     &mpi_status);
#endif

    /* now we have everything  - update the pointers */
    if (lcnt+seqbuf_n >= max_sql) {
      max_sql += max(MAXSQL/2,seqbuf_n);
      if ((seqpt=(struct sqs2 *)realloc(seqpt,max_sql*sizeof(struct sqs2)))
	  ==NULL)
	w_abort("cannot allocate seqpt(sqs2)","");
    }

    /* convert from offsets to pointers into buffer */
    /* ntx = 0; */
    for (i=0; i<seqbuf_n; i++,lcnt++) {
      seqpt[lcnt].n1 = n1_arr[i];
      seqpt[lcnt].m_seqnm = m_seqnm_arr[i];
      seqpt[lcnt].aa1 = &seq_buf[aa1i_arr[i]];
      /*      ntx += n1_arr[i]+1 + SEQ_PAD */

#ifdef DEBUG
      /* must have null's at both ends of sequence */
      if (seqpt[lcnt].aa1[-1]!= '\0') {
	fprintf(stderr,"Missing null at start: %d %d\n",
		lcnt,seqpt[lcnt].aa1[-1]);
	seqpt[lcnt].aa1[-1]='\0';
      }
      if (seqpt[lcnt].aa1[seqpt[lcnt].n1]!= '\0') {
	fprintf(stderr,"Missing null at end: %d %d\n",
		lcnt,seqpt[lcnt].aa1[seqpt[lcnt].n1]);
	seqpt[lcnt].aa1[seqpt[lcnt].n1]='\0';
      }
#endif
    }
  }
  /* all done - lcnt has the total number of library sequences */

#ifdef DEBUG
  if (lcnt > 0)
    for (i=0; i<10; i++) {
      for (j=0; j<10; j++) libstr[j]=pst.sq[seqpt[i].aa1[j]];
      libstr[10]='\0';
      fprintf(stderr,"[%d] n1: %d seqnm: %d aa1: %s\n",
	      worker,seqpt[i].n1,seqpt[i].m_seqnm,libstr);
    }
#endif

  /* send back the number of descriptions received */

#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&lcnt,1,1);
  pvm_send(hosttid,STARTTYPE0);
#endif
#ifdef MPI_SRC
/*  p4_dprintf(" have %d descriptions to send\n",lcnt); */
  MPI_Send(&lcnt,1,MPI_INT,hosttid,STARTTYPE0,MPI_COMM_WORLD);
#endif  

/*****************************************************************/
/* Library reads are finished, get ready to do searches          */
/*****************************************************************/

  /* get last set of numbers */
#ifdef PVM_SRC
  bufid = pvm_recv(hosttid,STARTTYPE0);
  pvm_upkint(last_msg_b,2,1);
  pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
  MPI_Recv(last_msg_b, 2, MPI_INT, hosttid, STARTTYPE0, MPI_COMM_WORLD,
	     &mpi_status);
#endif
  m_msg.nbr_seq = last_msg_b[0];
  qres_bufsize = last_msg_b[1];

#ifdef DEBUG
#ifdef PVM_SRC
  fprintf(stderr,"[%d] have nbr_seq %d qres_bufsize %d\n",worker,
	     m_msg.nbr_seq, qres_bufsize);
#endif
#ifdef MPI_SRC
  /*  p4_dprintf("[%d] have nbr_seq %d qres_bufsize %d\n",worker,
	     m_msg.nbr_seq, qres_bufsize);
  */;
#endif
#endif  
  /* If self search, receive sequence numbering data */
  if (m_msg.self) {
#ifdef PVM_SRC
    bufid = pvm_recv(hosttid,STARTTYPE1);
    pvm_upkint(&lend,1,1);
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv(&lend,1,MPI_INT,hosttid,STARTTYPE1,MPI_COMM_WORLD,&mpi_status);
#endif
  }
  
  /* allocate space for a_res flag array */

  if ((walign_done[0] = (int *)calloc(lcnt,sizeof(int)))==NULL) {
    w_abort("cannot allocate walign_done");
  }
  walign_cnt[0]=0;

  if ((walign_done[1] = (int *)calloc(lcnt,sizeof(int)))==NULL) {
    w_abort("cannot allocate walign_done");
  }
  walign_cnt[1]=0;

  /* was commented in for only FASTX/TFASTX, but do it always to
     simplify */
  aainit(pst.tr_type, pst.debug_lib);
  pst.maxlen = m_msg.maxn;

/*****************************************************************/
/* Main search loop, which calles do_work() repeatedly           */
/*****************************************************************/

  cur_n0 = 0;
  while (1)  {
/*
#ifdef DEBUG
#ifdef PVM_SRC
    fprintf(stderr," W: %d waiting MSEQTYPE\n",worker);
#endif
#ifdef MPI_SRC
    p4_dprintf(" W: %d waiting MSEQTYPE\n",worker);
#endif
#endif
*/

/*****************************************************************/
/* Wait for a query sequence from the manager                    */
/*****************************************************************/

#ifdef PVM_SRC
    bufid = pvm_recv(hosttid,MSEQTYPE);
    pvm_upkbyte((char *)&qm_msg,sizeof(qm_msg),1);
#endif
#ifdef MPI_SRC
    MPI_Recv(&qm_msg,sizeof(struct mngmsg),MPI_BYTE,hosttid,MSEQTYPE,
	     MPI_COMM_WORLD,&mpi_status);
#endif
#ifdef DEBUG
    fprintf(stderr,"[%d] have MSEQTYPE n0: %d s_func: %d slist: %d qf: %d\n",
	    worker,qm_msg.n0,qm_msg.s_func,qm_msg.slist,qm_msg.qshuffle);
#endif

/*****************************************************************/
/* New query sequence indicated by qm_msg.slist=0                */
/*****************************************************************/

    if (qm_msg.n0 > 0 && qm_msg.slist == 0) {

      if (cur_n0 > 0) {

/*****************************************************************/
/*    free everything associated with previous search            */
/*****************************************************************/

 	close_work (aa0[0], cur_n0, &pst, &f_str[0]);  
	free_ares(seqpt, 0, walign_done[0], walign_cnt[0], worker);
	walign_cnt[0] = 0;
	if (m_msg.ann_flg) free(m_msg.aa0a);


	if (m_msg.qframe == 2) {
	  close_work(aa0[1], cur_n0, &pst, &f_str[1]);
	  free_ares(seqpt, 1, walign_done[1], walign_cnt[1], worker);
	  walign_cnt[1] = 0;
	}
	if (old_shuffle) {
	  close_work(aa0s,cur_n0, &pst, &qf_str);
	  aa0s--;
	  free(aa0s);
	  old_shuffle = 0;
	}
	if (pst.pam_pssm) {
	  free_pam2p(pst.pam2p[0]);
	  free_pam2p(pst.pam2p[1]);
	}
      }

/*****************************************************************/
/*    Start allocating things for the next search                */
/*****************************************************************/

      pst.pam_pssm = qm_msg.pam_pssm;
      cur_n0 = qm_msg.n0;
      if (m_msg.ann_flg) {
	if ((m_msg.aa0a = calloc(qm_msg.n0+1,sizeof(char)))==NULL) {
	  w_abort(" cannot allocate aa0a");
	}
      }

/*****************************************************************/
/*    Get the next query sequence                                */
/*****************************************************************/

#ifdef PVM_SRC
      pvm_upkbyte((char *)aa0[0],qm_msg.n0+1+SEQ_PAD,1);
      if (m_msg.ann_flg) {
	pvm_upkbyte((char *)m_msg.aa0a,qm_msg.n0+1,1);
      }
#endif
#ifdef MPI_SRC
      MPI_Recv(aa0[0],qm_msg.n0+1+SEQ_PAD,MPI_BYTE,hosttid,
	       MSEQTYPE1,MPI_COMM_WORLD, &mpi_status);
      if (m_msg.ann_flg) {
	MPI_Recv(m_msg.aa0a,qm_msg.n0+1,MPI_BYTE,hosttid,
		 MSEQTYPE2,MPI_COMM_WORLD, &mpi_status);
      }
#endif

#ifdef DEBUG
      /* must have null's at both ends of sequence */
      if (aa0[0][-1]!= '\0') {
	fprintf(stderr,"Missing null at start: %s %d\n",
		qm_msg.libstr,aa0[0][-1]);
	aa0[0][-1]='\0';
      }
      if (aa0[0][qm_msg.n0]!= '\0') {
	fprintf(stderr,"Missing null at end: %s %d\n",
		qm_msg.libstr,aa0[0][qm_msg.n0]);
	aa0[qm_msg.n0]='\0';
      }

      /* This discovers most reasons for core dumps */
      if (pst.debug_lib)
	for (j=0; j<qm_msg.n0; j++)
	  if (aa0[0][j]>pst.nsq) {
	    fprintf(stderr,
		    "seq: %s residue[%d/%d] %d range (%d)\n",
		    qm_msg.libstr,j,qm_msg.n0,aa0[0][j],pst.nsq);
	    aa0[0][j]=0;
	    qm_msg.n0=j-1;
	    break;
	  }
#endif
      update_params(&qm_msg,&m_msg,&pst);
    }

/*****************************************************************/
/*    End of free()'s/ initialization for new sequence           */
/*****************************************************************/

#ifdef PVM_SRC
    pvm_freebuf(bufid);
#endif

    if (qm_msg.n0 == -1) {

/*****************************************************************/
/*   All done with searches                                      */
/*****************************************************************/
/*    printf(" %d: got n0 == -1\n",worker); */
      break;
    }

    /*    p4_dprintf(" W:%d n0:%d slist:%d s_func:%d (%d)\n",worker,qm_msg.n0,qm_msg.slist,qm_msg.s_func,qres_bufsize); */

/*****************************************************************/
/*   if qm_msg.slist > 0, search specific sequences, to be sent  */
/*****************************************************************/

    if (qm_msg.slist > 0) {  /* list search, not library search */
      if (liblist != NULL) free(liblist);

      /* get the list of sequences */
      if ((liblist=(struct stage2_str *)
	   calloc(qm_msg.slist,sizeof(struct stage2_str)))==NULL) {
	  sprintf(errstr,"sequence list %d",qm_msg.slist);
	  w_abort (errstr, "");
	}

#ifdef PVM_SRC
      bufid = pvm_recv(hosttid,LISTTYPE);
      pvm_upkbyte((char *)liblist,qm_msg.slist*sizeof(struct stage2_str),1);
      pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
      MPI_Recv(liblist,qm_msg.slist*sizeof(struct stage2_str),MPI_BYTE,
	       hosttid,LISTTYPE,MPI_COMM_WORLD, &mpi_status);
#endif
    }

/*****************************************************************/
/* have list of sequences to be compared/aligned                 */
/*****************************************************************/

    /* Initial stuff */
    if (qm_msg.slist == 0) {
/*****************************************************************/
/* New query - set up matrices and init_work()                   */
/*****************************************************************/
#ifdef DEBUG
/*
      fprintf(stderr,"n1: %d\t",qm_msg.n0);
      for (i=0; i<10; i++) fprintf(stderr,"%c",nt[aa0[0][i]]);
      fprintf(stderr,"\n");
*/
#endif
      if (pst.pam_pssm) {
	pst.pam2p[0] = alloc_pam2p(qm_msg.n0,nsq);
	pst.pam2p[1] = alloc_pam2p(qm_msg.n0,nsq);
      }

      init_work (aa0[0], qm_msg.n0, &pst, &f_str[0]);
      f_str[5]=f_str[4]=f_str[3]=f_str[2]=f_str[1]=f_str[0];

      if (qm_msg.qshuffle) {
	if ((aa0s=(unsigned char *)malloc((qm_msg.n0+2)*sizeof (char)))==NULL)
	  w_abort ("Unable to allocate aa0s array - exiting!","");
	*aa0s='\0';
	aa0s++;

	memcpy(aa0s,aa0[0],qm_msg.n0+1);
	qshuffle(aa0s,qm_msg.n0,qm_msg.nm0);
#ifdef DEBUG
	fprintf(stderr,"[%d] shuffle: %d\n",worker,qm_msg.n0);
	fputs("   ",stderr);
	for (i=0; i<5; i++) {fprintf(stderr,"%c",pst.sq[aa0s[i]]);}
	fputc('\n',stderr);
#endif

	init_work (aa0s, qm_msg.n0, &pst, &qf_str);
	old_shuffle=1;
      }

      if (m_msg.qframe == 2) {
	memcpy(aa0[1],aa0[0],qm_msg.n0+1);
	revcomp(aa0[1],qm_msg.n0,&pst.c_nt[0]);
	init_work (aa0[1], qm_msg.n0, &pst, &f_str[1]);
      }
#ifdef DEBUG
/*
	fprintf(stderr,"[%d] init_work qf: %d nf: %d\n",worker,m_msg.qframe,m_msg.nframe);
*/
#endif
    }

/*****************************************************************/
/* Finished with initialization,                                 */
/* start doing comparisons or alignments                         */
/*****************************************************************/

    bestcnt = 0;
    if (qm_msg.slist == 0) {	/* library search */

/*****************************************************************/
/* Start library search                                          */
/*****************************************************************/

      for (count=0; count < lcnt; count++) {

	for (itt=m_msg.revcomp; itt<=m_msg.nitt1; itt++) {

	  rst.score[0] = rst.score[1] = rst.score[2] = 0;
	  if (m_msg.self) {
	    lsn = lend + count;
	    if ((qm_msg.seqnm > lsn) && (((qm_msg.seqnm + lsn) % 2) != 0)) {
	      do_work (aa0[itt], qm_msg.n0,seqpt[count].aa1, seqpt[count].n1,
		       itt, &pst, f_str[itt], 0, &rst);
	    }
	    else if ((qm_msg.seqnm <= lsn) && (((qm_msg.seqnm+lsn)%2) == 0)) {
	      do_work (aa0[itt], qm_msg.n0, seqpt[count].aa1, seqpt[count].n1, 
		       itt, &pst, f_str[itt], 0, &rst);
	    }
	    else continue;
	  }
	  else {
	    do_work (aa0[itt], qm_msg.n0, seqpt[count].aa1, seqpt[count].n1,
		     itt, &pst, f_str[itt], 0, &rst);
	    if (qm_msg.qshuffle) {
	      do_work (aa0s, qm_msg.n0, seqpt[count].aa1, seqpt[count].n1,
		     itt, &pst, qf_str, 1, &qrst);
	    }
	  }
#ifdef DEBUG
/*
	  if (count < 10 || (count % 200 == 199)) {
	    fprintf(stderr,"[node %d] itt:%d/%d (%d) %3d %3d %3d - %d/%d\n",
		    worker,itt,m_msg.nitt1,count,
		    rst.score[0],rst.score[1],rst.score[2],
		    seqpt[count].m_seqnm,seqpt[count].n1);
	  }
*/
#endif
	  sw_score = -1;

	  bestr[bestcnt].seqnm  = count;
	  bestr[bestcnt].m_seqnm  = seqpt[count].m_seqnm;
	  bestr[bestcnt].score[0] = rst.score[0];
	  bestr[bestcnt].score[1] = rst.score[1];
	  bestr[bestcnt].score[2] = rst.score[2];
	  bestr[bestcnt].escore = rst.escore;
	  bestr[bestcnt].segnum = rst.segnum;
	  bestr[bestcnt].seglen = rst.seglen;
	  bestr[bestcnt].frame = itt;
	  bestr[bestcnt].comp = rst.comp;
	  bestr[bestcnt].H = rst.H;

	  bestr[bestcnt].qr_score = qrst.score[pst.score_ix];
	  bestr[bestcnt].qr_escore = qrst.escore;

	  if (pst.zsflag >= 10) {
	    if (pst.zs_win > 0) 
	      wshuffle(seqpt[count].aa1, aa1s,seqpt[count].n1,pst.zs_win,&ieven);
	    else 
	      shuffle(seqpt[count].aa1, aa1s,seqpt[count].n1);

	    do_work(aa0[itt],qm_msg.n0,aa1s,seqpt[count].n1,itt, &pst, 
		    f_str[itt], 0, &rst);
	    bestr[bestcnt].r_score = rst.score[pst.score_ix];
	  }

	  bestcnt++;
	  if (bestcnt >= qres_bufsize) {
#ifdef DEBUG
	    fprintf(stderr," worker: %d sending %d results\n",worker,qres_bufsize);
#endif
	    send_bestr(hosttid,curtype,bestr,qres_bufsize,bestcnt);
	    bestcnt = 0;
	  }
	}
      }	/* END - for count loop */
      send_bestr(hosttid, curtype, bestr,qres_bufsize, (bestcnt | FINISHED));
    }

/*****************************************************************/
/* End of library search section                                 */
/*****************************************************************/

/*****************************************************************/
/* Do do_opt() from list   s_func=DO_CALC_FLG                    */
/*****************************************************************/

    else if (qm_msg.s_func== DO_CALC_FLG) {  /* qm_msg.slist > 0 */

      bestcnt = 0;
      for (count=0; count < qm_msg.slist; count++) {
	rst.score[0] = rst.score[1] = rst.score[2] = 0;
	itt = liblist[count].frame;
	seqnm = bestr2[bestcnt].seqnm  = liblist[count].seqnm;
	bestr2[bestcnt].m_seqnm = seqpt[seqnm].m_seqnm;

	do_opt (aa0[itt], qm_msg.n0, seqpt[seqnm].aa1,
		 seqpt[seqnm].n1, itt,
		 &pst, f_str[itt], &rst);

	bestr2[bestcnt].score[0] = rst.score[0];
	bestr2[bestcnt].score[1] = rst.score[1];
	bestr2[bestcnt].score[2] = rst.score[2];
	bestr2[bestcnt].escore = rst.escore;
	bestr2[bestcnt].segnum = rst.segnum;
	bestr2[bestcnt].seglen = rst.seglen;
	bestr2[bestcnt].aln_code_n = 0;
	bestcnt++;

	if (bestcnt >= BFR2) {
	  send_bestr2(hosttid,bestr2,bestcnt);
	  bestcnt = 0;
	}
      }	/* END - for count loop */

      send_bestr2(hosttid,bestr2,(bestcnt|FINISHED));
    }

/*****************************************************************/
/* s_func=DO_OPT_FLG                                             */
/*                                                               */
/* from list:                                                    */
/* if (m_msg.stages > 1) do_opt()                                */
/* do_walign()                                                   */
/* calc_id or calc_code, no calcons                              */
/*****************************************************************/

    /* s_func == 1 means do_opt if necessary */
    else if (qm_msg.s_func== DO_OPT_FLG) {  /* qm_msg.slist > 0 */
#ifdef DEBUG
      fprintf(stderr," [%d] starting s_func:1 slist: %d\n",
	      worker,qm_msg.slist);
#endif
      /* get the buffer once - re-use it for the entire slist */
      if (m_msg.show_code == SHOW_CODE_ALIGN) {
	seqc_buff_len = (BFR2+5)*256;
	seqc = seqc_buff = (char *)calloc(seqc_buff_len,sizeof(char));
	seqc_buff_cnt = 0;
	if (seqc_buff == NULL) {
	  seqc_buff_cnt = seqc_buff_len = 0;
	}
      }

      bestcnt = 0;
      for (count=0; count < qm_msg.slist; count++) {
	rst.score[0] = rst.score[1] = rst.score[2] = 0;
	itt = liblist[count].frame;
	seqnm = liblist[count].seqnm;

	bestr2[bestcnt].seqnm  = seqnm;
	bestr2[bestcnt].m_seqnm = seqpt[seqnm].m_seqnm;
	if (m_msg.stages > 1) {
	  do_opt (aa0[itt], qm_msg.n0, seqpt[seqnm].aa1,
		  seqpt[seqnm].n1, itt,
		  &pst, f_str[itt], &rst);

	  bestr2[bestcnt].score[0] = rst.score[0];
	  bestr2[bestcnt].score[1] = rst.score[1];
	  bestr2[bestcnt].score[2] = rst.score[2];
	}

	if (m_msg.markx & MX_M9SUMM) {
#ifdef DEBUG
	  fprintf(stderr," [%d] starting do_walign seqnm: %d n1: %d\n",
		  worker,seqnm,seqpt[seqnm].n1);
#endif
	  aln_dp = &bestr2[bestcnt].aln_d;
	  memcpy(aln_dp, &aln,sizeof(struct a_struct));

	  sw_score = do_walign(aa0[itt], qm_msg.n0,
			       seqpt[seqnm].aa1, seqpt[seqnm].n1,
			       itt, &pst, f_str[itt],
			       &seqpt[seqnm].a_res[itt],
			       &have_walign);
	  seqpt[seqnm].sw_score[itt] = sw_score;

	  /* the a_res[itt] provided by do_walign is re-used - so it
	     must be copied to a valid location */

	  if (have_walign) {
	    if ((tres = calloc(seqpt[seqnm].a_res[itt].nres+1,sizeof(int)))==NULL) {
	      w_abort(" cannot allocate tres");
	    }
	    else {
	      memcpy(tres,seqpt[seqnm].a_res[itt].res,sizeof(int)*seqpt[seqnm].a_res[itt].nres);
	      seqpt[seqnm].a_res[itt].res = tres;
	      /*
		fprintf(stderr, " [%d] saving %d:%d[%d]:%o\n", worker, 
		walign_cnt[itt],seqnm,itt, seqpt[seqnm].a_res[itt].res);
	      */
	      if (walign_cnt[itt] < lcnt) walign_done[itt][walign_cnt[itt]++] = seqnm;
	      else w_abort(" walign_cnt overrun");
	      seqpt[seqnm].walign_dflg[itt] = 1;
	    }
	  }
	  aln_func_vals(itt, aln_dp);

#ifdef DEBUG
	  fprintf(stderr," [%d] starting calc_id sw_score: %d\n",
		  worker,sw_score);
	  fprintf(stderr,"bi: %d seqc_buff_cnt: %d - seqc_buff_len: %d\n",
		  bestcnt, seqc_buff_cnt, seqc_buff_len);
#endif
	  aln_code_n = 0;	/* must be set in case no seqc_code */
	  if (m_msg.show_code == SHOW_CODE_ALIGN) {
	    if (seqc_buff_cnt < seqc_buff_len - 256) {
	      lc=calc_code(aa0[itt],qm_msg.n0,
			   seqpt[seqnm].aa1, seqpt[seqnm].n1,
			   aln_dp,seqpt[seqnm].a_res[itt],pst,
			   seqc,seqc_buff_len-seqc_buff_cnt-10,
			   f_str[itt]);
	      aln_code_n = strlen(seqc);
	      seqc_buff_cnt += aln_code_n + 1;
/*
	      fprintf(stderr,"%d:%d:%d: %d/%d - [%d] %s\n",
		      worker,seqnm,bestcnt,aln_code_n,seqc_buff_cnt, seqc-seqc_buff,seqc);
*/
	      seqc += aln_code_n;
	      *seqc++ = '\0';
	    }
	  }
	  else {
	    lc=calc_id(aa0[itt],qm_msg.n0,
		       seqpt[seqnm].aa1, seqpt[seqnm].n1,
		       aln_dp,seqpt[seqnm].a_res[itt],pst,f_str[itt]);
	  }

	  nident = aln_dp->nident;
	  aln_dp->a_len = lc;

	  if (lc > 0) percent = (100.0*(float)nident)/(float)lc;
	  else percent = 0.0;

	  ngap = aln_dp->ngap_q + aln_dp->ngap_l;
#ifndef SHOWSIM
	  if (lc-ngap > 0) gpercent = (100.0*(float)nident)/(float)(lc-ngap);
#else
	  if (lc > 0) gpercent =(100.0*(float)aln_dp->nsim)/(float)lc;
#endif
	  else gpercent = -1.0;

	  bestr2[bestcnt].sw_score = sw_score;
	  bestr2[bestcnt].percent = percent;
	  bestr2[bestcnt].gpercent = gpercent;
	  bestr2[bestcnt].aln_code_n = aln_code_n;
	}
	bestcnt++;

	if (bestcnt >= BFR2) {
	  send_bestr2(hosttid,bestr2,bestcnt);
	  if (m_msg.show_code == SHOW_CODE_ALIGN) {
	    send_code(hosttid,seqc_buff,seqc_buff_cnt);
	    memset(seqc_buff,0,seqc_buff_len);
	    seqc = seqc_buff;
	    seqc_buff_cnt = 0;
	  }
	  bestcnt = 0;
	}
      }	/* END - for count loop */

      send_bestr2(hosttid,bestr2,(bestcnt|FINISHED));
      if (m_msg.show_code == SHOW_CODE_ALIGN) {
	send_code(hosttid,seqc_buff,seqc_buff_cnt);
	if (seqc_buff) free(seqc_buff);
      }
    }
    /* get alignments */

/*****************************************************************/
/* s_list >                                                      */
/* s_func=DO_ALIGN_FLG                                           */
/*                                                               */
/* from list:                                                    */
/* do_walign() if not done already                               */ 
/* calcons()                                                     */
/*****************************************************************/

    else if (qm_msg.s_func==DO_ALIGN_FLG) {
      for (count=0; count < qm_msg.slist; count++) {
	itt = liblist[count].frame;
	seqnm = liblist[count].seqnm;
/*
	fprintf(stderr,"worker: %d; %s, frame: %d\n",worker,qm_msg.libstr,itt);
*/
	if (!seqpt[seqnm].walign_dflg[itt]) {
	  seqpt[seqnm].sw_score[itt] = 
	    sw_score = do_walign (aa0[itt], qm_msg.n0,seqpt[seqnm].aa1,
				  seqpt[seqnm].n1, itt,
				  &pst, f_str[itt],
				  &seqpt[seqnm].a_res[itt],
				  &have_walign);
	}
	else {
	  sw_score = seqpt[seqnm].sw_score[itt];
	  pre_cons(seqpt[seqnm].aa1,seqpt[seqnm].n1,itt,f_str[itt]);
	}

	aln_func_vals(itt, &aln);

	if (aln.showall==1)
	  maxc = seqpt[seqnm].a_res[itt].nres + max(seqpt[seqnm].a_res[itt].min0,seqpt[seqnm].a_res[itt].min1)+
	    max((qm_msg.n0-seqpt[seqnm].a_res[itt].max0),
		(seqpt[seqnm].n1-seqpt[seqnm].a_res[itt].max1))+4;
	else  maxc = seqpt[seqnm].a_res[itt].nres + 4*aln.llen+4;

	initseq(&seqc0, &seqc0a, &seqc1, &seqca, maxc);

	if (!m_msg.ann_flg) {
	  nc=calcons(aa0[itt],qm_msg.n0,
		   seqpt[seqnm].aa1, seqpt[seqnm].n1,
		   &lc,&aln,seqpt[seqnm].a_res[itt],pst,
		   seqc0,seqc1,seqca,f_str[itt]);
	  memset(seqc0a,' ',nc);
	  seqc0a[nc]='\0';
	}
	else {
	  nc=calcons_a(aa0[itt],m_msg.aa0a,qm_msg.n0,
		       seqpt[seqnm].aa1, seqpt[seqnm].n1,
		       &lc,&aln,seqpt[seqnm].a_res[itt],pst,
		       seqc0,seqc0a,seqc1,seqca,
		       m_msg.ann_arr,f_str[itt]);
	}

	/*
	fprintf(stderr,"[%d] nident: %d nsim: %d lc: %d\n",aln.nident, aln.nsim, lc);
	*/

	maxc = max(strlen(seqc0),strlen(seqc1))+1;
	nident = aln.nident;
	percent = (100.0*(float)nident)/(float)lc;
	ngap = aln.ngap_q+aln.ngap_l;
#ifndef SHOWSIM
	if (lc-ngap > 0) gpercent = (100.0*(float)nident)/(float)(lc-ngap);
#else
	if (lc > 0) gpercent = (100.0*(float)aln.nsim)/(float)lc;
#endif
	else gpercent = -1.0;

#ifdef PVM_SRC
	pvm_initsend(PvmDataRaw);
	pvm_pkint(&nc,1,1);
	pvm_pkint(&lc,1,1);
	pvm_pkint(&maxc,1,1);
	pvm_pkfloat(&percent,1,1);
	pvm_pkfloat(&gpercent,1,1);
	pvm_pkint(&sw_score,1,1);
	pvm_pkbyte((char *)&aln,sizeof(struct a_struct),1);
	pvm_send(hosttid,ALN1TYPE);
#ifdef DEBUG
	fprintf(stderr,"[%d] ALN1TYPE sent: %d\n",worker,qm_msg.n0);
#endif
	pvm_initsend(PvmDataRaw);
	pvm_pkbyte(seqc0,maxc,1);
	if (m_msg.ann_flg) pvm_pkbyte(seqc0a,maxc,1);
	pvm_pkbyte(seqc1,maxc,1);
	pvm_pkbyte(seqca,maxc,1);
	pvm_send(hosttid,ALN2TYPE);
#endif
#ifdef MPI_SRC
	last_msg_b[0]=nc;
	last_msg_b[1]=lc;
	last_msg_b[2]=maxc;
	last_msg_b[3]=sw_score;
	MPI_Send(last_msg_b,4,MPI_INT,hosttid,ALN1TYPE,MPI_COMM_WORLD);
	MPI_Send(&percent,1,MPI_FLOAT,hosttid,ALN2TYPE,MPI_COMM_WORLD);
	MPI_Send(&gpercent,1,MPI_FLOAT,hosttid,ALN2TYPE,MPI_COMM_WORLD);

/*	p4_dprintf("[%d] sending aln\n",worker); */
	MPI_Send(&aln,sizeof(struct a_struct),MPI_BYTE,hosttid,
		 ALN3TYPE,MPI_COMM_WORLD);

	MPI_Send(seqc0,maxc,MPI_BYTE,hosttid,ALN2TYPE,MPI_COMM_WORLD);
	if (m_msg.ann_flg) MPI_Send(seqc0a,maxc,MPI_BYTE,hosttid,ALN2TYPE,MPI_COMM_WORLD);
	MPI_Send(seqc1,maxc,MPI_BYTE,hosttid,ALN3TYPE,MPI_COMM_WORLD);
	MPI_Send(seqca,maxc,MPI_BYTE,hosttid,ALN3TYPE,MPI_COMM_WORLD);
#endif
	freeseq(&seqc0,&seqc0a,&seqc1,&seqca);
      }	
    }
    
/* send back parameter settings */
    if (worker==FIRSTWORK && qm_msg.slist==0) {
      get_param(&pst, gstring2,gstring3);
#ifdef PVM_SRC
      pvm_initsend(PvmDataRaw);
      pvm_pkbyte(gstring2,sizeof(gstring2),1);
      pvm_pkbyte(gstring3,sizeof(gstring3),1);
      pvm_send(hosttid,PARAMTYPE);
#endif
#ifdef MPI_SRC
      MPI_Send(gstring2,sizeof(gstring2),MPI_BYTE,
	       hosttid,PARAMTYPE,MPI_COMM_WORLD);
      MPI_Send(gstring3,sizeof(gstring3),MPI_BYTE,
	       hosttid,PARAMTYPE,MPI_COMM_WORLD);
#endif
    }
    
    if (qm_msg.slist==0) {
      if (curtype == ONETYPE) curtype = TWOTYPE;
      else curtype = ONETYPE;
    }
  }	    /* END - while (1) loop */
#ifdef PVM_SRC
  pvm_exit();
#endif
#ifdef MPI_SRC
/*  MPI_Finalize(); */
#endif
}

void
send_bestr(int hosttid, int curtype, 
	   struct comstr *bestr, int buf_size, int lastcnt) {

  bestr[buf_size].seqnm = lastcnt;
#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkbyte((char *)&bestr[0],sizeof(struct comstr)*(buf_size+1),1);
  pvm_send(hosttid,curtype);
#endif
#ifdef MPI_SRC
  MPI_Send(bestr,sizeof(struct comstr)*(buf_size+1),MPI_BYTE,
	   hosttid,curtype,MPI_COMM_WORLD);
#endif
}

void
send_bestr2(int hosttid, struct comstr2 *bestr2,
	    int lastcnt)
{
  bestr2[BFR2].seqnm = lastcnt;
#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkbyte((char *)&bestr2[0],sizeof(struct comstr2)*(BFR2+1),1);
  pvm_send(hosttid,LISTRTYPE);
#endif
#ifdef MPI_SRC
  MPI_Send(&bestr2[0],sizeof(struct comstr2)*(BFR2+1),MPI_BYTE,
	   hosttid,LISTRTYPE,MPI_COMM_WORLD);
#endif
}

void
send_code(int hosttid, char *seqc_buff, int seqc_buff_len) {

#ifdef PVM_SRC
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&seqc_buff_len,1,1);
  if (seqc_buff_len > 0) pvm_pkbyte(seqc_buff,seqc_buff_len,1);
  pvm_send(hosttid,CODERTYPE);
#endif
#ifdef MPI_SRC
  MPI_Send(&seqc_buff_len,1,MPI_INT,
	   hosttid,CODERTYPE,MPI_COMM_WORLD);
  if (seqc_buff_len>0) MPI_Send(seqc_buff,seqc_buff_len,MPI_BYTE,
			       hosttid,CODERTYPE,MPI_COMM_WORLD);
#endif
}

#ifdef PVM_SRC
int tidtonode(tid)
     int tid;
{
  int i;
  for (i=FIRSTNODE; i< nnodes; i++) if (tid==pinums[i]) return i;
  fprintf(stderr," cannot find tid %d\n",tid);
  return -1;
}
#endif

void
free_ares(struct sqs2 *seqpt, int itt, int *walign_done, int walign_cnt, int worker) {

  int i, seqnm;

  for (i=0; i< walign_cnt; i++) {
    seqnm = walign_done[i];
    walign_done[i]=0;
    if (seqpt[seqnm].walign_dflg[itt]) {
      if (seqpt[seqnm].a_res[itt].nres > 0 ) {
	/*
	fprintf(stderr, "[%d] freeing %d:%d[%d]:%o\n",
		worker,i,seqnm,itt,seqpt[seqnm].a_res[itt].res);
	*/
	seqpt[seqnm].a_res[itt].nres = 0;
	free(seqpt[seqnm].a_res[itt].res);
      }
    }
    else {
      w_abort(" have walign_done but no walign_dflag");
    }
    seqpt[seqnm].walign_dflg[itt] = 0;
  }
}
