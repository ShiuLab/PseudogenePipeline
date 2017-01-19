/* Concurrent read version */

/* $Name: fa_34_26_5 $ - $Id: defs.h,v 1.26 2006/06/22 02:35:05 wrp Exp $ */

#ifdef SUNOS
#include <sys/stdtypes.h>
#endif

#ifndef IS_BIG_ENDIAN
#if defined(__BIG_ENDIAN__) || defined(_BIG_ENDIAN)
#define IS_BIG_ENDIAN
#else
#undef IS_BIG_ENDIAN
#endif
#endif

#if !defined(MAX_WORKERS) && !defined(PCOMPLIB)
#define MAX_WORKERS 1
#endif

/* 3-Oct-2003 - we can now have 2 nucleotide query types, DNA
   and RNA.  pst.dnaseq can also be SEQT_RNA.
   ldnaseq can only be DNA */

#define SEQT_DNA 1
#define SEQT_RNA 3	/* DNA and RNA seqtypes must be odd */

#define SEQT_PROT 0
#define SEQT_UNK -1
#define SEQT_OTHER 2

#ifndef DEF_NMLEN
#define DEF_NMLEN 6
#endif

/* unfortunately, there is an important relationship between MAXTRN and
   MAXTST+MAXLIB embedded here.  MAXTRN must be >= (MAXTST+MAXLIB)/3
   or it will be possible for a translated DNA sequence to be longer
   than the translation space available */

#define MAX_STR	512 /* standard label/message buffer */
#define MAX_SSTR 32 /* short string */
#define MAX_FN  120 /* maximum size of a file name */
#define MAX_CH	40 /* maximum number of library choices */
#ifndef SMALLMEM
#define MAX_LF  500 /* maximum numer of library files */
#else
#define MAX_LF  80 /* maximum numer of library files */
#endif

/* padding at the end of sequences for ALTIVEC, other vector
   processors */
#define SEQ_PAD 16

#define MAX_UID 20 /* length of libstr, used for character keys with SQL */

#define AVE_AA_LEN 400
#define AVE_NT_LEN 5000
#define MAX_AA_BUF 5000		/* 5000 later */
#define MAX_NT_BUF 1000		/* 2000 later */

#ifndef SMALLMEM
#define MAXTST	40000		/* longest query */
#define MAXLIB	120000		/* longest library */
#define MAXPLIB	600000		/* longest library with p_comp* */
#define MIN_RES 2000		/* minimum amount allocated for alignment */
#ifndef TFAST
#define MAXTRN  80000		/* buffer for fastx translation */
#else
#define MAXTRN 180000		/* buffer for tfastx translation */
#endif
#define SEQDUP	1200		/* future - overlap */
#ifndef PCOMPLIB
#ifndef MAXBEST
#define MAXBEST	60000	/* max number of best scores */
#endif
#define MAXSTATS 60000
#else
#ifndef MAXBEST
#define MAXBEST	60000	/* max number of best scores */
#endif
#define MAXSTATS 60000
#endif
#define BIGNUM  1000000000
#ifndef MAXINT
#define MAXINT 2147483647
#endif
#define MAXLN	120	/* size of a library name */
#else
#define MAXTST	1500
#define MAXLIB	10000
#define MAXPLIB	100000		/* longest library with p_comp* */
#define MIN_RES 1000
#ifndef TFAST
#define MAXTRN  6000
#else
#define MAXTRN 11500
#endif
#define SEQDUP	300
#define MAXBEST 2000
#define MAXSTATS 20000
#define BIGNUM  32767
#define MAXINT  32767
#define MAXLN	40	/* size of a library name */
#endif
#if !defined(TFAST)
#define MAXTOT (MAXTST+MAXLIB)
#define MAXDIAG	(MAXTST+MAXLIB)
#else
#define MAXTOT (MAXTST+MAXTRN)
#define MAXDIAG	(MAXTST+MAXTRN)
#endif

#define MAXPAM	600	/* maximum allowable size of the pam matrix */
#define PROF_MAX 500
#define ALF_MAX 30

#ifdef SUPERFAMNUM
#define NSFCHAR '!'
#endif

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

#define MX_ATYPE 7	/* markx==0,1,2 7=> no alignment */
#define MX_ASEP  8	/* markx==3  - separate lines */
#define MX_AMAP  16	/* markx==4,5 - graphic map */
#define MX_HTML  32	/* markx==6  - HTML */
#define MX_M9SUMM 64	/* markx==9(c) */
#define MX_M10FORM 128	/* markx==10 */

/* codes for -m 9 */
#define SHOW_CODE_ID	1	/* identity only */
#define SHOW_CODE_ALIGN 2	/* encoded alignment */
