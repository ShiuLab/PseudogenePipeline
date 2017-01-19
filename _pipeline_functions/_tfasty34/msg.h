/* Concurrent read version */

/* $Name: fa_34_26_5 $ - $Id: msg.h,v 1.9 2006/03/17 18:34:59 wrp Exp $ */

/* Cube definitions */

#ifdef PVM_SRC
#define FIRSTNODE	1
#define FIRSTWORK	1
#else
#define FIRSTNODE	1
#define FIRSTWORK	1
#endif

#define MAXNOD		128
#define ALLTYPES        -1
#ifdef IPSC2
#define HOSTPID		99
#define MANAGEPID 	100
#define WORKPID 	101
#else
#define HOSTPID		0
#define MANAGEPID 	0
#define WORKPID 	0
#endif
#define MANAGER		0
#define ALLNODES        -1
#define ALLPIDS         -1
#define STARTTYPE0 	0
#define STARTTYPE1	1
#define STARTTYPE2	2
#define STARTTYPE3	3
#define STARTTYPE4	4
#define STARTTYPE5	5
#define STARTTYPE6	6
#define PARAMTYPE	7
#define HSEQTYPE	3
#define MSEQTYPE	4
#define ONETYPE		5
#define TWOTYPE		6
#define MSEQTYPE0	7
#define MSEQTYPE1	8
#define MSEQTYPE2	8
#define LISTTYPE	10
#define LISTRTYPE	11
#define CODERTYPE	12
#define ALN1TYPE	21
#define ALN2TYPE	22
#define ALN3TYPE	23
#define FINISHED 	16384	/* this must be larger than BFR */

#define DO_SEARCH_FLG 0
#define DO_OPT_FLG 1
#define DO_ALIGN_FLG 2
#define DO_CALC_FLG 3
