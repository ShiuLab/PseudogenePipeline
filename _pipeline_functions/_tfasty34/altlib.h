
/* $Name: fa_34_26_5 $ - $Id: altlib.h,v 1.9 2006/02/07 17:52:06 wrp Exp $ */

/* #ifdef UNIX */
/* ncbi blast 1.3 format */
/*
#define NCBIBL13 11
extern int ncbl_getliba();
extern void ncbl_ranlib();
void ncbl_closelib();
*/
#define NCBIBL20 12
/* #endif */

#ifdef MYSQL_DB
#define MYSQL_LIB 16
#define LASTLIB MYSQL_LIB+1
#endif

#ifdef PGSQL_DB
#define PGSQL_LIB 17
#define LASTLIB PGSQL_LIB+1
#endif

#if !defined (LASTLIB) && defined(NCBIBL20)
#define LASTLIB NCBIBL20+1
#endif
#if !defined (LASTLIB)
#define LASTLIB 10
#endif

#define FASTA_F 0
#define DEFAULT 0
#define FULLGB 1
#define UNIXPIR 2
#define EMBLSWISS 3
#define INTELLIG 4
#define VMSPIR 5
#define GCGBIN 6
#define LASTTXT 6

int agetlib(); void aranlib();	/* pearson fasta format */
int lgetlib(); void lranlib();	/* full uncompressed GB FULLGB*/
int pgetlib(); void pranlib();	/* PIR UNIX protein UNIXPIR */
int egetlib(); void eranlib();	/* EMBL/SWISS-PROT EMBLSWISS */
int igetlib(); void iranlib();	/* Intelligenetics INTELLIG */
int vgetlib(); void vranlib();	/* PIR VMS format */
int gcg_getlib(); void gcg_ranlib();	/* GCG 2bit format */

#ifdef NCBIBL20
extern int ncbl2_getliba(); /* ncbi blast 2.0 format */
extern void ncbl2_ranlib();
void ncbl2_closelib();
#endif

#ifdef MYSQL_DB
extern int mysql_getlib();
extern void mysql_ranlib();
int mysql_closelib();
#endif

int (*getliba[LASTLIB])()={
  agetlib,lgetlib,pgetlib,egetlib,
  igetlib,vgetlib,gcg_getlib,agetlib,
  agetlib,agetlib
#ifdef UNIX
  ,agetlib
#ifdef NCBIBL13
  ,ncbl_getliba
#else
  ,ncbl2_getliba
#endif
#ifdef NCBIBL20
  ,ncbl2_getliba
#endif
#ifdef MYSQL_DB
  ,agetlib
  ,agetlib
  ,agetlib
  ,mysql_getlib
#endif
#endif
};

void (*ranliba[LASTLIB])()={
  aranlib,lranlib,pranlib,eranlib,
  iranlib,vranlib,gcg_ranlib,aranlib,
  aranlib,aranlib
#ifdef UNIX
  ,aranlib
#ifdef NCBIBL13
  ,ncbl_ranlib
#else
  ,ncbl2_ranlib
#endif
#ifdef NCBIBL20
  ,ncbl2_ranlib
#endif
#ifdef MYSQL_DB
  ,aranlib
  ,aranlib
  ,aranlib
  ,mysql_ranlib
#endif
#endif
};


/* mmap()ed functions */
#ifdef USE_MMAP
int agetlibm(); void aranlibm();
int lgetlibm(); void lranlibm();
void vranlibm();
int gcg_getlibm();

int (*getlibam[])()={
  agetlibm,lgetlibm, NULL, NULL,NULL,agetlibm,gcg_getlibm
};

void (*ranlibam[])()={
  aranlibm,lranlibm,NULL,NULL,NULL,vranlibm,vranlibm
};
#endif
