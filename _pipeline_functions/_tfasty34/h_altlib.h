
/* $Name: fa_34_26_5 $ - $Id: h_altlib.h,v 1.2 1999/12/30 01:26:57 wrp Exp $ */


#define LASTENTRY 10
#define LASTLIB 10
#define BINARYGB 9
#define DEFAULT 0
#define FULLGB 1
#define UNIXPIR 2
#define EMBLSWISS 3
#define INTELLIG 4
#define VMSPIR 5

int agetlib_h();	/* pearson fasta format */
int agetntlib_h();	/* pearson fasta format nucleotides */
int vgetlib_h();	/* PIR VMS format */

int (*h_getliba[LASTLIB])()={
	agetlib_h,agetlib_h,agetlib_h,agetlib_h,
	agetlib_h,vgetlib_h,agetlib_h,agetlib_h,
	agetlib_h,agetlib_h};

int (*h_getntliba[LASTLIB])()={
	agetntlib_h,agetntlib_h,agetntlib_h,agetntlib_h,
	agetntlib_h,agetntlib_h,agetntlib_h,agetntlib_h,
	agetntlib_h,agetntlib_h};

