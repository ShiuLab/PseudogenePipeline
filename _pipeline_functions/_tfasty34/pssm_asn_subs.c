/* pssm_asn_subs.c */


/* $Name: fa_34_26_5 $ - $Id: pssm_asn_subs.c,v 1.15 2007/04/02 18:08:11 wrp Exp $ */

/* copyright (C) 2005 by William R. Pearson and the U. of Virginia */

/* this code is designed to parse the ASN.1 binary encoded scoremat
   object produced by blastpgp -C file.ckpt_asn -u 2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"

int parse_pssm_asn();
int parse_pssm2_asn();

int
parse_pssm_asn_fa(FILE *afd, int *n_rows, int *n_cols,
		  unsigned char **query, double ***freqs,
		  char *matrix, int *gap_open, int *gap_extend,
		  double *lambda);



#define COMPO_NUM_TRUE_AA 20

/**positions of true characters in protein alphabet*/
/*
static int trueCharPositions[COMPO_NUM_TRUE_AA] = {
  1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22
};
*/

#define COMPO_LARGEST_ALPHABET 28

/*
static char ncbieaatoa[COMPO_LARGEST_ALPHABET] = {"-ABCDEFGHIJKLMNOPQRSTUVWXYZ"};

static int alphaConvert[COMPO_LARGEST_ALPHABET] = {
  (-1), 0, (-1), 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15,
  16, 19,   17, (-1), 18, (-1), (-1), (-1), (-1), (-1)
};
*/

int pssm_aa_order[20] = { 1,  /*A*/
			  16, /*R*/
			  13, /*N*/
			   4, /*D*/
			   3, /*C*/
			  15, /*Q*/
			   5, /*E*/
			   7, /*G*/
			   8, /*H*/
			   9, /*I*/
			  11, /*L*/
			  10, /*K*/
			  12, /*M*/
			   6, /*F*/
			  14, /*P*/
			  17, /*S*/
			  18, /*T*/
			  20, /*W*/
			  22, /*Y*/
			  19}; /*V*/


#define ASN_SEQ 48
#define ASN_SEQOF 49

#define ASN_PSSM_QUERY 166
#define ASN_PSSM2_QUERY 162

#define ASN_PSSM_IS_PROT 160
#define ASN_PSSM2_MATRIX 161
#define ASN_PSSM_NROWS 162
#define ASN_PSSM_NCOLS 163

#define ASN_PSSM2_NCOLS 163
#define ASN_PSSM2_NROWS 164
#define ASN_PSSM_BYCOL 165
#define ASN_PSSM_INTERMED_DATA 167
#define ASN_PSSM_FREQS 162
#define ASN_PSSM2_FREQS 165
#define ASN_PSSM2_LAMBDA 166

#define ASN_IS_STR 26
#define ASN_IS_INT  2
#define ASN_IS_BOOL 1
#define ASN_IS_ENUM 10

struct asn_bstruct {
  FILE *fd;
  unsigned char *buf;
  unsigned char *abp;
  unsigned char *buf_max;
  int len;
};

#define ASN_BUF 1024

unsigned char *
chk_asn_buf(struct asn_bstruct *asnp, int v) {
  int new_buf;
  
  if (v > ASN_BUF) {
    fprintf(stderr," attempt to read %d bytes ASN.1 data > buffer size (%d)\n",
	    v, ASN_BUF);
    exit(1);
  }

  if (asnp->abp + v > asnp->buf_max) {

    /* move down the left over stuff */
    asnp->len = asnp->buf_max - asnp->abp;

    memmove(asnp->buf, asnp->abp, asnp->len);

    asnp->abp = asnp->buf;
    new_buf = ASN_BUF - asnp->len;
    
    if (!feof(asnp->fd) && 
	(new_buf=fread(asnp->buf + asnp->len, sizeof(char), new_buf, asnp->fd)) != 0) {
      asnp->len += new_buf;
    }

    asnp->buf_max = asnp->buf + asnp->len;

    if (asnp->len < v) {
      fprintf(stderr, " Unable to read %d bytes\n",v);
      exit(1);
    }
  }
  /* otherwise, v bytes are currently in the buffer */

  return asnp->abp;
}

/* read_asn_dest reads v bytes into oct_str if v <= o_len */
/* read_asn_dest is required for ASN data entities that are longer than ASN_BUF (1024) */
unsigned char *
read_asn_dest(struct asn_bstruct *asnp, int v, unsigned char *oct_str, int o_len) {
  int new_buf;
  unsigned char *oct_ptr;
  

  if (v > o_len) {
    fprintf(stderr, " read_asn_dest - cannot read %d bytes into %d buffer\n",
	    v, o_len);
    exit(1);
  }

  if (asnp->abp + v <= asnp->buf_max) {
    memmove(oct_str, asnp->abp, v);
    return asnp->abp+v;
  }
  else {
    /* move down the left over stuff */

    asnp->len = asnp->buf_max - asnp->abp;

    memmove(oct_str, asnp->abp, asnp->len);
    oct_ptr = oct_str+asnp->len;
    v -= asnp->len;

    asnp->abp = asnp->buf;
    new_buf = ASN_BUF;
    
    while ((new_buf=fread(asnp->buf, sizeof(char), new_buf, asnp->fd)) != 0) {
      asnp->len = new_buf;
      asnp->buf_max = asnp->buf + asnp->len;
      if (v <= new_buf) {	/* we have it all this time */
	memmove(oct_ptr, asnp->buf, v);
	asnp->len -= v;
	asnp->abp = asnp->buf + v;
	break;
      }
      else {	/* we need to read some more */
	memmove(oct_ptr, asnp->buf, new_buf);
	v -= new_buf;
	new_buf = ASN_BUF;
      }
    }
  }
  return asnp->buf + v;
}

unsigned char *
get_astr_bool(struct asn_bstruct *asnp, int *val) {

  int v_len, v;

  asnp->abp = chk_asn_buf(asnp,5);

  v = 0;
  if (*asnp->abp++ != 1) { /* check for int */
    fprintf(stderr," bool missing\n");
  }
  else {
    v_len = *asnp->abp++;
    if (v_len != 1) {
      fprintf(stderr, "boolean length != 1 : %d\n", v_len);
      v = *asnp->abp++;
    }
    else { v = *asnp->abp++;}
  }
  asnp->abp += 2;	/* skip over null's */
  *val = v;
  return asnp->abp;
}

unsigned char *
get_astr_int(struct asn_bstruct *asnp,
	    int *val) {

  int v_len, v;

  v = 0;

  asnp->abp = chk_asn_buf(asnp,8);

  if (*asnp->abp++ != 2) { /* check for int */
    fprintf(stderr," int missing\n");
  }
  else {
    v_len = *asnp->abp++;
    while (v_len-- > 0) {
      v *= 256;
      v += *asnp->abp++;
    }
    asnp->abp += 2;	/* skip over null's */
  }
  *val = v;
  return asnp->abp;
}

unsigned char *
get_astr_enum(struct asn_bstruct *asnp, int *val) {

  int v_len, v;

  asnp->abp = chk_asn_buf(asnp,5);

  v = 0;
  if (*asnp->abp++ != ASN_IS_ENUM) { /* check for int */
    fprintf(stderr," enum missing\n");
  }
  else {
    v_len = *asnp->abp++;
    while (v_len-- > 0) { v *= 256;  v += *asnp->abp++; }
    asnp->abp += 2;	/* skip over null's */
  }
  *val = v;

  return asnp->abp;
}

unsigned char *
get_astr_packedfloat(struct asn_bstruct *asnp, double *val) {

  int v_len, v;
  char tmp_str[64];

  asnp->abp = chk_asn_buf(asnp,2);

  v = 0;
  if (*asnp->abp++ != 9) { /* check for packed float */
    fprintf(stderr," float missing\n");
    *val = 0;
    return asnp->abp;
  }
  else {
    v_len = *asnp->abp++;

    if (v_len > 63) {
      fprintf(stderr," real string too long: %d\n",v_len);
    }

    asnp->abp = chk_asn_buf(asnp,v_len);

    if (v_len == 2  && *asnp->abp == '\0' && *(asnp->abp+1)=='0') {
      asnp->abp += 2;
      *val = 0.0;
    }
    else {	/* copy and scan it */
      if (*asnp->abp != '\0') {
	fprintf(stderr, " packedfloat - expected 0, got %d\n", *asnp->abp);
	*val = -1.0;
	return asnp->abp;
      }
      asnp->abp++;
      strncpy(tmp_str, (char *)asnp->abp, sizeof(tmp_str)-1);
      tmp_str[v_len-1] = '\0';
      tmp_str[63] = '\0';
      sscanf(tmp_str,"%lg",val);
      asnp->abp += v_len-1;
    }
  }
  return asnp->abp;
}

unsigned char *
get_astr_str(struct asn_bstruct *asnp, char *text, int t_len) {

  int v_len;

  asnp->abp = chk_asn_buf(asnp,2);

  text[0] = '\0';
  if (*asnp->abp++ != ASN_IS_STR) { /* check for str */
    fprintf(stderr," str missing\n");
  }
  else {
    v_len = *asnp->abp++;
    if (v_len > 128) { /* need to read the length from the next bytes */
      t_len = v_len &0x7f;

      asnp->abp = chk_asn_buf(asnp,t_len);

      for (v_len =0; t_len; t_len--) { v_len = (v_len << 8) + *asnp->abp++; }
    }

    /* read v_len bytes */

    asnp->abp = read_asn_dest(asnp,v_len, (unsigned char *)text, t_len);
    asnp->abp += 2; /* skip over last nulls */
  }
  return asnp->abp;
}

#define ASN_BIOSEQ_SEQ 160
#define ASN_BIOSEQ_ID  160
#define ASN_BIOSEQ_ID_VAL 160

#define ASN_BIOSEQ_ID_LOCAL 161
#define ASN_BIOSEQ_ID_GIBBSQ 162
#define ASN_BIOSEQ_ID_GIBBMT 163
#define ASN_BIOSEQ_ID_GB 164
#define ASN_BIOSEQ_ID_EMBL 165
#define ASN_BIOSEQ_ID_PIR 166
#define ASN_BIOSEQ_ID_SP 167
#define ASN_BIOSEQ_ID_PATENT 168
#define ASN_BIOSEQ_ID_OTHER 169
#define ASN_BIOSEQ_ID_GEN 170
#define ASN_BIOSEQ_ID_GI 171

#define ASN_BIOSEQ_TEXTID_NAME 160
#define ASN_BIOSEQ_TEXTID_ACC 161
#define ASN_BIOSEQ_TEXTID_REL 162
#define ASN_BIOSEQ_TEXTID_VER 163

#define ASN_BIOSEQ_DESCR  161
#define ASN_BIOSEQ_INST  162
#define ASN_BIOSEQ_TITLE  164
#define ASN_BIOSEQ_INST_REPR  160
#define ASN_BIOSEQ_INST_MOL  161
#define ASN_BIOSEQ_INST_LEN  162
#define ASN_BIOSEQ_INST_TOPOL  166
#define ASN_BIOSEQ_INST_SEQD  167
#define ASN_OCTET_STR  65
#define ASN_NCBIeaa 65

unsigned char *
get_astr_seqdescr(struct asn_bstruct *asnp,
		 char *descr) {

  int end_seq=0;
  
  /* get seqof '1' */
  /* get 164/128 -  title */
  /* get string */
  /* pop nulls */

  asnp->abp = chk_asn_buf(asnp,6);

  if (*asnp->abp == ASN_SEQOF) {
    end_seq++;
    asnp->abp += 2;
  }
  else {
    fprintf(stderr, " missing ASN_SEQOF '1': %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }

  if (*asnp->abp == ASN_BIOSEQ_TITLE) { 
    asnp->abp+=2;
    asnp->abp = get_astr_str(asnp, descr, MAX_STR);
  }
  else {
    fprintf(stderr, " missing ASN_BIOSEQ_TITLE '1': %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }

  asnp->abp = chk_asn_buf(asnp,2);

  asnp->abp += 2; /* skip over nulls */

  return asnp->abp;
}

unsigned char *
get_astr_octstr(struct asn_bstruct *asnp,
	       unsigned char *oct_str,
	       int o_len) {

  int q_len, v_len;

  asnp->abp = chk_asn_buf(asnp,2);

  if (*asnp->abp++ == ASN_NCBIeaa) {
    /* get length  of length */
    if (*asnp->abp > 128) {
      v_len = *asnp->abp++ & 0x7f;

      asnp->abp = chk_asn_buf(asnp,v_len);

      q_len = 0;
      while (v_len-- > 0) {
	q_len *= 256;
	q_len += *asnp->abp++;
      }
    }
    else {
      q_len = *asnp->abp++ & 0x7f;
    }

    asnp->abp = read_asn_dest(asnp, q_len, oct_str, o_len);

    oct_str[min(q_len,o_len)]='\0';

    asnp->abp += 2;	/* skip characters and NULL's */
  }
  return asnp->abp;
}

unsigned char *
get_astr_seqinst(struct asn_bstruct *asnp,
		unsigned char **query,
		int *nq) {

  int end_seq=0, tmp;

  /* get sequence '0' */
  /* get 160/128/10/len/val -  repr enum raw val */
  /* get 161/128/10/len/val -  mol enum aa val */
  /* get 162/128/02/len/val -  length int val */
  /* get 166/128 - topology (empty) */
  /* get 167/128 - seq-data */
  /* get 65/len+128/len/octet_string */
  /* pop nulls */

  asnp->abp = chk_asn_buf(asnp,12);

  if (*asnp->abp == ASN_SEQ) {
    end_seq++;
    asnp->abp += 2;
  }
  else {
    fprintf(stderr, " missing ASN_SEQ '0': %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }

  if (*asnp->abp == ASN_BIOSEQ_INST_REPR && *(asnp->abp+1) == 128) {
    asnp->abp+=2;
    asnp->abp = get_astr_enum(asnp, &tmp);
  }
  else {
    fprintf(stderr, " missing ASN_BIOSEQ_INST_REPR 160: %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }

  if (*asnp->abp == ASN_BIOSEQ_INST_MOL && *(asnp->abp+1) == 128) {
    asnp->abp+=2;
    asnp->abp = get_astr_enum(asnp, &tmp);
  }
  else {
    fprintf(stderr, " missing ASN_BIOSEQ_INST_MOL 161: %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }

  if (*asnp->abp == ASN_BIOSEQ_INST_LEN) {
    asnp->abp+=2;
    asnp->abp = get_astr_int(asnp, nq);
  }
  else {
    fprintf(stderr, " missing ASN_BIOSEQ_INST_LEN 161: %0x %0x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }

  if ((*query = (unsigned char *)calloc(*nq + 1, sizeof(char)))==NULL) {
    fprintf(stderr, " cannot read %d char query\n", *nq+1);
  }

  if (*asnp->abp == ASN_BIOSEQ_INST_TOPOL && *(asnp->abp+1) == 128 ) {
    asnp->abp += 2;
  }

  if (*asnp->abp == ASN_BIOSEQ_INST_SEQD) {
    asnp->abp+=2;
    asnp->abp = get_astr_octstr(asnp, *query, *nq );
  }
  else {
    fprintf(stderr, " missing ASN_BIOSEQ_INST_SEQD 166: %0x %0x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }

  asnp->abp += 4; /* skip over nulls */

  return asnp->abp;
}


unsigned char *
get_astr_textid( struct asn_bstruct *asnp,
		char *name,
		char *acc) {
  int end_seq = 0;
  int ver;

  chk_asn_buf(asnp,16);

  if (*asnp->abp != ASN_SEQ) {
    fprintf(stderr, " Expected ASN_SEQ: %0x %0x\n",*asnp->abp, asnp->abp[1]);
  }
  else {asnp->abp += 2; end_seq++;}

  name[0] = acc[0] = '\0';
  
  while (*asnp->abp != '\0') {
    if (*asnp->abp == ASN_BIOSEQ_TEXTID_NAME) {
      asnp->abp+=2;
      asnp->abp = get_astr_str(asnp, name, MAX_SSTR);
    }
    if (*asnp->abp == ASN_BIOSEQ_TEXTID_ACC) {
      asnp->abp+=2;
      asnp->abp = get_astr_str(asnp, acc, MAX_SSTR);
    }
    if (*asnp->abp == ASN_BIOSEQ_TEXTID_VER) {
      asnp->abp+=2;
      asnp->abp = get_astr_int(asnp, &ver);
    }
  }
  asnp->abp += 4;
  while (end_seq-- > 0) { asnp->abp += 4; }
  return asnp->abp;
}

unsigned char *
get_astr_query(struct asn_bstruct *asnp,
	      int *gi,
	      char *name,
	      char *acc,
	      char *descr,
	      unsigned char **query,
	      int *nq
	      ) {

  int end_seq = 0;

  asnp->abp = chk_asn_buf(asnp,32);

  if (*asnp->abp != ASN_BIOSEQ_SEQ) {
    fprintf(stderr, "Bioseq - missing SEQ 1: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { asnp->abp += 2;}

  if (*asnp->abp != ASN_SEQ && *asnp->abp != ASN_SEQOF ) {
    fprintf(stderr, "Bioseq - missing SEQUENCE tag 1: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { 
    end_seq++;
    asnp->abp += 2;
  }

  if (*asnp->abp != ASN_BIOSEQ_ID) {
    fprintf(stderr, "Bioseq - missing ID tag: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else {
    asnp->abp += 2;
    if (*asnp->abp != ASN_SEQOF) {
      fprintf(stderr, "missing bioseq/id/SEQOF tag: %d\n",*asnp->abp);
      return asnp->abp;
    }
    else {
      asnp->abp += 2;
      if (*asnp->abp == ASN_BIOSEQ_ID_VAL && *(asnp->abp+1)==128) { asnp->abp += 2;}

      if (*asnp->abp == ASN_BIOSEQ_ID_GI ) {
	asnp->abp+=2;
	asnp->abp = get_astr_int(asnp, gi);
      }

      if (*asnp->abp == ASN_BIOSEQ_ID_LOCAL) {
	*gi = 0;
	acc[0] = '\0';

	asnp->abp+=2;
	asnp->abp = get_astr_str(asnp, name, MAX_SSTR);
	asnp->abp += 2;
      }
      else if (*asnp->abp == ASN_BIOSEQ_ID_SP || *asnp->abp == ASN_BIOSEQ_ID_EMBL ||
	  *asnp->abp == ASN_BIOSEQ_ID_GB || *asnp->abp == ASN_BIOSEQ_ID_PIR ||
	  *asnp->abp == ASN_BIOSEQ_ID_OTHER ) {

	asnp->abp+=2;
	asnp->abp = get_astr_textid(asnp, name, acc);
      }
    }
  }

  while (*asnp->abp == 0) asnp->abp += 2;

  if (*asnp->abp == ASN_BIOSEQ_DESCR) {
    asnp->abp+=2;
    asnp->abp = get_astr_seqdescr(asnp, descr);
    asnp->abp += 2; 		/* skip nulls */
  }
  else { descr[0] = '\0';}

  if (*asnp->abp != ASN_BIOSEQ_INST) {
    fprintf(stderr, "Bioseq - missing ID tag: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { 
    asnp->abp += 2;
    asnp->abp = get_astr_seqinst(asnp, query, nq);
    asnp->abp += 2; 		/* skip nulls */
  }
  return asnp->abp;
}

unsigned char *
get_astr_query2(struct asn_bstruct *asnp,
	      int *gi,
	      char *name,
	      char *acc,
	      char *descr,
	      unsigned char **query,
	      int *nq
	      ) {

  int end_seq = 0;

  asnp->abp = chk_asn_buf(asnp,32);

  if (*asnp->abp != ASN_BIOSEQ_SEQ) {
    fprintf(stderr, "Bioseq - missing SEQ 1: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { asnp->abp += 2;}

  if (*asnp->abp != ASN_SEQOF ) {
    fprintf(stderr, "Bioseq2 - missing SEQOF tag 1: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { 
    end_seq++;
    asnp->abp += 2;
  }

  if (*asnp->abp != ASN_BIOSEQ_ID) {
    fprintf(stderr, "Bioseq - missing ID tag: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else {
    asnp->abp += 2;
    if (*asnp->abp == ASN_SEQOF) {
      asnp->abp += 2;
    }

    if (*asnp->abp == ASN_BIOSEQ_ID_VAL && *(asnp->abp+1)==128) { asnp->abp += 2;}

    if (*asnp->abp == ASN_BIOSEQ_ID_GI ) {
      asnp->abp+=2;
      asnp->abp = get_astr_int(asnp, gi);
    }

    if (*asnp->abp == ASN_BIOSEQ_ID_LOCAL) {
      *gi = 0;
      acc[0] = '\0';

      asnp->abp+=2;
      asnp->abp = get_astr_str(asnp, name, MAX_SSTR);
      asnp->abp += 2;
      }
    else if (*asnp->abp == ASN_BIOSEQ_ID_SP || *asnp->abp == ASN_BIOSEQ_ID_EMBL ||
	     *asnp->abp == ASN_BIOSEQ_ID_GB || *asnp->abp == ASN_BIOSEQ_ID_PIR ||
	     *asnp->abp == ASN_BIOSEQ_ID_OTHER ) {

      asnp->abp+=2;
      asnp->abp = get_astr_textid(asnp, name, acc);
    }
  }

  while (*asnp->abp == 0) asnp->abp += 2;

  if (*asnp->abp == ASN_BIOSEQ_DESCR) {
    asnp->abp+=2;
    asnp->abp = get_astr_seqdescr(asnp, descr);
    asnp->abp += 2; 		/* skip nulls */
  }
  else { descr[0] = '\0';}

  if (*asnp->abp != ASN_BIOSEQ_INST) {
    fprintf(stderr, "Bioseq - missing ID tag: %2x %2x\n",*asnp->abp, asnp->abp[1]);
    return asnp->abp;
  }
  else { 
    asnp->abp += 2;
    asnp->abp = get_astr_seqinst(asnp, query, nq);
    asnp->abp += 2; 		/* skip nulls */
  }
  return asnp->abp;
}

unsigned char *
get_pssm_freqs(struct asn_bstruct *asnp,
	       double **freqs,
	       int n_rows,  
	       int n_cols,
	       int by_row) {

  int i_rows, i_cols;
  int in_seq = 0;

  double f_val;

  asnp->abp = chk_asn_buf(asnp,4);

  if (*asnp->abp == ASN_SEQ) {
    in_seq = 1;
    asnp->abp += 2;
    in_seq = 1;
  }

  if (!by_row) {
    for (i_cols = 0; i_cols < n_cols; i_cols++) {
      for (i_rows = 0; i_rows < n_rows; i_rows++) {
	asnp->abp = get_astr_packedfloat(asnp, &f_val);
	freqs[i_cols][i_rows] = f_val;
      }
    }
  }
  else {
    for (i_rows = 0; i_rows < n_rows; i_rows++) {
      for (i_cols = 0; i_cols < n_cols; i_cols++) {
	asnp->abp = get_astr_packedfloat(asnp, &f_val);
	freqs[i_cols][i_rows] = f_val;
      }
    }
  }
  if (in_seq) {asnp->abp +=2;}	/* skip nulls */
  asnp->abp += 2;
  return asnp->abp;
}

unsigned char *
get_pssm_intermed(struct asn_bstruct *asnp,
		  double **freqs,
		  int n_rows,
		  int n_cols,
		  int by_row) {

  asnp->abp = chk_asn_buf(asnp,4);

  if (*asnp->abp == ASN_SEQ) {
    asnp->abp += 2;
    if (*asnp->abp == ASN_PSSM_FREQS) {
      asnp->abp+=2;
      asnp->abp = get_pssm_freqs(asnp, freqs, n_rows, n_cols, by_row);
    }
    asnp->abp +=2;	/* skip nulls */
  }
  asnp->abp += 2;
  return asnp->abp;
}


#define ASN_PSSM_PARAMS 161
#define ASN_PSSM_PARAMS_PSEUDOCNT 160
#define ASN_PSSM_PARAMS_RPSPARAMS 161
#define ASN_PSSM_RPSPARAMS_MATRIX 160
#define ASN_PSSM_RPSPARAMS_GAPOPEN 161
#define ASN_PSSM_RPSPARAMS_GAPEXT 162

unsigned char *
get_pssm_rpsparams(struct asn_bstruct *asnp,
	       char *matrix,
	       int *gap_open,
	       int *gap_ext) {

  int end_seq=0;

  asnp->abp = chk_asn_buf(asnp,4);

  if (*asnp->abp == ASN_SEQ) {
    asnp->abp += 2;
    end_seq++;
  }

  if (*asnp->abp == ASN_PSSM_RPSPARAMS_MATRIX) {
    asnp->abp+=2;
    asnp->abp = get_astr_str(asnp, matrix, MAX_SSTR);
  }

  if (*asnp->abp == ASN_PSSM_RPSPARAMS_GAPOPEN) {
    asnp->abp+=2;
    asnp->abp = get_astr_int(asnp, gap_open);
  }

  if (*asnp->abp == ASN_PSSM_RPSPARAMS_GAPEXT) {
    asnp->abp+=2;
    asnp->abp = get_astr_int(asnp, gap_ext);
  }

  if (end_seq) { chk_asn_buf(asnp,end_seq * 2); }
  while (end_seq-- > 0) { asnp->abp += 2; }
  return asnp->abp;
}

unsigned char *
get_pssm_params(struct asn_bstruct *asnp,
	       int *pseudo_cnts,
	       char *matrix,
	       int *gap_open,
	       int *gap_ext) {

  int end_seq=0;

  asnp->abp = chk_asn_buf(asnp,6);

  if (*asnp->abp == ASN_SEQ) {
    asnp->abp += 2;
    end_seq++;
  }

  if (*asnp->abp == ASN_PSSM_PARAMS_PSEUDOCNT) {
    asnp->abp+=2;
    asnp->abp = get_astr_int(asnp, pseudo_cnts);
  }

  if (*asnp->abp == ASN_PSSM_PARAMS_RPSPARAMS) {
    asnp->abp+=2;
    asnp->abp = get_pssm_rpsparams(asnp, matrix, gap_open, gap_ext);
    asnp->abp += 2;
  }
  while (end_seq-- > 0) { asnp->abp+=2; }
  return asnp->abp;
}


unsigned char *
get_pssm2_intermed(struct asn_bstruct *asnp,
		   double ***freqs,
		   int n_rows,
		   int n_cols) {

  int i;
  double **my_freqs;

  if ((my_freqs = (double **) calloc(n_cols, sizeof(double *)))==NULL) {
    fprintf(stderr, " cannot allocate freq cols - %d\n", n_cols);
    exit(1);
  }

  if ((my_freqs[0] = (double *) calloc(n_cols * n_rows, sizeof(double)))==NULL) {
    fprintf(stderr, " cannot allocate freq rows * cols - %d * %d\n", n_rows, n_cols);
    exit(1);
  }

  for (i=1; i < n_cols; i++) {
    my_freqs[i] = my_freqs[i-1] + n_rows;
  }

  *freqs = my_freqs;

  chk_asn_buf(asnp, 8);

  return get_pssm_freqs(asnp, my_freqs, n_rows, n_cols, 0);
}

int
parse_pssm2_asn(struct asn_bstruct *asnp,
		int *gi,
		char *name,
		char *acc,
		char *descr,
		unsigned char **query, 
		int *nq,
		int *n_rows,
		int *n_cols,
		double ***freqs,
		int *pseudo_cnts,
		char *matrix, 
		double *lambda_p) {

  int is_protein;
  int have_rows, have_cols;

  chk_asn_buf(asnp, 32);

  if (memcmp(asnp->abp, "\241\2000\200",4) != 0) {
    fprintf(stderr, "improper PSSM2 start\n");
    return -1;
  }
  else {asnp->abp+=4;}

  if (*asnp->abp == ASN_BIOSEQ_SEQ ) {
    asnp->abp = get_astr_query2(asnp, gi, name, acc, descr, query, nq);
  }

  /* finish up the nulls */
  while (*asnp->abp == '\0') { asnp->abp += 2;}

  if (*asnp->abp == ASN_PSSM2_QUERY && 
	asnp->abp[2] != ASN_SEQ ) {
      fprintf(stderr, "improper PSSM2 start\n");
      return -1;
  }
  else {asnp->abp += 4;}

  while (*asnp->abp != '\0' ) {

    switch (*asnp->abp) {
    case ASN_PSSM_IS_PROT :
      asnp->abp+=2;
      asnp->abp = get_astr_bool(asnp, &is_protein);
      break;

    case ASN_PSSM2_MATRIX :
      asnp->abp+=2;
      asnp->abp = get_astr_str(asnp, matrix, MAX_SSTR);
      break;

    case ASN_PSSM2_NROWS :
      asnp->abp+=2;
      asnp->abp = get_astr_int(asnp, n_rows);
      
      if (*n_rows > 0) { have_rows = 1; }
      else {
	fprintf(stderr, " bad n_row count\n");
	exit(1);
      }
      break;

    case  ASN_PSSM2_NCOLS :
      asnp->abp+=2;
      asnp->abp = get_astr_int(asnp, n_cols);
      if (*n_cols > 0) {
	have_cols = 1;
      }
      else {
	fprintf(stderr, " bad n_row count\n");
	exit(1);
      }
      break;
    
    case ASN_PSSM2_FREQS :
      asnp->abp += 4;
      if (*asnp->abp == '\0') { asnp->abp += 4;}
      break;

    case ASN_PSSM2_LAMBDA :
      asnp->abp += 2;
      asnp->abp = get_astr_packedfloat(asnp,lambda_p);
      asnp->abp +=2;	/* skip over end of ASN_PSSM2_LAMBDA */
      break;

    case ASN_PSSM_INTERMED_DATA :
      asnp->abp += 2;
      asnp->abp = get_pssm2_intermed(asnp, freqs, *n_rows, *n_cols);
      asnp->abp += 4;
      break;

    default: asnp->abp += 2;
    }
  }


  return 1;
}

int
parse_pssm_asn(FILE *afd,
	       int *gi,
	       char *name,
	       char *acc,
	       char *descr,
	       unsigned char **query, 
	       int *nq,
	       int *n_rows,
	       int *n_cols,
	       double ***freqs,
	       int *pseudo_cnts,
	       char *matrix,
	       int *gap_open,
	       int *gap_ext,
	       double *lambda_p) {

  int is_protein, pssm_version;
  int i;
  int have_rows, have_cols, by_col;
  double **my_freqs;

  struct asn_bstruct asn_str;

  if ((asn_str.buf = (unsigned char *)calloc(ASN_BUF, sizeof(char))) == NULL ) {
    fprintf(stderr, " cannot allocate asn_buf (%d)\n",ASN_BUF);
    exit(1);
  }

  asn_str.fd = afd;
  asn_str.len = ASN_BUF;
  asn_str.abp = asn_str.buf_max = asn_str.buf + ASN_BUF;

  chk_asn_buf(&asn_str, 32);

  if (memcmp(asn_str.abp, "0\200\240\200",4) != 0) {
    fprintf(stderr, "improper PSSM header -");
    return -1;
  }
  else {asn_str.abp+=4;}

  if (*asn_str.abp == ASN_IS_INT) {
    asn_str.abp = get_astr_int(&asn_str, &pssm_version);
    if (pssm_version != 2) {
      fprintf(stderr, "PSSM2 version mismatch: %d\n",pssm_version);
      return -1;
    }
    *gap_open = *gap_ext = 0;
    return parse_pssm2_asn(&asn_str, gi, name, acc, descr,
			   query, nq,
			   n_rows, n_cols, freqs,
			   pseudo_cnts, matrix,
			   lambda_p);
  }

  if (*asn_str.abp == ASN_SEQ) { asn_str.abp += 2;  }

  if (*asn_str.abp == ASN_PSSM_IS_PROT ) {
    asn_str.abp+=2;
    asn_str.abp = get_astr_bool(&asn_str, &is_protein);
  }

  if (*asn_str.abp == ASN_PSSM_NROWS ) {
    asn_str.abp+=2;
    asn_str.abp = get_astr_int(&asn_str, n_rows);

    if (*n_rows > 0) { have_rows = 1; }
    else {
      fprintf(stderr, " bad n_row count\n");
      exit(1);
    }
  }

  if (*asn_str.abp == ASN_PSSM_NCOLS ) {
    asn_str.abp+=2;
    asn_str.abp = get_astr_int(&asn_str, n_cols);
    if (*n_cols > 0) {
      have_cols = 1;
    }
    else {
      fprintf(stderr, " bad n_row count\n");
      exit(1);
    }
  }

  if (*asn_str.abp == ASN_PSSM_BYCOL ) {
    asn_str.abp+=2;
    asn_str.abp = get_astr_bool(&asn_str, &by_col);
  }

  /* we have read everything up to the query 

     n_cols gives us the query length, which we can allocate;
  */

  if (*asn_str.abp == ASN_PSSM_QUERY ) {
    asn_str.abp+=2;
    asn_str.abp = get_astr_query(&asn_str, gi, name, acc, descr, query, nq);
    *nq = *n_cols;
  }

  /* finish up the nulls */
  while (*asn_str.abp == '\0') { asn_str.abp += 2;}

  if (*asn_str.abp == ASN_PSSM_INTERMED_DATA) {

    if (!have_rows || !have_cols) {
      fprintf(stderr, " cannot allocate freq - missing rows/cols - %d/%d\n",
	      have_rows, have_cols);
      return -1;
    }

    if ((my_freqs = (double **) calloc(*n_cols, sizeof(double *)))==NULL) {
      fprintf(stderr, " cannot allocate freq cols - %d\n", *n_cols);
      return -1;
    }

    if ((my_freqs[0] = (double *) calloc(*n_cols * *n_rows, sizeof(double)))==NULL) {
      fprintf(stderr, " cannot allocate freq rows * cols - %d * %d\n", *n_rows, *n_cols);
      return -1;
    }
    for (i=1; i < *n_cols; i++) {
      my_freqs[i] = my_freqs[i-1] + *n_rows;
    }

    *freqs = my_freqs;

    asn_str.abp+=2;
    asn_str.abp = get_pssm_intermed(&asn_str, my_freqs, *n_rows, *n_cols, by_col);
    asn_str.abp += 4;
  }

  if (*asn_str.abp == ASN_PSSM_PARAMS ) {
      asn_str.abp+=2;
      asn_str.abp = get_pssm_params(&asn_str, pseudo_cnts, matrix, gap_open, gap_ext);
  }
  else if (*asn_str.abp == 0) {asn_str.abp+=2;}
  return 1;
}

int
parse_pssm_asn_fa( FILE *fd, 
		   int *n_rows_p, int *n_cols_p,
		   unsigned char **query,
		   double ***freq2d,
		   char *matrix,
		   int *gap_open_p,
		   int *gap_extend_p,
		   double *lambda_p) {

  int qi, rj;
  int gi;
  double tmp_freqs[COMPO_LARGEST_ALPHABET];
  char name[MAX_SSTR], acc[MAX_SSTR], descr[MAX_STR];
  int  nq;
  int pseudo_cnts, ret_val;

  /* parse the file */

  ret_val = parse_pssm_asn(fd, &gi, name, acc, descr, query, &nq,
			   n_rows_p, n_cols_p, freq2d,
			   &pseudo_cnts, matrix, gap_open_p, gap_extend_p,
			   lambda_p);

  if (ret_val <=0) return ret_val;

  /* transform the frequencies */

  for (qi = 0; qi < *n_cols_p; qi++) {
    for (rj = 0; rj < *n_rows_p; rj++) { tmp_freqs[rj] = (*freq2d)[qi][rj];}

    for (rj = 0; rj < COMPO_NUM_TRUE_AA; rj++) {
      (*freq2d)[qi][rj] = tmp_freqs[pssm_aa_order[rj]];
    }
  }
  return 1;
}
