/*
  copyright (c) 1999 William R. Pearson
*/

/* $Name: fa_34_26_5 $ - $Id: mm_file.h,v 1.26 2006/10/05 18:20:40 wrp Exp $ */

/*
  mm_file.h - defines m_file_str for mmap()ed files 
*/

#include <sys/types.h>

#ifndef USE_FSEEKO
#define FSEEK fseek
#define FTELL ftell
typedef long fseek_t;
#else
#define FSEEK fseeko
#define FTELL ftello
typedef off_t fseek_t;
#endif
#define FSEEK_T_DEF

#ifdef HAS_INTTYPES
#include <inttypes.h>
#else
#ifdef WIN32
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#else
typedef long int64_t;
typedef unsigned long uint64_t;
#endif
#endif
#ifdef BIG_LIB64
typedef int64_t MM_OFF;
#else
typedef long MM_OFF;
#endif

#ifdef MYSQL_DB
#include <mysql.h>
#endif
#ifdef PGSQL_DB
#include <libpq-fe.h>
#endif

struct lmf_str {
  FILE *libf;		/* sequence file being read */
  FILE *hfile;		/* BLAST2.0 description file */
  unsigned int *oid_list;	/* oid list for subsets */
  int oid_seqs;		/* start offset for mask array */
  int pref_db;		/* preferred database */
  unsigned int max_oid;	/* start offset for mask array */

  char lb_name[120];	/* file name */
  int lb_type;		/* library type */
  int *sascii;		/* ascii -> sq mapping */

  /* used by flat files */
  char *lline;		/* last line read */
  unsigned char *cpsave;	/* position in line for lgetlib() */
  fseek_t lpos;			/* position in file */

  /* Genbank Flat files */
  int lfflag;		/* flag for CRLF in EMBL CDROM files */

  /* stuff for GCG format files (5,6) */
  int gcg_binary;	/* flag for binary gcg format */
  long gcg_len;		/* length of GCG sequence */

  int bl_lib_pos;	  /* for ncbl2 */
  int bl_format_ver;	  /* blast formatdb version */
  char opt_text[MAX_FN];	  /* text after filename */

  /* used when memory mapping */
  int mm_flg;		/* mmap worked */
  int mmap_fd;		/* mmap_fd */
  char *mmap_base;	/* base */
  char *mmap_addr;	/* current pos */
  long st_size;		/* file size */

  MM_OFF *d_pos_arr;	/* pointer to desc. offsets */
  MM_OFF *s_pos_arr;	/* pointer to seq. offsets */
  MM_OFF *a_pos_arr;	/* pointer to aux offsets */

  /* currently available only for memory mapped files */
  int max_cnt;		/* # database entries */
  int64_t tot_len;	/* total residue length */
  long max_len;		/* maximum sequence lengh */
  int lib_aa;		/* 0 = DNA, 1 = prot */
  char *tmp_buf;	/* temporary buffer */
  int tmp_buf_max;	/* max size */

  /* used for SQL database queries */
  char *sql_db, *sql_query, *sql_getdesc, *sql_getseq;
  int sql_reopen;
  char **sql_uid_arr;	/* indexed by lpos */
  /* used to get sequence data */
  char *sql_seqp;

#ifdef MYSQL_DB
  /* used to open the database */
  MYSQL *mysql_conn;
  MYSQL_RES *mysql_res;
  MYSQL_ROW mysql_row;
#endif

#ifdef PGSQL_DB
  /* used to open the database */
  PGconn *pgsql_conn;
  PGresult *pgsql_res;
#endif

  int (*getlib)(unsigned char *seq, int maxs,
		char *libstr, int n_libstr,
		fseek_t *libpos,
		int *lcont,
		struct lmf_str *lm_fd,
		long *l_off);

  void (*ranlib)(char *str, int cnt,
		 fseek_t libpos, char *libstr,
		 struct lmf_str *lm_fd);
};

