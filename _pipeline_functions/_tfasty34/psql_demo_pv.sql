xdb.wrplab seqdb_demo wrplab gstmu;
SELECT acc, protein.seq, 'sp|'||acc||'|'||sp_name||' '||descr
 FROM annot INNER JOIN protein USING(prot_id) WHERE annot.db='sp' LIMIT 50000;
SELECT acc, descr FROM annot WHERE acc='#' AND db='sp';
SELECT acc,protein.seq FROM protein INNER JOIN annot USING(prot_id)
 WHERE annot.acc='#' AND db='sp';

