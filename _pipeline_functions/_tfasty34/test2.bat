rem ""
rem "starting fasta34_t - protein on win32"
rem ""
fasta34_t -q -m 6 -Z 100000 mgstm1.aa:1-100 q > test_m1.ok2_t.html
fasta34_t -S -q -z 11 -O test_m1.ok2_t_p25 -s P250 mgstm1.aa:100-218 q
rem "done"
rem "starting fastxy34_t"
fastx34_t -m 9c -S -q mgtt2_x.seq q 1 > test_t2.xk1_t
fasty34_t -S -q mgtt2_x.seq q > test_t2.yk2_t
fastx34_t -m 9c -S -q -z 2 mgstm1.esq a > test_m1.xk2_tz2
fasty34_t -S -q -z 2 mgstm1.esq a > test_m1.yk2_tz2
rem "done"
rem "starting fastxy34_t rev"
fastx34_t -m 9c -q -m 5 mgstm1.rev q > test_m1.xk2r_t
fasty34_t -q -m 5 -M 200-300 -z 2 mgstm1.rev q > test_m1.yk2r_tz2
fasty34_t -q -m 5 -z 11 mgstm1.rev q > test_m1.yk2rz11_t
rem "done"
rem "starting ssearch34_t"
ssearch34_t -m 9c -S -z 3 -q mgstm1.aa  q > test_m1.ss_tz3
ssearch34_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_t_p25
rem "starting ssearch34_t"
ssearch34sse2_t -m 9c -S -z 3 -q mgstm1.aa  q > test_m1.ss_tz3sse2
ssearch34sse2_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_t_p25sse2
rem "done"
rem "starting prss34"
prss34_t -q -k 1000 -A mgstm1.aa xurt8c.aa  > test_m1.rss
prfx34_t -q -k 1000 -A mgstm1.esq xurt8c.aa > test_m1.rfx
rem "done"
rem "starting fasta34_t - DNA"
fasta34_t -S -q -z 2 mgstm1.seq %M 4 > test_m1.ok4_tz2
fasta34_t -S -q mgstm1.rev %M 4 > test_m1.ok4r_t
rem "done"
rem "starting tfastxy34_t"
tfastx34_t -m 9c -q -i -3 -m 6 mgstm1.aa %p > test_m1.tx2_t.html
tfasty34_t -q -i -3 -N 5000 mgstm1.aa %p > test_m1.ty2_t
rem "done"
rem "starting fastf34_t"
fastf34_t -q m1r.aa q > test_mf.ff_t
fastf34 -q m1r.aa q > test_mf.ff_s
rem "done"
rem "starting tfastf34_t"
tfastf34_t -q m1r.aa %r > test_mf.tf_tr
rem "done"
rem "starting fasts34_t"
fasts34_t -q -V '*?@' ngts.aa q > test_m1.fs1_t
fasts34_t -q ngt.aa q > test_m1.fs_t
fasts34_t -q -n mgstm1.nts m > test_m1.nfs_t
rem "done"
rem "starting tfasts34_t"
tfasts34_t -q n0.aa %r > test_m1.ts_r
rem "done"
rem "starting fasta34 - protein"
fasta34 -q -z 2 mgstm1.aa q 1 > test_m1.ok1z2
fasta34 -q -s P250 mgstm1.aa q > test_m1.ok2_p25 
rem "done"
rem "starting fastx3"
fastx34 -m 9c -q mgstm1.esq q > test_m1.ok2x 
rem "done"
rem "starting fasty3"
fasty34 -q mgstm1.esq q > test_m1.ok2y 
rem "done"
rem "starting fasta34 - DNA "
fasta34 -m 9c -q mgstm1.seq M 4 > test_m1.ok4 
rem "done"
rem "starting ssearch3"
ssearch34 -S -q -z 2 mgstm1.aa q > test_m1.ss_z2
ssearch34 -q -s P250 mgstm1.aa q > test_m1.ss_p25 
ssearch34sse2 -S -q -z 2 mgstm1.aa q > test_m1.ss_z2_sse2
ssearch34sse2 -q -s P250 mgstm1.aa q > test_m1.ss_p25_sse2 
rem "done"
rem "starting tfastxy3"
tfastx34 -q mgstm1.aa M > test_m1.tx2 
tfasty34 -m 9c -q mgstm1.aa M > test_m1.ty2 
rem "done"
rem "starting fasts34"
fasts34 -q -V '@?*' ngts.aa q > test_m1.fs1
fasts34 -q ngt.aa q > test_m1.fs
rem "done"
