#!/bin/csh -f
echo ""
echo "starting fasta34_t - protein" `date` "on" `hostname`
echo `uname -a`
echo ""
fasta34_t -q -m 6 -Z 100000 mgstm1.aa:1-100 q > test_m1.ok2_t.html
fasta34_t -S -q -z 11 -O test_m1.ok2_t_p25 -s P250 mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy34_t" `date`
fastx34_t -m 9 -S -q mgtt2_x.seq q > test_t2.xk2_t
fasty34_t -S -q mgtt2_x.seq q > test_t2.yk2_t
fastx34_t -m 9 -S -q -z 2 mgstm1.esq a > test_m1.xk2_tz2
fasty34_t -S -q -z 2 mgstm1.esq a > test_m1.yk2_tz2
echo "done"
echo "starting fastxy34_t rev" `date`
fastx34_t -m 9 -q -m 5 mgstm1.rev q > test_m1.xk2r_t
fasty34_t -q -m 5 -M 200-300 -z 2 mgstm1.rev q > test_m1.yk2r_tz2
fasty34_t -q -m 5 -z 11 mgstm1.rev q > test_m1.yk2rz11_t
echo "done"
echo "starting ssearch34_t" `date`
ssearch34_t -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.ss_tz3
ssearch34_t -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_t_p25
echo "done"
echo "starting fasta34_t - DNA" `date`
fasta34_t -q -z 2 mgstm1.seq %M 4 > test_m1.ok4_tz2
fasta34_t -q mgstm1.rev %M 4 > test_m1.ok4r_t
echo "done"
echo "starting tfasta34_t" `date`
tfasta34_t -q mgstm1.aa %M > test_m1.tk2_t
echo "done"
echo "starting tfastxy34_t" `date`
tfastx34_t -m 9 -q -i -3 -m 6 mgstm1.aa %m > test_m1.tx2_t.html
tfasty34_t -q -3 -N 5000 mgstm1.aa %m > test_m1.ty2_t
echo "done"
echo "starting fastf34_t" `date`
fastf34_t -q m1r.aa q > test_mf.ff_s
echo "done"
echo "starting tfastf34_t" `date`
tfastf34_t -q m1r.aa %m > test_mf.tf_r
echo "done"
echo "starting fasts34_t" `date`
fasts34_t -q n0.aa q > test_m1.fs_s
echo "done"
echo "starting tfasts34_t" `date`
tfasts34_t -q n0.aa %m > test_m1.ts_r
echo "done"
echo "starting fasta34 - protein" `date`
fasta34 -q -z 2 mgstm1.aa q > test_m1.ok2z2
fasta34 -q -s P250 mgstm1.aa q > test_m1.ok2_p25 
echo "done"
echo "starting fastx3" `date`
fastx34 -m 9 -q mgstm1.esq q > test_m1.ok2x 
echo "done"
echo "starting fasty3" `date`
fasty34 -q mgstm1.esq q > test_m1.ok2y 
echo "done"
echo "starting fasta34 - DNA " `date`
fasta34 -m 9 -q mgstm1.seq %m 4 > test_m1.ok4 
echo "done"
echo "starting ssearch3" `date`
ssearch34 -S -q -z 2 mgstm1.aa q > test_m1.ss_z2
ssearch34 -q -s P250 mgstm1.aa q > test_m1.ss_p25 
echo "done"
echo "starting tfasta3" `date`
tfasta34 -q mgstm1.aa %m > test_m1.tk2 
echo "done"
echo "starting tfastxy3" `date`
tfastx34 -q mgstm1.aa %m > test_m1.tx2 
tfasty34 -m 9 -q mgstm1.aa %m > test_m1.ty2 
echo "done" `date`
