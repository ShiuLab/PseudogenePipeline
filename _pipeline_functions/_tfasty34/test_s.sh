#!/bin/csh -f
echo ""
echo "starting fasta34 - protein" `date` "on" `hostname`
echo `uname -a`
echo ""
fasta34 -q -m 6 -Z 100000 mgstm1.aa:1-100 q > test_m1.ok2.html
fasta34 -S -q -z 11 -O test_m1.ok2_p25 -s P250 mgstm1.aa:100-218 q
echo "done"
echo "starting fastxy34" `date`
fastx34 -m 9 -S -q mgtt2_x.seq q > test_t2.xk2
fasty34 -S -q mgtt2_x.seq q > test_t2.yk2
fastx34 -m 9 -S -q -z 2 mgstm1.esq a > test_m1.xk2z2
fasty34 -S -q -z 2 mgstm1.esq a > test_m1.yk2z2
echo "done"
echo "starting fastxy34 rev" `date`
fastx34 -m 9 -q -m 5 mgstm1.rev q > test_m1.xk2r
fasty34 -q -m 5 -M 200-300 -z 2 mgstm1.rev q > test_m1.yk2rz2
fasty34 -q -m 5 -z 11 mgstm1.rev q > test_m1.yk2rz11
echo "done"
echo "starting ssearch34" `date`
ssearch34 -m 9 -S -z 3 -q mgstm1.aa  q > test_m1.ssz3
ssearch34 -q -M 200-300 -z 2 -Z 100000 -s P250 mgstm1.aa q > test_m1.ss_p25
echo "done"
echo "starting fasta34 - DNA" `date`
fasta34 -q -z 2 mgstm1.seq %RMB 4 > test_m1.ok4z2
fasta34 -q mgstm1.rev %RMB 4 > test_m1.ok4r
echo "done"
echo "starting tfasta34" `date`
tfasta34 -q mgstm1.aa %RMB > test_m1.tk2
echo "done"
echo "starting tfastxy34" `date`
tfastx34 -m 9 -q -i -3 -m 6 mgstm1.aa %p > test_m1.tx2.html
tfasty34 -q -i -3 -N 5000 mgstm1.aa %p > test_m1.ty2
echo "done"
echo "starting fastf34" `date`
fastf34 -q m1r.aa q > test_mf.ff_s
echo "done"
echo "starting tfastf34" `date`
tfastf34 -q -E 0.0001 m1r.aa %r > test_mf.tf_r
echo "done"
echo "starting fasts34" `date`
fasts34 -q n0.aa q > test_m1.fs_s
echo "done"
echo "starting tfasts34" `date`
tfasts34 -q n0.aa %r > test_m1.ts_r
echo "done"
echo "done" `date`
