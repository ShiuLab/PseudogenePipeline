#!/bin/csh -f
echo "starting fasta34_t - protein" `date`
foreach z ( 1 2 3 6 11 )
fasta34_t -q  -z $z mgstm1.aa a > test_m1_a.ok2_t_${z}
fasta34_t -q  -z $z oohu.aa a > test_m1_b.ok2_t_${z}
fasta34_t -q -S -z $z prio_atepa.aa a > test_m1_c.ok2S_t_${z}
fasta34_t -q -S -z $z h10_human.aa a > test_m1_d.ok2S_t_${z}
end
echo "done"
echo "starting ssearch34_t" `date`
foreach z ( 1 2 3 6 11 )
ssearch34_t -q  -z $z mgstm1.aa a > test_m1_a.ssS_t_${z}
ssearch34_t -q  -z $z oohu.aa a > test_m1_b.ssS_t_${z}
ssearch34_t -q -sBL62 -S -f -11 -z $z prio_atepa.aa a > test_m1_c.ssSbl62_t_${z}
ssearch34_t -q -sBL62 -S -f -11 -z $z h10_human.aa a > test_m1_d.ssSbl62_t_${z}
end
echo "done"
echo "starting fasta34 - protein" `date`
foreach z ( 1 2 3 6 11 )
fasta34 -q  -z $z mgstm1.aa a > test_m1_a.ok2_${z}
fasta34 -q  -z $z oohu.aa a > test_m1_b.ok2_${z}
fasta34 -q -S -sBL62 -f -11 -z $z prio_atepa.aa a > test_m1_c.ok2Sbl62_${z}
fasta34 -q -S -sBL62 -f -11 -z $z h10_human.aa a > test_m1_d.ok2Sbl62_${z}
end
echo "done"
echo "starting ssearch3" `date`
foreach z ( 1 2 3 6 11 )
ssearch34 -q  -z $z mgstm1.aa a > test_m1_a.ssS_${z}
ssearch34 -q  -z $z oohu.aa a > test_m1_b.ssS_${z}
ssearch34 -q -S -z $z prio_atepa.aa a > test_m1_c.ssS_${z}
ssearch34 -q -S -z $z h10_human.aa a > test_m1_d.ssS_${z}
end
echo "done" `date`
