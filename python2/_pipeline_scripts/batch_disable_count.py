import os, sys

if len(sys.argv) < 5:
    print "This program is using script_step6_0321_jpl.py or to generate disable_count" 
    print "Usage: batch_disable_count.py     name_of_the_sw.out_file    start_flage_size    final_flage_size    increase_by"
  
    sys.exit()
    

name  = str(sys.argv[1])
start = int(sys.argv[2])
end   = int(sys.argv[3])
add   = int(sys.argv[4])

c = start
while c <= end:
    print "The flank is %i" %c
    os.system( "python ~/codes/script_step6_0321_jpl.py ~/codes/blosum50.matrix %s %i" % (name,c))
    c += add
    