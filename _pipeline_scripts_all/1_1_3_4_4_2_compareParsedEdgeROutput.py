#This script is designed to compare 6 edgeR output files (comparing 4 groups)
#to classify and organized genes into 4 single groups (one gene is upregulated
#and the others are not), 6 double groups (two genes are upregulated and the
#other two are downregulated), 4 triple groups (three genes are upregulated
#and the other is downregulated) and a miscelaneous group where I haven't yet
#quantified the relationship that's going on.
#Created by David E. Hufnagel on Nov 21, 2012
#WARNING: ANALYSIS IS DESIGNED FOR HAVING 4 INPUT GROUPS (6 comparisons)

import sys

inpAvB = open(sys.argv[1])   #input parsed edgeR output file comparing group A to group B
inpAvC = open(sys.argv[2])   #input parsed edgeR output file comparing group A to group C
inpAvD = open(sys.argv[3])   #input parsed edgeR output file comparing group A to group D
inpBvC = open(sys.argv[4])   #input parsed edgeR output file comparing group B to group C
inpBvD = open(sys.argv[5])   #input parsed edgeR output file comparing group B to group D
inpCvD = open(sys.argv[6])   #input parsed edgeR output file comparing group C to group D
out = open(sys.argv[7], "w") #output file with info in the first line of each category (starting with "##") and all the genes in the category after that line in the full edgeR output format.




#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))




inpAvB.close()
inpAvC.close()
inpAvD.close()
inpBvC.close()
inpBvD.close()
inpCvD.close()
out.close()
