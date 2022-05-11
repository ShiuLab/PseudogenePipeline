#This script is designed to change all versions of a name for a species to
#another in a tab delimited file Ex: In all instances of
#ATG203932.1 --> ATG203932
#Created by David E. Hufnagel on July 5, 2012

import sys

inp = open(sys.argv[1])      #input file with old names
out = open(sys.argv[2], "w") #output file with new names
spe = sys.argv[3]            #2 char species ID Ex: At





#Write the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

#Go through inp and make the necessary changes
for line in inp:
    if not line.startswith("#"):
        lineLst = line.split("\t")
        
        #identify keyword for species of interest
        if spe == "Rr":
            keyword = spe
        else:
            print "\n***ERROR: need new lines here!***\n"

        #go through each part of lineLst and make necessary changes
        newLineLst = []
        for word in lineLst:
            newWord = word #in case newWord is not changed new it remains the old word
            if word.startswith(keyword):
                if spe == "Rr":
                    newWord = "_".join(word.split("_")[:-1])
                else:
                    print "\n***ERROR: need new lines here!***\n"
            newLineLst.append(newWord)

        #output new lines
        ind = 0  #index
        while ind < len(newLineLst):
            if ind < len(newLineLst) - 1:
                out.write("%s\t" % (newLineLst[ind]))
            else:
                out.write("%s\n" % (newLineLst[ind].strip()))
            ind += 1
            
                
            




inp.close()
out.close()
