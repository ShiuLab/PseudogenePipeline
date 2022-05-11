#This script is a wrapper for pseudogene timing and related analysis
#The overall idea:
"""
1_3_C_GetPseudogenizationTime.py- -
                                  - - > 1_3_C_JoinPseuDupDataPerGene.py - - > 1_3_C_condenseTimingData.py
1_3_C_GetDupTimeFromKs.py - - - - -
"""

import sys, os

retDupKaKs = sys.argv[1] #input file with Ka/Ks analysis info for retain duplicates
inOutKaKs = sys.argv[2]  #input Ka/Ks info file for in vs. out gene
psInKaKs = sys.argv[3]   #input Ka/Ks info file for pseudogene vs. in gene
psWGD = sys.argv[4]      #an info file having WGD derived pseudogenes in the first col
orthos = sys.argv[5]     #input orthologs file to get Br-At relationships 
scriptDir = "/" + sys.argv[6].strip("/") #The directory holding all the scripts in the wrapper



#WARNING: THINK ABOUT WHETHER YOU WANT THE _MOD VERSION OF THE FIRST SCRIPT OR NOT!!!
#ALSO CONSIDER WHETHER KA OR KS HAS BEEN SET TO NOT BE GREATER THAN A CERTAIN #
#ALSO MAKE SURE THE PS-IN FILE HAS ONLY PS-IN
#os.system("python %s/1_3_C_GetPseudogenizationTime_mod.py %s %s %s %s pseudoTiming.info .007" % (scriptDir, inOutKaKs, psInKaKs, psWGD, orthos))
os.system("python %s/1_3_C_GetPseudogenizationTime.py %s %s %s %s pseudoTiming.info .007" % (scriptDir, inOutKaKs, psInKaKs, psWGD, orthos))
#os.system("python %s/1_3_C_GetPseudogenizationTimeRr.py %s %s %s %s pseudoTiming.info .007" % (scriptDir, inOutKaKs, psInKaKs, psWGD, orthos))
os.system("python %s/1_3_C_GetDupTimeFromKs.py %s dupTime.2col" % (scriptDir, retDupKaKs))
os.system("python %s/1_3_C_JoinPseuDupDataPerGene.py pseudoTiming.info dupTime.2col pseudoSpeDupTiming.info" % (scriptDir))
os.system("python %s/1_3_C_condenseTimingData.py pseudoSpeDupTiming.info pseudoSpeDupTiming.info.cond" % (scriptDir))
