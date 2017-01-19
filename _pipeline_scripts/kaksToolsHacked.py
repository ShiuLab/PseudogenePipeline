__version__='031103_3'
__author__='Runsun Pan'


__history__='''
 '

031103_3: c_KaKsCalculator()

031103_2:

	rough version of load for all formats tested ok

031103_1:

	c_2naSeqs.load( ... )	 

031101_1:

	c_naSeq:
	   change some attributes to properties;
	   in degenList, if a stop codon is found, it's degenList will be
	   [0, 0, 0]
	   
	   [1, 1, 2,   1, 1, 2,	0, 0, 0,  ...] means
			the 3rd codon is a stop codon

	   NOTE: the property 'codons' can be read and written

	   new :
	   stopCodons (= ['TAG', 'TGA', 'TAA'] (universal) or
					 ['AGG', 'AGA', 'TAG', 'TAA'] (mito ..))
	   stopCodonsAt() : return a list of integers which are the
						codon indices where stop codons occur.

	   __call__(						

031031_2:
	made c_naSeq and c_2naSeqs derived from Object (new class)
	
	c_naSeq: use property and delete _getattr_

031031_1

	c_naSeq: allow load(seq, format='10aas')
			 checkAttrs()
	c_2naSeqs:
		new method:
		  load(self, seqs, format='fasta', divider='\n')
		  allow loading seq pairs in several formats
		  
		tested: plain format
		  
	

031030_1:

	c_naSeq: added new argument 'format' to method load()
		accepted formats are:
		['fasta', 'mega', 'embl', 'gcg', 'genbank', 'ig', 'plain']
		

030613_2:

	bugfix for batchRun

030613_1:

	More stingent guard to catch error (dividedByZero, etc). Now
	as long as one error occurs, the variabled is assigned a string
	value to indicate the error --- no longer assigning a '0' to
	an error variable.

030610_1:

	The function 'batchRun' now return
	[results, validResults, sums, counts, means, sds]

	where
		results: a dict in which ALL output data are saved:
				 {'KaByKs':[...], 'Ka':[...], 'Ks':[...], 'vKa':[...], 'vKs':[...]}

		validResults: a dict in which VALID output data are saved:
				 {'KaByKs':[...], 'Ka':[...], 'Ks':[...], 'vKa':[...], 'vKs':[...]}

		sums, counts, means, sds are in format of :
			{'KaByKs':f, 'Ka':f, 'Ks':f, 'vKa':f, 'vKs':f}
			where f's are a float type number

	The sum, count, means and sds are also saved to the output file.			
	
030609_3:

	fix the 'divided by zero' error in the newly added vK[d]


030609_2:

	fix the 'divided by zero' error in the newly added vK[d]
	
030609_1:

	add vK[d] such that the formula Li85 can give vKa and vKs.

030605_1:

	bugfix to accomodate the situations when a[d] or b[d] <=0
	
030603_3:

	fine tune to make it easier to use
 
030603_2:

	fine tune to make it easier to use
	
030603_1:

	successful batch run for c_2naSeqs
 
030602_2:

	bug fix
 
030602_1:

	Add .argCorrection (see description in transList()) 

030601_3:

	746 function calls in 0.106 CPU seconds

030601_2:

	c_2naSeqs: An attemp to change functions degenLists(), transList()
	and paraDict() into attributes and called with self.__getattr__().
	It's successful but seems to slow down a bit:

	030601_1: 955 function calls in 0.213 CPU seconds
	030601_2: 1421 function calls (1237 primitive calls) in 0.276 CPU seconds

	will revert back to 030601_1 in the next version.

030601_1:

	add attributes 'treatDegen3As', 'geneCodeType' to c_2naSeqs
 
030531_5:

	more detailed testing routines for c_2naSeqs
 
030531_4:

	c_2naSeqs tested; Ka, Ks, vKa, vKs printed ok but don't know if
	the values are correct.
	Will remove old 'non-class' functions in the nexe version.

030531_3:

	c_naSeq tested successful. One bug: the geneCodeDict doesn't contain
	the correct info for the code with degenerate level=3
	
030531_2:

	complete the first draft of untested c_naSeq, c_2naSeqs

030531_1:

	c_naSeq, c_2naSeqs
 
030530_2:

	codonToAasDict_TCAG
	
030530_1:

	geneCodeDict and codonDegenDict built-in to this file
 
030529_2:

	action93: batch operation for Li's 1993 method
	action85: modified formula of calculating Ks by adding B4 to Ks

030529_1:

	action85: batch operation for Li's 1985 method

030528_5:

	first Ka, Ks result
	
030528_4:

	more functions
	
030528_3:

	more functions
	
030528_2:

	more functions

030528_1:

	complete description of Li's method for calc distance

030527_1:

	adding description of Li's method for calc distance

'''
import sys, os, time, math, string
from pprint import pprint
	

#=============================================================
#
# Setup and data match tables
#
#=============================================================

geneCodeDict = {'universal':{'TTT':'Phe','TCT':'Ser','TAT':'Tyr','TGT':'Cys',
							 'TTC':'Phe','TCC':'Ser','TAC':'Tyr','TGC':'Cys',
							 'TTA':'Leu','TCA':'Ser','TAA':'Stop','TGA':'Stop',
							 'TTG':'Leu','TCG':'Ser','TAG':'Stop','TGG':'Trp',
							 
							 'CTT':'Leu','CCT':'Pro','CAT':'His','CGT':'Arg',
							 'CTC':'Leu','CCC':'Pro','CAC':'His','CGC':'Arg',
							 'CTA':'Leu','CCA':'Pro','CAA':'Gln','CGA':'Arg',
							 'CTG':'Leu','CCG':'Pro','CAG':'Gln','CGG':'Arg',

							 'ATT':'Ile','ACT':'Thr','AAT':'Asn','AGT':'Ser',
							 'ATC':'Ile','ACC':'Thr','AAC':'Asn','AGC':'Ser',
							 'ATA':'Ile','ACA':'Thr','AAA':'Lys','AGA':'Arg',
							 'ATG':'Met','ACG':'Thr','AAG':'Lys','AGG':'Arg',

							 'GTT':'Val','GCT':'Ala','GAT':'Asp','GGT':'Gly',
							 'GTC':'Val','GCC':'Ala','GAC':'Asp','GGC':'Gly',
							 'GTA':'Val','GCA':'Ala','GAA':'Glu','GGA':'Gly',
							 'GTG':'Val','GCG':'Ala','GAG':'Glu','GGG':'Gly'
							 },
			 'mitochondria':{'TTT':'Phe','TCT':'Ser','TAT':'Tyr','TGT':'Cys',
							 'TTC':'Phe','TCC':'Ser','TAC':'Tyr','TGC':'Cys',
							 'TTA':'Leu','TCA':'Ser','TAA':'Stop','TGA':'Trp',
							 'TTG':'Leu','TCG':'Ser','TAG':'Stop','TGG':'Trp',
							 
							 'CTT':'Leu','CCT':'Pro','CAT':'His','CGT':'Arg',
							 'CTC':'Leu','CCC':'Pro','CAC':'His','CGC':'Arg',
							 'CTA':'Leu','CCA':'Pro','CAA':'Gln','CGA':'Arg',
							 'CTG':'Leu','CCG':'Pro','CAG':'Gln','CGG':'Arg',

							 'ATT':'Ile','ACT':'Thr','AAT':'Asn','AGT':'Ser',
							 'ATC':'Ile','ACC':'Thr','AAC':'Asn','AGC':'Ser',
							 'ATA':'Met','ACA':'Thr','AAA':'Lys','AGA':'Stop',
							 'ATG':'Met','ACG':'Thr','AAG':'Lys','AGG':'Stop',

							 'GTT':'Val','GCT':'Ala','GAT':'Asp','GGT':'Gly',
							 'GTC':'Val','GCC':'Ala','GAC':'Asp','GGC':'Gly',
							 'GTA':'Val','GCA':'Ala','GAA':'Glu','GGA':'Gly',
							 'GTG':'Val','GCG':'Ala','GAG':'Glu','GGG':'Gly'
							 }
				}

#===============================================================================
# codonToAasDict_TCAG['universal']['TTT']=
#   [
#	['Phe','Leu','Ile','Val'],  #=> the aas by changing the 1st code to TCAG
#	['Phe','Ser','Tyr','Cys'],  #=> the aas by changing the 2nd code to TCAG
#	['Phe','Phe','Leu','Leu']   #=> the aas by changing the 3rd code to TCAG
#   ]
#
# The original purpose of this table is to calculate the degen list:
#
# ['Phe','Phe','Leu','Leu'] means the 3rd position is 2-fold deneracy
#	  So the degen list for 'TTT' is [1,1,2]
#
# 'TCT':[['Ser','Pro','Thr','Ala'],['Phe','Ser','Tyr','Cys'],['Ser','Ser','Ser','Ser']],
# ==> degen list = [1,1,4]
#
# 'TTG':[['Leu','Leu','Met','Val'],['Leu','Ser','Stop','Trp'],['Phe','Phe','Leu','Leu']],
# ==> degen list = [2,1,2]
#
#===============================================================================
codonToAasDict_TCAG={
 'universal':{
	'TTT':[['Phe','Leu','Ile','Val'],['Phe','Ser','Tyr','Cys'],['Phe','Phe','Leu','Leu']],
	'TCT':[['Ser','Pro','Thr','Ala'],['Phe','Ser','Tyr','Cys'],['Ser','Ser','Ser','Ser']],
	'TAT':[['Tyr','His','Asn','Asp'],['Phe','Ser','Tyr','Cys'],['Tyr','Tyr','Stop','Stop']],
	'TGT':[['Cys','Arg','Ser','Gly'],['Phe','Ser','Tyr','Cys'],['Cys','Cys','Stop','Trp']],
	'TTC':[['Phe','Leu','Ile','Val'],['Phe','Ser','Tyr','Cys'],['Phe','Phe','Leu','Leu']],
	'TCC':[['Ser','Pro','Thr','Ala'],['Phe','Ser','Tyr','Cys'],['Ser','Ser','Ser','Ser']],
	'TAC':[['Tyr','His','Asn','Asp'],['Phe','Ser','Tyr','Cys'],['Tyr','Tyr','Stop','Stop']],
	'TGC':[['Cys','Arg','Ser','Gly'],['Phe','Ser','Tyr','Cys'],['Cys','Cys','Stop','Trp']],
	'TTA':[['Leu','Leu','Ile','Val'],['Leu','Ser','Stop','Stop'],['Phe','Phe','Leu','Leu']],
	'TCA':[['Ser','Pro','Thr','Ala'],['Leu','Ser','Stop','Stop'],['Ser','Ser','Ser','Ser']],
	#'TAA':[['Stop','Gln','Lys','Glu'],['Leu','Ser','Stop','Stop'],['Tyr','Tyr','Stop','Stop']],
	#'TGA':[['Stop','Arg','Arg','Gly'],['Leu','Ser','Stop','Stop'],['Cys','Cys','Stop','Trp']],
	'TTG':[['Leu','Leu','Met','Val'],['Leu','Ser','Stop','Trp'],['Phe','Phe','Leu','Leu']],
	'TCG':[['Ser','Pro','Thr','Ala'],['Leu','Ser','Stop','Trp'],['Ser','Ser','Ser','Ser']],
	#'TAG':[['Stop','Gln','Lys','Glu'],['Leu','Ser','Stop','Trp'],['Tyr','Tyr','Stop','Stop']],
	'TGG':[['Trp','Arg','Arg','Gly'],['Leu','Ser','Stop','Trp'],['Cys','Cys','Stop','Trp']],
	'CTT':[['Phe','Leu','Ile','Val'],['Leu','Pro','His','Arg'],['Leu','Leu','Leu','Leu']],
	'CCT':[['Ser','Pro','Thr','Ala'],['Leu','Pro','His','Arg'],['Pro','Pro','Pro','Pro']],
	'CAT':[['Tyr','His','Asn','Asp'],['Leu','Pro','His','Arg'],['His','His','Gln','Gln']],
	'CGT':[['Cys','Arg','Ser','Gly'],['Leu','Pro','His','Arg'],['Arg','Arg','Arg','Arg']],
	'CTC':[['Phe','Leu','Ile','Val'],['Leu','Pro','His','Arg'],['Leu','Leu','Leu','Leu']],
	'CCC':[['Ser','Pro','Thr','Ala'],['Leu','Pro','His','Arg'],['Pro','Pro','Pro','Pro']],
	'CAC':[['Tyr','His','Asn','Asp'],['Leu','Pro','His','Arg'],['His','His','Gln','Gln']],
	'CGC':[['Cys','Arg','Ser','Gly'],['Leu','Pro','His','Arg'],['Arg','Arg','Arg','Arg']],
	'CTA':[['Leu','Leu','Ile','Val'],['Leu','Pro','Gln','Arg'],['Leu','Leu','Leu','Leu']],
	'CCA':[['Ser','Pro','Thr','Ala'],['Leu','Pro','Gln','Arg'],['Pro','Pro','Pro','Pro']],
	'CAA':[['Stop','Gln','Lys','Glu'],['Leu','Pro','Gln','Arg'],['His','His','Gln','Gln']],
	'CGA':[['Stop','Arg','Arg','Gly'],['Leu','Pro','Gln','Arg'],['Arg','Arg','Arg','Arg']],
	'CTG':[['Leu','Leu','Met','Val'],['Leu','Pro','Gln','Arg'],['Leu','Leu','Leu','Leu']],
	'CCG':[['Ser','Pro','Thr','Ala'],['Leu','Pro','Gln','Arg'],['Pro','Pro','Pro','Pro']],
	'CAG':[['Stop','Gln','Lys','Glu'],['Leu','Pro','Gln','Arg'],['His','His','Gln','Gln']],
	'CGG':[['Trp','Arg','Arg','Gly'],['Leu','Pro','Gln','Arg'],['Arg','Arg','Arg','Arg']],
	'ATT':[['Phe','Leu','Ile','Val'],['Ile','Thr','Asn','Ser'],['Ile','Ile','Ile','Met']],
	'ACT':[['Ser','Pro','Thr','Ala'],['Ile','Thr','Asn','Ser'],['Thr','Thr','Thr','Thr']],
	'AAT':[['Tyr','His','Asn','Asp'],['Ile','Thr','Asn','Ser'],['Asn','Asn','Lys','Lys']],
	'AGT':[['Cys','Arg','Ser','Gly'],['Ile','Thr','Asn','Ser'],['Ser','Ser','Arg','Arg']],
	'ATC':[['Phe','Leu','Ile','Val'],['Ile','Thr','Asn','Ser'],['Ile','Ile','Ile','Met']],
	'ACC':[['Ser','Pro','Thr','Ala'],['Ile','Thr','Asn','Ser'],['Thr','Thr','Thr','Thr']],
	'AAC':[['Tyr','His','Asn','Asp'],['Ile','Thr','Asn','Ser'],['Asn','Asn','Lys','Lys']],
	'AGC':[['Cys','Arg','Ser','Gly'],['Ile','Thr','Asn','Ser'],['Ser','Ser','Arg','Arg']],
	'ATA':[['Leu','Leu','Ile','Val'],['Ile','Thr','Lys','Arg'],['Ile','Ile','Ile','Met']],
	'ACA':[['Ser','Pro','Thr','Ala'],['Ile','Thr','Lys','Arg'],['Thr','Thr','Thr','Thr']],
	'AAA':[['Stop','Gln','Lys','Glu'],['Ile','Thr','Lys','Arg'],['Asn','Asn','Lys','Lys']],
	'AGA':[['Stop','Arg','Arg','Gly'],['Ile','Thr','Lys','Arg'],['Ser','Ser','Arg','Arg']],
	'ATG':[['Leu','Leu','Met','Val'],['Met','Thr','Lys','Arg'],['Ile','Ile','Ile','Met']],
	'ACG':[['Ser','Pro','Thr','Ala'],['Met','Thr','Lys','Arg'],['Thr','Thr','Thr','Thr']],
	'AAG':[['Stop','Gln','Lys','Glu'],['Met','Thr','Lys','Arg'],['Asn','Asn','Lys','Lys']],
	'AGG':[['Trp','Arg','Arg','Gly'],['Met','Thr','Lys','Arg'],['Ser','Ser','Arg','Arg']],
	'GTT':[['Phe','Leu','Ile','Val'],['Val','Ala','Asp','Gly'],['Val','Val','Val','Val']],
	'GCT':[['Ser','Pro','Thr','Ala'],['Val','Ala','Asp','Gly'],['Ala','Ala','Ala','Ala']],
	'GAT':[['Tyr','His','Asn','Asp'],['Val','Ala','Asp','Gly'],['Asp','Asp','Glu','Glu']],
	'GGT':[['Cys','Arg','Ser','Gly'],['Val','Ala','Asp','Gly'],['Gly','Gly','Gly','Gly']],
	'GTC':[['Phe','Leu','Ile','Val'],['Val','Ala','Asp','Gly'],['Val','Val','Val','Val']],
	'GCC':[['Ser','Pro','Thr','Ala'],['Val','Ala','Asp','Gly'],['Ala','Ala','Ala','Ala']],
	'GAC':[['Tyr','His','Asn','Asp'],['Val','Ala','Asp','Gly'],['Asp','Asp','Glu','Glu']],
	'GGC':[['Cys','Arg','Ser','Gly'],['Val','Ala','Asp','Gly'],['Gly','Gly','Gly','Gly']],
	'GTA':[['Leu','Leu','Ile','Val'],['Val','Ala','Glu','Gly'],['Val','Val','Val','Val']],
	'GCA':[['Ser','Pro','Thr','Ala'],['Val','Ala','Glu','Gly'],['Ala','Ala','Ala','Ala']],
	'GAA':[['Stop','Gln','Lys','Glu'],['Val','Ala','Glu','Gly'],['Asp','Asp','Glu','Glu']],
	'GGA':[['Stop','Arg','Arg','Gly'],['Val','Ala','Glu','Gly'],['Gly','Gly','Gly','Gly']],
	'GTG':[['Leu','Leu','Met','Val'],['Val','Ala','Glu','Gly'],['Val','Val','Val','Val']],
	'GCG':[['Ser','Pro','Thr','Ala'],['Val','Ala','Glu','Gly'],['Ala','Ala','Ala','Ala']],
	'GAG':[['Stop','Gln','Lys','Glu'],['Val','Ala','Glu','Gly'],['Asp','Asp','Glu','Glu']],
	'GGG':[['Trp','Arg','Arg','Gly'],['Val','Ala','Glu','Gly'],['Gly','Gly','Gly','Gly']]
	},
 'mitochondria':{
	'TTT':[['Phe','Leu','Ile','Val'],['Phe','Ser','Tyr','Cys'],['Phe','Phe','Leu','Leu']],
	'TCT':[['Ser','Pro','Thr','Ala'],['Phe','Ser','Tyr','Cys'],['Ser','Ser','Ser','Ser']],
	'TAT':[['Tyr','His','Asn','Asp'],['Phe','Ser','Tyr','Cys'],['Tyr','Tyr','Stop','Stop']],
	'TGT':[['Cys','Arg','Ser','Gly'],['Phe','Ser','Tyr','Cys'],['Cys','Cys','Trp','Trp']],
	'TTC':[['Phe','Leu','Ile','Val'],['Phe','Ser','Tyr','Cys'],['Phe','Phe','Leu','Leu']],
	'TCC':[['Ser','Pro','Thr','Ala'],['Phe','Ser','Tyr','Cys'],['Ser','Ser','Ser','Ser']],
	'TAC':[['Tyr','His','Asn','Asp'],['Phe','Ser','Tyr','Cys'],['Tyr','Tyr','Stop','Stop']],
	'TGC':[['Cys','Arg','Ser','Gly'],['Phe','Ser','Tyr','Cys'],['Cys','Cys','Trp','Trp']],
	'TTA':[['Leu','Leu','Met','Val'],['Leu','Ser','Stop','Trp'],['Phe','Phe','Leu','Leu']],
	'TCA':[['Ser','Pro','Thr','Ala'],['Leu','Ser','Stop','Trp'],['Ser','Ser','Ser','Ser']],
	#'TAA':[['Stop','Gln','Lys','Glu'],['Leu','Ser','Stop','Trp'],['Tyr','Tyr','Stop','Stop']],
	'TGA':[['Trp','Arg','Stop','Gly'],['Leu','Ser','Stop','Trp'],['Cys','Cys','Trp','Trp']],
	'TTG':[['Leu','Leu','Met','Val'],['Leu','Ser','Stop','Trp'],['Phe','Phe','Leu','Leu']],
	'TCG':[['Ser','Pro','Thr','Ala'],['Leu','Ser','Stop','Trp'],['Ser','Ser','Ser','Ser']],
	#'TAG':[['Stop','Gln','Lys','Glu'],['Leu','Ser','Stop','Trp'],['Tyr','Tyr','Stop','Stop']],
	'TGG':[['Trp','Arg','Stop','Gly'],['Leu','Ser','Stop','Trp'],['Cys','Cys','Trp','Trp']],
	'CTT':[['Phe','Leu','Ile','Val'],['Leu','Pro','His','Arg'],['Leu','Leu','Leu','Leu']],
	'CCT':[['Ser','Pro','Thr','Ala'],['Leu','Pro','His','Arg'],['Pro','Pro','Pro','Pro']],
	'CAT':[['Tyr','His','Asn','Asp'],['Leu','Pro','His','Arg'],['His','His','Gln','Gln']],
	'CGT':[['Cys','Arg','Ser','Gly'],['Leu','Pro','His','Arg'],['Arg','Arg','Arg','Arg']],
	'CTC':[['Phe','Leu','Ile','Val'],['Leu','Pro','His','Arg'],['Leu','Leu','Leu','Leu']],
	'CCC':[['Ser','Pro','Thr','Ala'],['Leu','Pro','His','Arg'],['Pro','Pro','Pro','Pro']],
	'CAC':[['Tyr','His','Asn','Asp'],['Leu','Pro','His','Arg'],['His','His','Gln','Gln']],
	'CGC':[['Cys','Arg','Ser','Gly'],['Leu','Pro','His','Arg'],['Arg','Arg','Arg','Arg']],
	'CTA':[['Leu','Leu','Met','Val'],['Leu','Pro','Gln','Arg'],['Leu','Leu','Leu','Leu']],
	'CCA':[['Ser','Pro','Thr','Ala'],['Leu','Pro','Gln','Arg'],['Pro','Pro','Pro','Pro']],
	'CAA':[['Stop','Gln','Lys','Glu'],['Leu','Pro','Gln','Arg'],['His','His','Gln','Gln']],
	'CGA':[['Trp','Arg','Stop','Gly'],['Leu','Pro','Gln','Arg'],['Arg','Arg','Arg','Arg']],
	'CTG':[['Leu','Leu','Met','Val'],['Leu','Pro','Gln','Arg'],['Leu','Leu','Leu','Leu']],
	'CCG':[['Ser','Pro','Thr','Ala'],['Leu','Pro','Gln','Arg'],['Pro','Pro','Pro','Pro']],
	'CAG':[['Stop','Gln','Lys','Glu'],['Leu','Pro','Gln','Arg'],['His','His','Gln','Gln']],
	'CGG':[['Trp','Arg','Stop','Gly'],['Leu','Pro','Gln','Arg'],['Arg','Arg','Arg','Arg']],
	'ATT':[['Phe','Leu','Ile','Val'],['Ile','Thr','Asn','Ser'],['Ile','Ile','Met','Met']],
	'ACT':[['Ser','Pro','Thr','Ala'],['Ile','Thr','Asn','Ser'],['Thr','Thr','Thr','Thr']],
	'AAT':[['Tyr','His','Asn','Asp'],['Ile','Thr','Asn','Ser'],['Asn','Asn','Lys','Lys']],
	'AGT':[['Cys','Arg','Ser','Gly'],['Ile','Thr','Asn','Ser'],['Ser','Ser','Stop','Stop']],
	'ATC':[['Phe','Leu','Ile','Val'],['Ile','Thr','Asn','Ser'],['Ile','Ile','Met','Met']],
	'ACC':[['Ser','Pro','Thr','Ala'],['Ile','Thr','Asn','Ser'],['Thr','Thr','Thr','Thr']],
	'AAC':[['Tyr','His','Asn','Asp'],['Ile','Thr','Asn','Ser'],['Asn','Asn','Lys','Lys']],
	'AGC':[['Cys','Arg','Ser','Gly'],['Ile','Thr','Asn','Ser'],['Ser','Ser','Stop','Stop']],
	'ATA':[['Leu','Leu','Met','Val'],['Met','Thr','Lys','Stop'],['Ile','Ile','Met','Met']],
	'ACA':[['Ser','Pro','Thr','Ala'],['Met','Thr','Lys','Stop'],['Thr','Thr','Thr','Thr']],
	'AAA':[['Stop','Gln','Lys','Glu'],['Met','Thr','Lys','Stop'],['Asn','Asn','Lys','Lys']],
	#'AGA':[['Trp','Arg','Stop','Gly'],['Met','Thr','Lys','Stop'],['Ser','Ser','Stop','Stop']],
	'ATG':[['Leu','Leu','Met','Val'],['Met','Thr','Lys','Stop'],['Ile','Ile','Met','Met']],
	'ACG':[['Ser','Pro','Thr','Ala'],['Met','Thr','Lys','Stop'],['Thr','Thr','Thr','Thr']],
	'AAG':[['Stop','Gln','Lys','Glu'],['Met','Thr','Lys','Stop'],['Asn','Asn','Lys','Lys']],
	#'AGG':[['Trp','Arg','Stop','Gly'],['Met','Thr','Lys','Stop'],['Ser','Ser','Stop','Stop']],
	'GTT':[['Phe','Leu','Ile','Val'],['Val','Ala','Asp','Gly'],['Val','Val','Val','Val']],
	'GCT':[['Ser','Pro','Thr','Ala'],['Val','Ala','Asp','Gly'],['Ala','Ala','Ala','Ala']],
	'GAT':[['Tyr','His','Asn','Asp'],['Val','Ala','Asp','Gly'],['Asp','Asp','Glu','Glu']],
	'GGT':[['Cys','Arg','Ser','Gly'],['Val','Ala','Asp','Gly'],['Gly','Gly','Gly','Gly']],
	'GTC':[['Phe','Leu','Ile','Val'],['Val','Ala','Asp','Gly'],['Val','Val','Val','Val']],
	'GCC':[['Ser','Pro','Thr','Ala'],['Val','Ala','Asp','Gly'],['Ala','Ala','Ala','Ala']],
	'GAC':[['Tyr','His','Asn','Asp'],['Val','Ala','Asp','Gly'],['Asp','Asp','Glu','Glu']],
	'GGC':[['Cys','Arg','Ser','Gly'],['Val','Ala','Asp','Gly'],['Gly','Gly','Gly','Gly']],
	'GTA':[['Leu','Leu','Met','Val'],['Val','Ala','Glu','Gly'],['Val','Val','Val','Val']],
	'GCA':[['Ser','Pro','Thr','Ala'],['Val','Ala','Glu','Gly'],['Ala','Ala','Ala','Ala']],
	'GAA':[['Stop','Gln','Lys','Glu'],['Val','Ala','Glu','Gly'],['Asp','Asp','Glu','Glu']],
	'GGA':[['Trp','Arg','Stop','Gly'],['Val','Ala','Glu','Gly'],['Gly','Gly','Gly','Gly']],
	'GTG':[['Leu','Leu','Met','Val'],['Val','Ala','Glu','Gly'],['Val','Val','Val','Val']],
	'GCG':[['Ser','Pro','Thr','Ala'],['Val','Ala','Glu','Gly'],['Ala','Ala','Ala','Ala']],
	'GAG':[['Stop','Gln','Lys','Glu'],['Val','Ala','Glu','Gly'],['Asp','Asp','Glu','Glu']],
	'GGG':[['Trp','Arg','Stop','Gly'],['Val','Ala','Glu','Gly'],['Gly','Gly','Gly','Gly']],
	}
	}
codonDegenDict={'universal':
				  { 'TTT':('Phe',[1, 1, 2]), 'TCT':('Ser',[1, 1, 4]),
					'TAT':('Tyr',[1, 1, 2]), 'TGT':('Cys',[1, 1, 2]),
					'TTC':('Phe',[1, 1, 2]), 'TCC':('Ser',[1, 1, 4]),
					'TAC':('Tyr',[1, 1, 2]), 'TGC':('Cys',[1, 1, 2]),
					'TTA':('Leu',[2, 1, 2]), 'TCA':('Ser',[1, 1, 4]),
					#'TAA':('Stop',[1, 2, 2]),'TGA':('Stop',[1, 2, 1]),
					'TTG':('Leu',[2, 1, 2]), 'TCG':('Ser',[1, 1, 4]),
					#'TAG':('Stop',[1, 1, 2]),
											 'TGG':('Trp',[1, 1, 1]),
					'CTT':('Leu',[1, 1, 4]), 'CCT':('Pro',[1, 1, 4]),
					'CAT':('His',[1, 1, 2]), 'CGT':('Arg',[1, 1, 4]),
					'CTC':('Leu',[1, 1, 4]), 'CCC':('Pro',[1, 1, 4]),
					'CAC':('His',[1, 1, 2]), 'CGC':('Arg',[1, 1, 4]),
					'CTA':('Leu',[2, 1, 4]), 'CCA':('Pro',[1, 1, 4]),
					'CAA':('Gln',[1, 1, 2]), 'CGA':('Arg',[2, 1, 4]),
					'CTG':('Leu',[2, 1, 4]), 'CCG':('Pro',[1, 1, 4]),
					'CAG':('Gln',[1, 1, 2]), 'CGG':('Arg',[2, 1, 4]),
					'ATT':('Ile',[1, 1, 3]), 'ACT':('Thr',[1, 1, 4]),
					'AAT':('Asn',[1, 1, 2]), 'AGT':('Ser',[1, 1, 2]),   
					'ATC':('Ile',[1, 1, 3]), 'ACC':('Thr',[1, 1, 4]),
					'AAC':('Asn',[1, 1, 2]), 'AGC':('Ser',[1, 1, 2]),		
					'ATA':('Ile',[1, 1, 3]), 'ACA':('Thr',[1, 1, 4]),
					'AAA':('Lys',[1, 1, 2]), 'AGA':('Arg',[2, 1, 2]),	 
					'ATG':('Met',[1, 1, 1]), 'ACG':('Thr',[1, 1, 4]),
					'AAG':('Lys',[1, 1, 2]), 'AGG':('Arg',[2, 1, 2]),		
					'GTT':('Val',[1, 1, 4]), 'GCT':('Ala',[1, 1, 4]),
					'GAT':('Asp',[1, 1, 2]), 'GGT':('Gly',[1, 1, 4]),	   
					'GTC':('Val',[1, 1, 4]), 'GCC':('Ala',[1, 1, 4]),
					'GAC':('Asp',[1, 1, 2]), 'GGC':('Gly',[1, 1, 4]),	 
					'GTA':('Val',[1, 1, 4]), 'GCA':('Ala',[1, 1, 4]),
					'GAA':('Glu',[1, 1, 2]), 'GGA':('Gly',[1, 1, 4]),	   
					'GTG':('Val',[1, 1, 4]), 'GCG':('Ala',[1, 1, 4]),
					'GAG':('Glu',[1, 1, 2]), 'GGG':('Gly',[1, 1, 4])
					},
				'mitochondria':
				  { 'TTT':('Phe',[1, 1, 2]),'TCT':('Ser',[1, 1, 4]),
					'TAT':('Tyr',[1, 1, 2]),'TGT':('Cys',[1, 1, 2]),
					'TTC':('Phe',[1, 1, 2]),'TCC':('Ser',[1, 1, 4]),
					'TAC':('Tyr',[1, 1, 2]),'TGC':('Cys',[1, 1, 2]),
					'TTA':('Leu',[2, 1, 2]),'TCA':('Ser',[1, 1, 4]),
					#'TAA':('Stop',[1, 1, 2]),
											'TGA':('Trp',[1, 1, 2]),
					'TTG':('Leu',[2, 1, 2]),'TCG':('Ser',[1, 1, 4]),
					#'TAG':('Stop',[1, 1, 2]),
											'TGG':('Trp',[1, 1, 2]),
					'CTT':('Leu',[1, 1, 4]),'CCT':('Pro',[1, 1, 4]),
					'CAT':('His',[1, 1, 2]),'CGT':('Arg',[1, 1, 4]),
					'CTC':('Leu',[1, 1, 4]),'CCC':('Pro',[1, 1, 4]),
					'CAC':('His',[1, 1, 2]),'CGC':('Arg',[1, 1, 4]),
					'CTA':('Leu',[2, 1, 4]),'CCA':('Pro',[1, 1, 4]),
					'CAA':('Gln',[1, 1, 2]),'CGA':('Arg',[1, 1, 4]),
					'CTG':('Leu',[2, 1, 4]),'CCG':('Pro',[1, 1, 4]),
					'CAG':('Gln',[1, 1, 2]),'CGG':('Arg',[1, 1, 4]),
					'ATT':('Ile',[1, 1, 2]),'ACT':('Thr',[1, 1, 4]),
					'AAT':('Asn',[1, 1, 2]),'AGT':('Ser',[1, 1, 2]),
					'ATC':('Ile',[1, 1, 2]),'ACC':('Thr',[1, 1, 4]),
					'AAC':('Asn',[1, 1, 2]),'AGC':('Ser',[1, 1, 2]),
					'ATA':('Met',[1, 1, 2]),'ACA':('Thr',[1, 1, 4]),
					'AAA':('Lys',[1, 1, 2]),#'AGA':('Stop',[1, 1, 2]),
					'ATG':('Met',[1, 1, 2]),'ACG':('Thr',[1, 1, 4]),
					'AAG':('Lys',[1, 1, 2]),#'AGG':('Stop',[1, 1, 2]),
					'GTT':('Val',[1, 1, 4]),'GCT':('Ala',[1, 1, 4]),
					'GAT':('Asp',[1, 1, 2]),'GGT':('Gly',[1, 1, 4]),
					'GTC':('Val',[1, 1, 4]),'GCC':('Ala',[1, 1, 4]),
					'GAC':('Asp',[1, 1, 2]),'GGC':('Gly',[1, 1, 4]),
					'GTA':('Val',[1, 1, 4]),'GCA':('Ala',[1, 1, 4]),
					'GAA':('Glu',[1, 1, 2]),'GGA':('Gly',[1, 1, 4]),
					'GTG':('Val',[1, 1, 4]),'GCG':('Ala',[1, 1, 4]),
					'GAG':('Glu',[1, 1, 2]),'GGG':('Gly',[1, 1, 4]),
				  }  
				}

						  
transCheckDict={   ('A','G'):'s',
				   ('G','A'):'s',
				   ('T','C'):'s',
				   ('C','T'):'s',
				   ('A','A'):'|',
				   ('G','G'):'|',
				   ('T','T'):'|',
				   ('C','C'):'|',
				   ('A','C'):'v',
				   ('A','T'):'v',
				   ('G','C'):'v',
				   ('G','T'):'v',
				   ('C','A'):'v',
				   ('T','A'):'v',
				   ('C','G'):'v',
				   ('T','G'):'v'
					}
formula={
 'Ka':
	{
	'Li85':'( L[2]*B[2] + L[0]*K[0] ) / ( 2*L[2]/3.0 + L[0] )',
	'Li93':'A[1] + ( L[1]*B[1] + L[2]*B[2] ) / ( L[1] + L[2] )'
	},
 'Ks':
	{	
	'Li85': '( L[2]*A[2] + L[4]*K[4] ) / ( L[2]/3.0 + L[4] )',
	'Li93':'( L[2]*A[2] + L[4]*A[4] ) / ( L[2]+ L[4] ) + B[4]',
	'yhAddB4': '( L[2]*A[2] + L[4]*K[4] ) / ( L[2]/3.0 + L[4] ) +B[4]'
	},
 'vKa':
	{
	'Li85':'9* ( (L[2]**2)*vB[2] + (L[0]**2)* vK[0] ) / ( 2*L[2] + 3*L[0] )*2',
	'Li93':'vA[1] + (( L[1]**2)*vB[1] + (L[2]**2)*vB[2] ) / ( L[1]+L[2] )**2 -2*b[1]*Q[1]*( a[1]*P[1] - c[1]*(1-Q[1]) ) / ( L[1]+L[2] )'
	},
 'vKs':
	{
	'Li85':'9* ( (L[2]**2) *vA[2] + (L[4]**2)* vK[4] ) / ( L[2] + 3* L[4] )**2',
	'Li93':'vB[4] + ( (L[2]**2)*vA[2] + (L[4]**2)*vA[4] ) /( L[2]+L[4] )**2-2*b[4]*Q[4]*( a[4]*P[4] - c[4]*(1-Q[4]) ) /( L[2]+L[4] )'
	}
	}
seqs=[
	'ACTCTCACCCTAGTGTATTTGAGAGAGTTC',
	'TCAACTGAGATGTGTTTAATGGGGGGA',
	'TCGACAGGGATATATCTAATGGGTATA',
	'''	   2	  300
	seq.1						   GTC CGA AGC CAT ATA CTA TTA CTA CTC TTA ACC ATA CCT TTA GTC CTA GTA CCA AGC ATA TAC CTA AGC ACC CAC AAA ACT ACT TTA ACA GTA GGC TTA CCA CTC ATA GAC GTT GTC CTA CTA CTG AAC ATT AAC GTC CTC ATA CTA GGC AAC TAT ACA ATA TTC CTG TCA ATA ACA TGC TCT CAC CAT CAT CTA CAA TAC AGC GCT CTC TTT GTC GTA ATA ACC GTC GCT ATT AAA GCT CTT GTA GCT AAA ATC ATC GGA CTA ACA CCA TCC ATC GAA ATT TAC TTC GCC CTC AAC TTC 
	seq.2						   GTC CGA AGC CAT ATA CTA CTA CTA CTC TTA ACC ATA CCC TCA GTC CTA CTA CCA AGC ACA TAT CTA AGC CCA CAC AAA ATA ACT TTA ACA GTT GGC CTT CCA CCT ATC GAT GTT ATC CTA CTC CTA AAT ATA AAT GTC CTC ATA CTG GAC AAC TAC ACA GTA TTC CTG TCA ATC ATA TGC TCT CAC CAT CAC TTA CAA TAT AGC GCT CTC TTT GTC GTA ATC ACC GTT GCT ATT AAA GCT CTT GTC GCT AAA ATC ATT CGA CTA ACA CCA ACC ATC ATA ATC TAC TTC CCC CTC AAT TTC
	''',
	'CGA',
	'AGT'
	]

LiMethod='''
 Lis distance method:

G1, G2 : sequence notation
N[g]   : length of sequence g where g = 1 or 2
L[d,g] : # of n.a. sites with degenerate level = d for sequence g
		 where d in (1,2,3,4) and g in (1,2)
		 Note:
				 N[g] = L[1,g] + L[2,g] + L[3,g] + L[4,g]

L[d]   : Average L[d,g] (over sequence (1,2)) for degenerate level d
		 where d in (1,2,3,4)

		 L[d] = ( L[d,1] + L[d,2] ) / 2

C[c,g] : codon c for sequence g where c = index of codon
W[c,w] : A pathways to convert C[c,1] to C[c,2] (or vice versa)

		 Let C[c,1]= AAT

			C[c,1] --> C[c,2]  W[c,p]			
		 ----------   ------  --------	 
		 (1)  AAT  -->  AAT	no W
		 (2)  AAT  -->  AAG	= W[c,1]
		 (3)  AAT  -->  ACG	:
		 
			===========================
			  W[c,1]		W[c,2]
			-----------	 -----------
			 AAT(Asn)	   AAT(Asn)	   
				| (asy)		 | (asy)
			 ACT(Thr)	   AAG(Lys)
				| (syn)		 | (asy)
			 ACG(Thr)	   ACG(Thr)
			===========================


NW[c] : Number of W[c,w] for codon c
		 = max(w)
		 = 0  for case (1)
		   1  for case (2)
		   2  for case (3)

PW[c,w]: Probability of W[c,w]. So SUM{w:1~NW[c]} of PW[c,w] =1

		 For case (3):

		 Let PW[c,1] = 0.7
			 PW[c,2] = 0.3

AW[c,w]: # of asynonymous substutitions for pathway W[c,w] for codon c
SW[c,w]: # of synonymous substutitions for pathway W[c,w] for codon c

		  For case (3):
			  --------------------
					   w=1 w=2
			  --------------------  
			  AW[c,w]  1   2
			  SW[c,w]  1   0
			  --------------------
			  PW[c,w] 0.7  0.3
			  --------------------
			  
			
AW[c]  : Average # of asynonymous substutitions for all W[c,w] for codon c
SW[c]  : Average # of synonymous substutitions for all W[c,w] for codon c

			AW[c] = SUM ( PW[c,w] * AW[c,w] ) for w = 1~ NW[c]
			SW[c] = SUM ( PW[c,w] * SW[c,w] ) for w = 1~ NW[c]

			AW[c]= 0.7 x 1 + 0.3 x 2 = 1.3
			SW[c]= 0.7 x 1 + 0.3 x 0 = 0.7


S[d]   : # of transitions (A-G, T-C) for degeneracy level d
V[d]   : # of transversions for degeneracy level d

		 NOTE: for all L[2]:
		 
			-- in mito: all S[2] are syn and all V[2] are nonsyn
			-- in universal: true except: the 1st code of Arg codon
										  (CGA, CGG, AGA, and AGG)
			-- all syn.	subn. are included in S[2]  ?????
			-- all nonsyn. subn. are included in V[2]  ????? 
										  
P[d]   : proportions of transitional differences per site of degenerate level d.
			= S[d]/L[d]
			
Q[d]   : proportions of transversional differences per site of degenerate level d.
			= V[d]/L[d]

A[d]   : mean numbers of transitional subns per site of degenerate level d.
B[d]   : mean numbers of transersional subns per site of degenerate level d.

		A[d] = (1/2) ln (a[d]) - (1/4) ln (b[d])
		B[d] = (1/2) ln (b[d])

		vA[d] = [ a[d]^2P[d] + c[d]^2Q[d] - (a[d]P[d]+c[d]Q[d])^2 ] / L[d]
		vB[d] = b[d]^2Q[d] (1-Q[d]) / L[d]
						
		where:
		 
		 a[d] = 1/(1-2P[d] - Q[d])
		 b[d] = 1/(1-2Q[d])
		 c[d] = (a[d]-b[d])/2

K[d]   : total # of subns per site of degenerate level d.		 

		K[d] = A[d] + B[d]

vK[d]  : variance

		vK[d]=(a[d]^2*P[d]+(b[d]+c[d])^2*Q[d]-(a[d]*P[d]+(b[d]+c[d])*Q[d])^2)/L[d]

Ks[d]  : # of subn per synonymous site for degenerate level d

		Ks[d] = L[d]A[d]
		
		Ks[2] = L[2]A[2] 
		Ks[4] = L[4]A[4]

Ks	 : # of subn per synonymous site

		Ks85  = ( L[2]A[2] + L[4]K[4] ) / ( L[2]/3 + L[4] )
		vKs85 = 9 { L[2]^2 vA[2] + L[4]^2 vK[4] } / ( L[2] + 3 L[4] )^2

		Ka85  = ( L[2]B[2] + L[0]K[0] ) / ( 2L[2]/3 + L[0] )
		vKa85 = 9 { L[2]^2 vB[2] + L[0]^2 vK[0] } / ( 2L[2] + 3L[0] )^2

1993 version:

		Ks93  = ( L[2]A[2] + L[4]A[4] ) / ( L[2] + L[4] ) + B[4]
		vKs93 = vB[4] + { L[2]^2 vA[2] + L[4]^2 vA[4] } / ( L[2] + L[4] )^2 
				- 2b[4]Q[4] { a[4]P[4] - c[4](1-Q[4]) } / ( L[2] + L[4] )

		Ka93  = A[0] + ( L[0]B[0] + L[2]B[2] ) / ( L[0] + L[2] )
		vKa93 = vA[0] + { L[0]^2 vB[0] + L[2]^2 vB[2] } / ( L[0] + L[2] )^2 
				- 2b[0]Q[0] { a[0]P[0] - c[0](1-Q[0]) } / ( L[0] + L[2] )

'''

ln = math.log
Li93={
	'Ka':formula['Ka']['Li93'],
	'Ks':formula['Ks']['Li93'],
	'vKa':formula['vKa']['Li93'],
	'vKs':formula['vKs']['Li93']}

Li85={
	'Ka':formula['Ka']['Li85'], 
	'Ks':formula['Ks']['Li85'],
	'vKa':formula['vKa']['Li85'],
	'vKs':formula['vKs']['Li85']}

Li85_KaAddB4={
	'Ka':formula['Ka']['Li85'] + '+B[4]',
	'Ks':formula['Ks']['Li85'],
	'vKa':formula['vKa']['Li85'],
	'vKs':formula['vKs']['Li85']}


''
def getListStartWith(aList, startsWith, isStrip=1):
	''' for a list:   L= ['abcdef', 'kkddff', 'xyz', '0wer'...],

		getListStartWith(L, 'kk') will return:
		   ['kkddff', 'xyz', '0wer'...],

		getListStartWith(L, 'xy') will return:
		   ['xyz', '0wer'...],

		if isStrip: any item '  xyz' will be considered 'xyz'
		else:	   the spaces in '  xyz' count.
		   
	'''
##			print 'in getListStartWith()'
##			print 'len(aList)=', len(aList)
##			print 'startwith=', startsWith
	tmp = aList[:]
	if isStrip: tmp = [x.strip() for x in tmp]
	startLineIndex = 0
	for i in range(len(tmp)):
##				print '[%i]='%i, tmp[i]
		if tmp[i].startswith(startsWith):
##					print 'startsWitch matched'
			startLineIndex = i
##			print 'startLineIndex=', startLineIndex
	return aList[startLineIndex:]
		
def ifor(a,b,c):
	if a: return b
	else: return c

#=============================================================
#
# classes and tests
#
#=============================================================
class c_naSeq(object):

	def __init__(self, seq='', format='plain', geneCodeType='universal'):
##		print 'in init: seq=', seq
		
		if seq:
##			print 'loading seq'
			self.load(seq, format)
		else  : self.__dict__['__codons']=[]
		self.__dict__['__geneCodeType']=''
		self.__dict__['treatDegen3As']=3
		self.__dict__['remark'] = ''
		self.__dict__['format']=format
##		print 'in init: before setgct, codons = ', self.codons
		
		self.setGeneCodeType(geneCodeType) # This will reload seq
##		print 'in init: after setgct, codons = ', self.codons
		
##		print "leaving init: self.__dict__['__codons'] = ", self.__dict__['__codons']
##		print 'leaving init: self.format = ', self.format
##		print 'leaving init: self.codons = ', self.codons

	def __get_codons(self):	   return self.__dict__['__codons']
	def __set_codons(self, val): self.__dict__['__codons']=val
	def __get_geneCodeType(self): return self.__dict__['__geneCodeType']
	def __get_aas(self):

		gcd = geneCodeDict[self.geneCodeType]
		codons = self.codons
		stopCodons = self.stopCodons
		return [gcd[x] for x in codons if x not in stopCodons]
##		return [geneCodeDict[self.__dict__['__geneCodeType']][x] \
##				for x in self.__dict__['__codons']]

	def __get_geneCodeDict(self):
		return geneCodeDict
	def __get_stopCodons(self):
		gcd = geneCodeDict[self.geneCodeType]
		return [k for k in gcd if gcd[k].lower()=='stop']
	
	stopCodons   = property(__get_stopCodons)	
	geneCodeDict = property(__get_geneCodeDict)				
	geneCodeType = property(__get_geneCodeType)
	codons	   = property(__get_codons, __set_codons)
	aas		  = property(__get_aas)

	def stopCodonsAt(self):
		''' return a list of integers being the indices of codons
			where the codon == a stop codon.			
		'''
		stopCodons = self.stopCodons
		codons	 = self.codons
		tmp = []
		for i in range(len(codons)):
			if codons[i] in stopCodons: tmp.append(i)
		return tmp

	def setGeneCodeType(self, geneCodeType='universal'):
		geneCodeType=geneCodeType.lower()
		if geneCodeType not in ['universal', 'mitochondria']:
			raise ValueError, '''The argument geneCodeType of .setGeneCodeType(),
			entered as %s, should have been either '%s' or '%s'.''' %(
				geneCodeType, 'universal', 'mitochondria')
		self.__dict__['__geneCodeType']=geneCodeType
		self.load(self(), format= self.format)
		return self	

	def load(self, seq, format=None):
		'''
		DNA sequence formats:
		http://www.genomatix.de/online_help/help/sequence_formats.html

		============================
		plain
		============================

		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCC


		============================
		mega
		============================

			   2	  300
		seq.1						   ATG ATC ACA TTT GCC AAC GCA CAA CAA 
		TCC GGA CTA ACA TGA CCA CTT ATC GTC ACT GCC GCC CAC TTA GCC AGC TCC 
		GAT AGC CCA CTC GTC CCC TTA TTA ATT GCC GTT ATT ACA ATC 
		seq.2						   ATC ATC ACA TTT ACA ACC GCA CAA CAA 
		TCC GGA CTA ACA TGA CCA CTT ATC GTC ACC ACC GCC CAC CTA GCA AGC TCC 
		GAA AGT CCA CTA GTC CTC TTA TTA ATC GCC GTC ATA ACA ATC 
		  
		============================
		fasta
		============================

		> Randseq1 first randomly generated seq 
		GGTGGTTACTAACCGTAAGAGATGATGTCGCCGTGGTCGCGTGGCGCCGCGGACCCAGAT
		TGTACTTCTCTGAGTCGTTCTAGATCGACCAGTCTTCTAGCTTGCCCG


		============================
		EMBL
		============================
		
		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237


		============================
		GCGformat
		============================

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
		AA03518  Length: 237  Check: 4514  ..

			   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			 181  tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc


		============================
		GenBank format
		============================

		LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
		DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
					rRNA and 5.8S rRNA genes, partial sequence.
		ACCESSION   U03518
		BASE COUNT	   41 a	 77 c	 67 g	 52 t
		ORIGIN	  
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
		//

		============================
		IG format
		============================

		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		
		============================
		10aas format
		============================
	
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agtt		

		or :

			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237

		'''
		

		def load_mega(seq):
			'''
			============================
			mega
			============================
				   2	  300
			seq.1						   ATG ATC ACA TTT GCC AAC GCA CAA CAA 
			TCC GGA CTA ACA TGA CCA CTT ATC GTC ACT GCC GCC CAC TTA GCC AGC TCC 
			GAT AGC CCA CTC GTC CCC TTA 
			'''
##			seq=seq.split('\n')
			seq=getListStartWith(seq, startsWith= 'seq.', isStrip=1)
			seq= ' '.join(seq)
			seq=seq.upper().split()[1:]
			self.__dict__['__codons']=seq
			
			


		def load_fasta(seq):
			'''
			============================
			fasta
			============================

			> Randseq1 first randomly generated seq 
			GGTGGTTACTAACCGTAAGAGATGATGTCGCCGTGGTCGCGTGGCGCCGCGGACCCAGAT
			TGTACTTCTCTGAGTCGTTCTAGATCGACCAGTCTTCTAGCTTGCCCG
			'''
##			seq=seq.strip().split('\n')
##			for i in range(2): print '\n\n' + seq[i]
			seq = [ line for line in seq if line ][1:]
			seq = ''.join(seq).upper() #.replace(' ','')
			load_plain(seq)
			

		def load_embl(seq):
			'''
			============================
			EMBL
			============================
			
			ID   AA03518	standard; DNA; FUN; 237 BP.
			XX
			AC   U03518;
			XX
			DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
			DE   rRNA and 5.8S rRNA genes, partial sequence.
			XX
			SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
				 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
				 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
				 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
				 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237
			'''
##			seq=seq.strip().split('\n')
			seq = getListStartWith(seq, startsWith= 'SQ', isStrip=1)[1:]
			load_10aas(seq)			

		def load_10aas(seq):
			''' NOTE: seq is entered as a list:
			
			Load formats like:			 

			GeneBank, GCG :
			
					1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
				   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
				  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgca
				  
			EMBL:
			
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237				
			
			'''
##			seq=seq.strip().split('\n')
			seq = ' '.join(seq).upper()
			seq = seq.split(' ')
			newSeq=[]
			for x in seq:
				try:
					x = int(x)   
				except:
					newSeq.append(x)   #append only when it is NOT an int
			 
			seq = ''.join(newSeq)
			load_plain(seq)

			
		def load_genbank(seq):
			'''
			============================
			GenBank format
			============================

			LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
			DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
						rRNA and 5.8S rRNA genes, partial sequence.
			ACCESSION   U03518
			BASE COUNT	   41 a	 77 c	 67 g	 52 t
			ORIGIN	  
					1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
				   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
				  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
				  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
		
			'''
##			seq=seq.strip().split('\n')
			seq = getListStartWith(seq, startsWith= 'ORIGIN', isStrip=1)[1:]
			load_10aas(seq)
						

		def load_gcg(seq):
			'''
			============================
			GCGformat
			============================

			ID   AA03518	standard; DNA; FUN; 237 BP.
			XX
			AC   U03518;
			XX
			DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
			DE   rRNA and 5.8S rRNA genes, partial sequence.
			XX
			SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			AA03518  Length: 237  Check: 4514  ..

				   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
				  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
				 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
				 181  tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
			'''
##			seq=seq.strip().split('\n')
			seq= [line for line in seq if line.strip()] # remove Empty lines
			seq = getListStartWith(seq, startsWith='SQ', isStrip=1)[2:]
			load_10aas(seq)

			
		def load_plain(seq):
##			seq=seq.strip().replace('\n','')
			
			seq = ''.join(seq).strip()
##			print 'in load_plain: seq = ', seq
			r = range(len(seq)-2)
			self.__dict__['__codons']=[]
			[self.__dict__['__codons'].append(seq[i:i+3]) for i in r if i%3==0]
##			print 'leaving load_plain: codons= ', self.codons
##		def load_codons(seq):
##			seq.strip().replace('\n','')
##			self.__dict__['__codons']=seq.split()

		def load_ig(seq):
			'''
			============================
			IG format
			============================

			; comment
			; comment
			U03518
			AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
			TGTTGCTTCGGCGGG
			'''
			
			
##			seq=seq.strip().split('\n')
			seq = ''.join([ x for x in seq if x.strip() and x[0]!=';'][1:])
			load_plain(seq)

					
		seq = seq.split('\n')
		if format==None: format=self.format
		format = format.lower()
		if format in ['fasta', 'mega', 'embl', 'gcg', 'genbank', 'ig',
					  'plain', '10aas']:
##			print 'in load: format = ', format
			self.format = format
			eval('load_%s(seq)' %format)
		else:

			raise ValueError, '''The argument "format", entered as '%s',
			in the function "load()" of c_naSeq must be one of :
			['fasta', 'mega', 'embl', 'gcg', 'genbank', 'ig', 'plain',
			'10aas'] ''' %format
	
		return self

	def degenList(self, treatDegen3As=None):
		''' Return a list of degenerate levels :
			[1, 1, 2,   1, 1, 2,	1, 1, 4,  ...]
			If stop codon, [0,0,0] 

			[1, 1, 2,   1, 1, 2,	0, 0, 0,  ...] means
			the 3rd codon is a stop codon
		'''
##		print '--- in degenList of c_naSeq --- Entering'
		if treatDegen3As==None: treatDegen3As = self.treatDegen3As
		codons = self.__dict__['__codons']
		#print 'codons=', codons
##		print 'na count = ', len(codons)*3
		gct	= self.geneCodeType
		tmp	= []
##		[ tmp.extend(codonDegenDict[gct][codon][1]) for codon in codons \
##		  if codonDegenDict[gct].has_key(codon)]
		for i in range(len(codons)):
			cdn = codons[i]
			if codonDegenDict[gct].has_key(cdn):
				tmp.extend(codonDegenDict[gct][cdn][1])
			else:
				tmp.extend([0,0,0])
				#print 'codon [%i]: %s was not found in codonDegenDict' %(i,cdn)
		#print 'treatDegen3As=', treatDegen3As
		tmp = [ (x==3) and treatDegen3As or x for x in tmp ]
##		print 'degenList count = ', len(tmp)
##		print '--- in degenList of c_naSeq --- Leaving'
		return tmp
		'''
		codonDegenDict={'universal':
				  { 'TTT':('Phe',[1, 1, 2]), 'TCT':('Ser',[1, 1, 4]),
					'TAT':('Tyr',[1, 1, 2]), 'TGT':('Cys',[1, 1, 2]),
					...
		'''
	def degenCountList(self, treatDegen3As=None):
		localDegenList = self.degenList(treatDegen3As)
		c= localDegenList.count
		return [c(1), c(2), c(3), c(4)]		

	def __call__(self, endAtStopCodon=0,
					   excludeStopCodons=0
					   ):
		codons = self.codons
##		print 'in __call__: codons = ', codons
		stopCodons = self.stopCodons
		ra = range(len(codons))

		if endAtStopCodon:

			for i in ra:
				if codons[i] in stopCodons:
					break
			outCodons = codons[:i]
		elif excludeStopCodons:
			outCodons = [cdn for cdn in codons if cdn not in stopCodons]
		else:
			outCodons = codons
		return ''.join(outCodons)
	
	def checkAttrs(self, html=0):
		seq=self()
		codons = self.codons
		t=[]
		t.append('-------------------------------------')
		t.append('seq()	  : %s' % seq)
		t.append('len(seq()) : %i' % len(seq))
		t.append('seq.codons : %s' % codons)
		t.append('codon count: %i' % len(codons))
		t.append('stopCodons : %s' % str(self.stopCodons))
		t.append('stopCodonsAt(): %s' % str(self.stopCodonsAt()))
		t.append('seq.aas   : %s' % self.aas)
		t.append(' aas count: %i' % len(self.aas))
		t.append('seq.degenList()	 : %s' % self.degenList())
		t.append('seq.geneCodeType	: %s' % self.geneCodeType)
		t.append('seq.degenCountList(): %s '% self.degenCountList())
		t.append('-------------------------------------\n')
		t = '\n'.join(t)
		if html:
			t = t.replace('\n','<br>').replace(' ', '&nbsp;')
		return t



def c_naSeq_test(fileName=''):
	
	def testOut():
		t.append('-------------------------------------')
		t.append('seq.geneCodeType: %s' % seq.geneCodeType)
		t.append('seq()	 : %s' % seq())
		t.append('seq.codons: %s' % seq.codons)
		t.append('seq.aas   : %s' % seq.aas)
		t.append('seq.degenList()	 : %s' % seq.degenList())
		t.append('seq.degenCountList(): %s '%seq.degenCountList())
		t.append('-------------------------------------\n')
	t=['#'*70]
	t.append('Testing c_naSeq')
	t.append('#'*70)
	t.append('seqs[0]=%s' %seqs[0])
	t.append('seq = c_naSeq(seqs[0])')
	seq = c_naSeq(seqs[0])
#	testOut()
	t.append(seq.checkAttrs())
	
	t.append("seq.setGeneCodeType('mitochondria')")
	seq.setGeneCodeType('mitochondria')
#	testOut()
	t.append( seq.checkAttrs() )
	
	t.append("seq.setGeneCodeType('universal')")
	seq.setGeneCodeType('universal')
	t.append('seq.load(seqs[1])')
	t.append('seqs[1]=%s' %seqs[1])
	seq.load(seqs[1])
	t.append( seq.checkAttrs() )

	t.append('seqs[2]=%s' %seqs[2])	
	t.append('seq.load(seqs[2])')
	seq.load(seqs[2])
	t.append( seq.checkAttrs() )
	
	t.append('seq.load(seqs[2])')
	seq.load(seqs[2])
	t.append( seq.checkAttrs() )

	if fileName:
		f=open(fileName, 'w')
		f.writelines([x+'\n' for x in t])
		f.close()
		
	print '\n'.join(t)

	

def c_naSeq_test_load():
	samples={
	 'plain':'''AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCC''',
	 'mega':'''
seq.1						   ATG ATC ACA TTT GCC AAC GCA CAA CAA 
TCC GGA CTA ACA TGA CCA CTT ATC GTC ACT GCC GCC CAC TTA GCC AGC TCC 
GAT AGC CCA CTC GTC CCC TTA TTA ATT GCC GTT ATT ACA ATC''',
	
	 'fasta':'''
> Randseq1 first randomly generated seq 
GGTGGTTACTAACCGTAAGAGATGATGTCGCCGTGGTCGCGTGGCGCCGCGGACCCAGAT
TGTACTTCTCTGAGTCGTTCTAGATCGACCAGTCTTCTAGCTTGCCCG''',
	
	 'embl':'''
ID   AA03518	standard; DNA; FUN; 237 BP.
XX
AC   U03518;
XX
DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
DE   rRNA and 5.8S rRNA genes, partial sequence.
XX
SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
	 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
	 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
	 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
	 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237
		''',

	 'gcg':'''
ID   AA03518	standard; DNA; FUN; 237 BP.
XX
AC   U03518;
XX
DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
DE   rRNA and 5.8S rRNA genes, partial sequence.
XX
SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
AA03518  Length: 237  Check: 4514  ..

	   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
	  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
	 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
	 181  tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
		''',
	
	 'genbank':'''
LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
			rRNA and 5.8S rRNA genes, partial sequence.
ACCESSION   U03518
BASE COUNT	   41 a	 77 c	 67 g	 52 t
ORIGIN	  
		1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
	   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
	  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
	  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
		''',

	  'ig':'''
; comment
; comment
U03518
AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		''',
	  '10aas1':'''
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agtt
	  ''',
	 '10aas2':'''
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagtt
	 '''
	}

	def testLoadingFormat(format):
		print '='*80
		print 'Loading [%s] format:' %format
		if format.startswith('10aas'):
			seq.load(samples[format], format='10aas')
		else:
			seq.load(samples[format], format=format)
		sample = samples[format].split('\n')
		for line in sample: print line
		print seq.checkAttrs()

	seq = c_naSeq()

	testLoadingFormat('plain')
	testLoadingFormat('fasta')
	testLoadingFormat('genbank')
	testLoadingFormat('gcg')
	testLoadingFormat('embl')
	testLoadingFormat('mega')
	testLoadingFormat('ig')
	testLoadingFormat('10aas1')
	testLoadingFormat('10aas2')
	
class c_2naSeqs(object):
	'''====================================================================
	 %s, a python class to deal with 2 nucleic acid sequences.

	Attributes:
		.seq1 : class of c_naSeq
		.seq2 : class of c_naSeq
		.geneCodeType : either 'universal' of 'mitochondria'
		.treatDegen3As: an integer between 1~4 standing for the degen level a
						n.a. site of degenerate level 3 should be treated with.
	
	''' %('c_2naSeqs')
	def __init__(self,
				 seq1='',
				 seq2='',
				 geneCodeType  = 'universal',
				 treatDegen3As = 2,
				 argCorrection = 1):
		self.__dict__['seq1']= c_naSeq(seq1,geneCodeType=geneCodeType)
		self.__dict__['seq2']= c_naSeq(seq2,geneCodeType=geneCodeType)
		self.__dict__['__KaKsFormula']={\
		  'Ka' :'(L[2]*B[2]+ L[0]*K[0]) / (2*L[2]/3+ L[0])',
		  'Ks' :'(L[2]*A[2]+ L[4]*K[4]) / (  L[2]/3+ L[4])',
		  'vKa':'9*(vA[2]*L[2]^2+ vK[4]*L[4]^2) / (  L[2]+ 3*L[4])^2',
		   'vKs':'9*(vB[2]*L[2]^2+ vK[0]*L[0]^2) / (2*L[2]+ 3*L[0])^2'}
		
		self.setGeneCodeType ( geneCodeType)
		self.setTreatDegen3As(treatDegen3As)
		self.setArgCorrection(argCorrection)	 
			# argCorrection: see description in transList()'

	def __get_geneCodeDict(self):
		return geneCodeDict

	def __get_stopCodons(self):
		gcd = geneCodeDict[geneCodeType]
		return [k for k in gcd if gcd[k].lower()=='stop']

	geneCodeDict = property(__get_geneCodeDict)	
	stopCodons   = property(__get_stopCodons)	

	def setGeneCodeType(self, geneCodeType):
		if geneCodeType in ['universal','mitochondria']:
			self.__dict__['geneCodeType']= geneCodeType
		return self

	def setTreatDegen3As(self, treatDegen3As):
		if treatDegen3As in [1,2,3,4]:
			self.__dict__['treatDegen3As']=treatDegen3As	
		return self

	def setArgCorrection(self, argCorrection):
		if argCorrection in [0,1]:
			self.argCorrection=argCorrection		
		return self
	
	def degenLists(self):
		#-------------------------------------------------
		# [[1,1,2,   1,1,2,  1,1,4, ...],
		#  [1,1,3,   2,1,2,  2,1,4, ...]]
		#-------------------------------------------------
		return [ self.seq1.degenList(self.treatDegen3As),
				 self.seq2.degenList(self.treatDegen3As) ]

	

	def checkStopCodons(self, endAtStopCodon=0):
		'''
		A seq containing stop codons will cause error when calculating
		Ka, Ks. (The calculation of Ka, Ks beyond a stop codon doesn't
		make sense anyway). This 
		So, if 
		'''
		
	
	def transList(self):
		#-------------------------------------------------
		# ['|', '|', 's', '|', '|', 'v', '|', 's' ...]
		#-------------------------------------------------
		if self.argCorrection:
			'''Correct the transformation of Arg (CGA, CGG, AGA, AGG)
			due to its first code change.
			When the first code in CGA, CGG change to A, or
				 the first code in AGA, AGG change to C,
			they are both 'transitions', but has to be considered
			'transversions'.
			This correction can be turned on/off by changing '.argCorrection'
			030602_1
			'''
			args = ['CGA', 'AGA', 'CGG', 'AGG']			
			t = []
			s1= self.seq1()
			s2= self.seq2()
			for i in range(len(s1)):
				if (i%3==0):  # if it is the first code in the codon
				   #print '[%i] s1=' %i, s1[i:i+3], ' ,s2=', s2[i:i+3]
				   
				   cdn1=s1[i:i+3]
				   cdn2=s2[i:i+3]
				   if (cdn1 in args) and (cdn2 in args) and (cdn1[0]!=cdn2[0]):
					   t.append('s')
				   elif ((cdn1 in ['CGA','CGG']) and s2[i]=='A') \
					or ((cdn1 in ['AGA','AGG']) and s2[i]=='C') \
					or ((cdn2 in ['CGA','CGG']) and s1[i]=='A') \
					or ((cdn2 in ['AGA','AGG']) and s1[i]=='C') :
						t.append('x')
						#print 't=', t
				   else:
						t.append(transCheckDict[(s1[i],s2[i])])
				else:
					t.append(transCheckDict[(s1[i],s2[i])])
			#print 'transList=', len(t), ',' , t
			return t		
		else:
			t= [transCheckDict[x] for x in zip(self.seq1(), self.seq2())]
			#print 'transList=\n', len(t), t
			return t
		
	def paraDict(self):	 # [S[d], V[d], P[d], Q[d]]
		''' Return a dictionary containing values of:
		{'S':S, 'V':V, 'L':L, 'P':P, 'Q':Q,
				'a':a, 'b':b, 'c':c, 'A':A, 'B':B,
				'K':K, 'vA':vA, 'vB':vB }
		'''
##		print 'entering paraDict'
		#-------------------------------------------------
		# [('|',1,1),
		#  ('|',1,1),
		#  ('s',2,3),
		#  ('|',1,2),
		#  ('|',1,1),
		#  ('v',2,2),...]
		#-------------------------------------------------
		zipList = zip(self.transList(),
					  self.degenLists()[0],
					  self.degenLists()[1])

				
		S=[0,0,0,0]
		V=[0,0,0,0]
##		print 'in paraDict 1'
		
		for i in range(len(zipList)):
			y = zipList[i]
			if y[0] =='s':
				S[y[1]-1] += 0.5
				S[y[2]-1] += 0.5
			elif y[0] =='v':	
				V[y[1]-1] += 0.5
				V[y[2]-1] += 0.5
			
			elif y[0] =='x':
				S[2-1] += 0.5
				if y[1]==2:
					V[y[2]-1]+=0.5
				elif y[2] ==2:
					V[y[1]-1]+=0.5		   

##		print 'in paraDict 2'
		L=[ (x+y)/2.0 for x,y in \
			zip(self.seq1.degenCountList(self.treatDegen3As),
				self.seq2.degenCountList(self.treatDegen3As)) ]

##		print 'in paraDict 3'
##		print 'L = ', [str(x) for x in L]
##		P = [ifor((y>0.0) ,(1.0 * x/y),'L=0') for x, y in zip(S,L)]
		P=[]
		for x, y in zip(S,L):
			if y>0.0: P.append(1.0*x/y)
			else	: P.append('L=0')
##		print 'P= ', P
##		
##		print 'in paraDict 4'
##		Q = [ifor((y>0.0),(1.0 * x/y),'L=0') for x, y in zip(V,L)]
		Q=[]
		for x, y in zip(V,L):
			if y>0.0: Q.append(1.0*x/y)
			else	: Q.append('L=0')
##		print 'Q= ', Q
		'''
		print 'S', S
		print 'V', V
		print 'L', L
		print 'P0', [(x,y) for x, y in zip(S,L)]
		print 'Q0', [(x,y) for x, y in zip(V,L)]
		print 'P', [(y>0) and (x/y) or 0 for x, y in zip(S,L)]
		print 'Q', [(y>0) and (x/y) or 0 for x, y in zip(V,L)]
		print 'a', [(1-2*pp-qq) for pp,qq in zip(P,Q)]
		'''

		a=[]
		for pp,qq in zip(P,Q):
			if pp=='L=0' or qq=='L=0':
				a.append('PQ_divBy0(L=0)')
			else:
				if (1-2*pp-qq)==0:
					a.append('a_divBy0(P=%1.3f,Q=%1.3f)'%(pp,qq))
				else:
					a.append(1.0/(1-2*pp-qq))
					
		#a = [ (1-2*pp-qq)==0 and 0 or 1/(1-2*pp-qq) for pp,qq in zip(P,Q)]

		b = []
		''''''
		for qq in Q:
			if qq=='L=0':
				b.append('L=0')
			elif qq==0.5:
				b.append('b_divBy0(Q=0.5)')
			else:
				b.append(1.0/(1-2*qq))
				
		#c = [(aa-bb)/2.0 for aa,bb in zip(a,b)]

		c=[]
		for aa,bb in zip(a,b):

			 # Checking if string first (means an Error value is carrived over
			 # from previous steps :
			 
			 tmp =''
			 if type(aa)==str: tmp =aa
			 if type(bb)==str:
				  if tmp: tmp+='/' + bb
				  else: tmp =bb
				  
			 if tmp:
				 c.append(tmp)
			 else:
				 c.append((aa-bb)/2.0)
				 

		#A = [ ln(aa)/2.0 - ln(bb)/4.0 for aa,bb in zip(a,b) ]
		#A = [ ((aa<=0) or (bb<=0)) and 'n/a' or (ln(aa)/2.0 - ln(bb)/4.0) for aa,bb in zip(a,b) ]
		A=[]
		for aa, bb in zip(a,b):

			 # Checking if string first (means an Error value is carrived over
			 # from previous steps :
			 
			 tmp =''
			 if type(aa)==str: tmp =aa
			 if type(bb)==str:
				  if tmp: tmp+='/' + bb
				  else: tmp =bb
				  
			 if tmp:
				 A.append(tmp)

			 else:
				 tmp=''
				 if aa<=0: tmp = 'ln(a=%1.4f)'%aa
				 if bb<=0:
					 if tmp: tmp+= '/' + 'ln(b=%1.4f)'%bb
					 else:   tmp = 'ln(b=%1.4f)'%bb
				 if tmp:
					 A.append(tmp)
				 else:
					 A.append(ln(aa)/2.0 - ln(bb)/4.0)
				 
			
		#B = [ (bb<=0) and 'n/a' or ln(bb)/2.0 for bb in b ]
		B=[]
		for bb in b:
			
			 if type(bb)==str:
				 B.append( bb )
			 elif bb<=0:
				 B.append( 'ln(b=%1.4f)'%bb )
			 else:
				 B.append(ln(bb)/2.0)		   


		#K = [ ((AA=='n/a') or (BB=='n/a')) and 'n/a' or AA + BB for AA, BB in zip(A,B) ]
		K=[]
		for AA, BB in zip(A,B):
			 # Checking if string first (means an Error value is carrived over
			 # from previous steps :
			 
			 tmp =''
			 if type(AA)==str: tmp =AA
			 if type(BB)==str:
				  if tmp: tmp+='/' + BB
				  else: tmp =BB
				  
			 if tmp:
				 K.append(tmp)
			 else:
				 K.append( AA + BB )

 

		#vA= [ (ll>0) \
		#	  and ((aa**2)*pp + (cc**2)*qq - ( aa*pp+cc*qq)**2)/ll \
		#	  or 0 \
		#	  for aa,cc,pp,qq,ll in zip(a,b,P,Q,L)]
		vA=[]
		for aa,cc,pp,qq,ll in zip(a,c,P,Q,L):
			if type(ll)==str:
				vA.append(ll)
			elif type(aa)==str:
				vA.append(aa)
			elif type(cc)==str:
				vA.append(cc)
			elif type(pp)==str:
				vA.append(pp)
			elif type(qq)==str:
				vA.append(qq)
			else:
				vA.append( ( (aa**2)*pp + (cc**2)*qq - ( aa*pp+cc*qq)**2)/ll )
				 
		#vB= [ (ll>0) \
		#	  and (bb**2)*qq*(1-qq)/ll \
		#	  or 0 \
		#	  for bb, qq, ll in zip(b,Q,L)]
		vB=[]
		for bb,qq,ll in zip(b,Q,L):
			if type(ll)==str:
				vB.append(ll)
			elif type(bb)==str:
				vB.append(bb)
			elif type(qq)==str:
				vB.append(qq)
			else:
				vB.append( (bb**2)*qq*(1-qq)/ll )
		
	   
	  
		vK=[]
		#r = range(4)
		#for i in r:
		#	if L[i]==0.0:
		#		vK.append(0)				
		#	else:
		#		vK.append( (a[i]**2*P[i]+(b[i]+c[i])**2*Q[i]- \
		#				  (a[i]*P[i]+(b[i]+c[i])*Q[i])**2)/L[i])

		for aa,bb,cc,pp,qq,ll in zip(a,b,c,P,Q,L):
			if type(ll)==str:
				vK.append(ll)
			elif type(aa)==str:
				vK.append(aa)
			elif type(bb)==str:
				vK.append(bb)
			elif type(cc)==str:
				vK.append(cc)
			elif type(pp)==str:
				vK.append(pp)
			elif type(qq)==str:
				vK.append(qq)
			else:
				vK.append((aa**2*pp+(bb+cc)**2*qq- \
						  (aa*pp+(bb+cc)*qq)**2)/ll )
		'''

		
		vK=[ ((LL==0.0) or (LL==0)) \
			 and 0 \
			 or (aa**2*pp+(bb+cc)**2*qq-(aa*pp+(bb+cc)*qq)**2)/LL \
			 for aa,bb,cc,pp,qq,LL in zip(a,b,c,P,Q,L)]
	 
		print 'xxx' 
		''' 
		
		S = [S[0]]+S
		V = [V[0]]+V
		L = [L[0]]+L
		P = [P[0]]+P
		Q = [Q[0]]+Q
		a = [a[0]]+a
		b = [b[0]]+b
		c = [c[0]]+c
		A = [A[0]]+A
		B = [B[0]]+B
		vA= [vA[0]]+vA
		vB= [vB[0]]+vB
		K = [K[0]]+K
		vK= [vK[0]]+vK

			 
		
		return {'S':S, 'V':V, 'L':L, 'P':P, 'Q':Q,
				'a':a, 'b':b, 'c':c, 'A':A, 'B':B,
				'K':K, 'vA':vA, 'vB':vB, 'vK':vK }

	def checkVars(self, br='\n'):
		t=[]
		#print 'para1'
		para=self.paraDict()
		#print 'para2'
		
		for k in ['S','V','L','P','Q','a','b','c',
				  'A','B','K','vA','vB']:
			#tmp = ',\t'.join([ (type(x)==float) and ('%1.3f' %x) or str(x) \
			#				   for x in para[k]]).replace("'",'')
			tmp = ',\t'.join([ ifor( (type(x)==float), ('%1.3f' %x), str(x) ) \
							   for x in para[k]]).replace("'",'')
			
			t.append('  %2s = [%s] '%(k, tmp))		
		
		return br.join(t)

	def displayAlign(self):
		t=['Degen1 : ',
		   '  seq1 : ',
		   '  SorV : ',
		   '  seq2 : ',
		   'Degen2 : ']
		n = len(self.seq1())
##		print '--- in displayAlign() ---'
##		print 'len(self.seq1()) = ', len(self.seq1())
##		print 'len(self.seq2()) = ', len(self.seq2())
##		print 'len(self.degenLists()[0])=', len(self.degenLists()[0])
##		print 'len(self.degenLists()[1])=', len(self.degenLists()[1])
##		
##		print 'seq1:\n', self.seq1()
##		print 'seq2:\n', self.seq2()
##		print 'degenLists:\n', self.degenLists()
##		print '--- in displayAlign() ---2'
		for i in range(n):
			t[0]+= str(self.degenLists()[0][i])
			t[1]+= self.seq1()[i]
			t[2]+= self.transList()[i]
			t[3]+= self.seq2()[i]
			t[4]+= str(self.degenLists()[1][i])						
			if i%3==2:
				t=[x+' ' for x in t]
		return t		
					
	def getKaKs(self,
				formula={'Ka':'',
						 'Ks':'',
						 'vKa':'',
						 'vKs':''},
				treatDegen3As='',
				geneCodeType ='',
				argCorrection=None):

		if treatDegen3As in [1,2,3,4]:
			self.treatDegen3As=treatDegen3As
		if geneCodeType.lower() in ['universal','mitochondria']:
			self.geneCodeType =geneCodeType
		if argCorrection in [0,1]:
##			print 'argCorrection=', argCorrection
			self.argCorrection=argCorrection

		para = self.paraDict()

		S = para['S']
		V = para['V']
		L = para['L']		
		P = para['P']
		Q = para['Q']
		a = para['a']
		b = para['b']
		c = para['c']
		A = para['A']
		B = para['B']
		vA= para['vA']
		vB= para['vB']
		K = para['K']
		vK= para['vK']
		
		errorText={
		 'formulaType':
			'''FormulaTypeError:
			The argument "formula" in .getKaKs() of class %s should have been a
			dictionary containing the following 4 keys: 'Ka', 'Ks', 'vKa', 'vKs'
			''' % self.__class__.__name__ ,

		 'formula':
			'''FormulaError:
			The formula for calculating Ka, Ks, vKa & vKs can have the following
			variables, each contains 4 sets of values: S = [s1,s2,s3,s4],
			corresponding to the S's at different degenerate levels. So to get
			the S at degenerate level of 2, you use: S[2].
			
			S : Counts of transitions between seq1 and seq2
			V : Counts of transverions between seq1 and seq2
			L : Average of L1, L2 where L1, L2 are the count of degen. levels.
				Example: L1=[2,3,0,5], L2= [2,2,2,4], L= [1,1.5,1,4.5]
			P : P[d]= S[d]/L[d]
			Q : Q[d]= V[d]/L[d]

			a : a[d] = 1/(1-2P[d] - Q[d])
			b : b[d] = 1/(1-2Q[d])
			c : c[d] = (a[d]-b[d])/2
			
			A : Numbers of transitional subns per site of degenerate level d.
				 A[d] = (1/2) ln (a[d]) - (1/4) ln (b[d])

			B : Numbers of transersional subns per site of degenerate level d.
				 B[d] = (1/2) ln (b[d])

			vA: Variance of A
				vA[d] = [ a[d]^2P[d] + c[d]^2Q[d]- (a[d]P[d]+c[d]Q[d])^2 ] /L[d]
			vB: Variance of B
				vB[d] = b[d]^2Q[d] (1-Q[d]) / L[d]

			K : total # of subns per site of degenerate level d.		 
				K[d] = A[d] + B[d]


			And some operators:

			x/y  : x divided by y 
			x**y : x times to the power of yth
			ln(x): natural log of x
			
			Note: when doing x/y, at least one of x, y should be a float #:
					3/2.0 ==> 1.5
					3.0/2 ==> 1.5
					3/2   ==> 1
			'''
			}
	   
		if type(formula)!= dict:
			raise TypeError, errorText['formulaType']

		if (not formula.has_key('Ka'))  or (not formula.has_key('Ks')) or \
		   (not formula.has_key('vKa')) or (not formula.has_key('vKs')):
			raise TypeError, errorText['formulaType']

		result={'KaByKs':None,'Ka':None, 'Ks':None, 'vKa':None, 'vKs':None}

		#
		# if formula not given, then Li85 is used:
		#
		if formula['Ka'].strip()=='':
		   formula ={'Ka':'(L[2]*B[2]+ L[0]*K[0]) / (2*L[2]/3+ L[0])',
				 'Ks':'(L[2]*A[2]+ L[4]*K[4]) / (  L[2]/3+ L[4])',
				 'vKa':'9*(vA[2]*L[2]^2+ vK[4]*L[4]^2) / (  L[2]+ 3*L[4])^2',
				 'vKs':'9*(vB[2]*L[2]^2+ vK[0]*L[0]^2) / (2*L[2]+ 3*L[0])^2'}

		for k in ['Ka','Ks','vKa','vKs']:
			f = formula[k].replace('^', '**')
##			print 'f = ', f
			if (f!=None) and (f.strip()!=''):
				f=f.replace('[0]','[1]').replace('\n','')
				try:					
					result[k] = eval(f)
				except:
					deno = '/'.join(f.split('/')[1:])
##					print 'deno = ', deno
					result[k]= 'FormulaError'
##			print 'f2= ', f
		try:
			result['KaByKs']= result['Ka']/result['Ks']
		except:
			result['KaByKs']= 'Error'
		return result

	def help(self, keyWord=''):
		
		helpText={}

	def _get_stopCodons(self):pass
		
	def load(self, seqs, format='fasta', divider='\n'):
		'''
		The divider should be in a line all by itself and should be in
		the beginning of that line. So the actual divider taht is used
		to divide the seqs are  \n + divider
		
		============================
		plain
		============================

		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCAACCTGCGGAAGGA
		TCATTACCGAGTGCGGGTCCTTTGGGCC

		AACCTGCGGAATCATTACCGAGTGCGGGTCCTTTGGGCCGGATCATTACCGAGTGC
		GGGTCCTTTGGGCCAACCTGCGGAAGGA
		

		============================
		mega
		============================

			   2	  300
		seq.1						   ATG ATC ACA TTT GCC AAC GCA CAA CAA 
		TCC GGA CTA ACA TGA CCA CTT ATC GTC ACT GCC GCC CAC TTA GCC AGC TCC 
		GAT AGC CCA CTC GTC CCC TTA TTA ATT GCC GTT ATT ACA ATC 
		seq.2						   ATC ATC ACA TTT ACA ACC GCA CAA CAA 
		TCC GGA CTA ACA TGA CCA CTT ATC GTC ACC ACC GCC CAC CTA GCA AGC TCC 
		GAA AGT CCA CTA GTC CTC TTA TTA ATC GCC GTC ATA ACA ATC 
		  
		============================
		fasta
		============================

		> Randseq1 first randomly generated seq 
		GGTGGTTACTAACCGTAAGAGATGATGTCGCCGTGGTCGCGTGGCGCCGCGGACCCAGAT
		TGTACTTCTCTGAGTCGTTCTAGATCGACCAGTCTTCTAGCTTGCCCG
		> Randseq1 first randomly generated seq 
		GGTGGTTACTAACCGTAAGAGATGATGTCGCCTCTTCTAGCTTGCCCGGTGGTCGCGTGG
		CGCCGCGGACCCAGATTGTACTTCTCTGAGTCGTTCTAGATCGACCAG

		============================
		EMBL
		============================
		
		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc		  237

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agt
			 
		============================
		GCGformat
		============================

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
		AA03518  Length: 237  Check: 4514  ..

			   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			 181  tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
		AA03518  Length: 237  Check: 4514  ..

			   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			 181  tgagttgatt gaatg
			 
		============================
		GenBank format
		============================

		LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
		DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
					rRNA and 5.8S rRNA genes, partial sequence.
		ACCESSION   U03518
		BASE COUNT	   41 a	 77 c	 67 g	 52 t
		ORIGIN	  
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc

		LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
		DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
					rRNA and 5.8S rRNA genes, partial sequence.
		ACCESSION   U03518
		BASE COUNT	   41 a	 77 c	 67 g	 52 t
		ORIGIN	  
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
			  
		============================
		IG format
		============================
		
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		'''
		divider = '\n' + divider
		format = format.lower()
		
		def load_plain(seqs):
			'''
			============================
			plain
			============================

			AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCAACCTGCGGAAGGA
			TCATTACCGAGTGCGGGTCCTTTGGGCC

			AACCTGCGGAATCATTACCGAGTGCGGGTCCTTTGGGCCGGATCATTACCGAGTGC
			GGGTCCTTTGGGCCAACCTGCGGAAGGA
			'''
			seqs = seqs.split(divider)

			if len(seqs) !=2:
				raiseSeqParsingError(seqsCount=len(seqs), divider=divider)
				return 
				
			s = c_naSeq(geneCodeType=self.geneCodeType)
			s.load(seqs[0], format='plain')
##			print 's1 codons:\n', s.codons
			self.__dict__['seq1']= s
			
			s = c_naSeq(geneCodeType=self.geneCodeType)
			s.load(seqs[1], format='plain')
##			print 's2 codons:\n', s.codons
			self.__dict__['seq2']= s


		def load_mega(seqs):
			'''
			============================
			mega
			============================

				   2	  300
			seq.1						   ATG ATC ACA TTT GCC AAC GCA CAA CAA 
			TCC GGA CTA ACA TGA CCA CTT ATC GTC ACT GCC GCC CAC TTA GCC AGC TCC 
			GAT AGC CCA CTC GTC CCC TTA TTA ATT GCC GTT ATT ACA ATC 
			seq.2						   ATC ATC ACA TTT ACA ACC GCA CAA CAA 
			TCC GGA CTA ACA TGA CCA CTT ATC GTC ACC ACC GCC CAC CTA GCA AGC TCC 
			GAA AGT CCA CTA GTC CTC TTA TTA ATC GCC GTC ATA ACA ATC
			'''
			seqs = seqs.split('seq.')
			seqs = [x for x in seqs if x.strip()]
			if len(seqs)>2: seqs = seqs[1:] # remove the header line
			
			s1 = ''.join(seqs[0].split(' ')[1:])
			s2 = ''.join(seqs[1].split(' ')[1:])

			s = c_naSeq(geneCodeType=self.geneCodeType)
			s.load(s1, format='plain')
			self.__dict__['seq1']= s
			
			s = c_naSeq(geneCodeType=self.geneCodeType)
			s.load(s2, format='plain')
			self.__dict__['seq2']= s

		def load_fasta(seqs):			
			'''
			============================
			fasta
			============================

			> Randseq1 first randomly generated seq 
			GGTGGTTACTAACCGTAAGAGATGATGTCGCCGTGGTCGCGTGGCGCCGCGGACCCAGAT
			TGTACTTCTCTGAGTCGTTCTAGATCGACCAGTCTTCTAGCTTGCCCG
			> Randseq2 2nd randomly generated seq 
			GGTGGTTACTAACCGTAAGAGATGATGTCGCCTCTTCTAGCTTGCCCGGTGGTCGCGTGG
			CGCCGCGGACCCAGATTGTACTTCTCTGAGTCGTTCTAGATCGACCAG
			'''
			seqs = seqs.split('\n')
			seqs = [x for x in seqs if x.strip()]
			for i in range(len(seqs)):
				seqs[i]= seqs[i].strip()
				if seqs[i].startswith('>'): seqs[i] = ''
			seqs = '\n'.join(seqs)
##			print 'seqs for fasta = ', seqs
			load_plain(seqs)


		def load_embl(seqs):
			'''
		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agt		  237

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
			 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc		60
			 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg	   120
			 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc	   180
			 tgagttgatt gaatgcaatc agt'''
			seqs = seqs.strip().split('\n\n')
			self.seq1.load(seqs[0], format= 'embl')			
			self.seq2.load(seqs[1], format= 'embl')



			

		def load_gcg(seqs):
			'''
		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
		AA03518  Length: 237  Check: 4514  ..

			   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			 181  tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc

		ID   AA03518	standard; DNA; FUN; 237 BP.
		XX
		AC   U03518;
		XX
		DE   Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
		DE   rRNA and 5.8S rRNA genes, partial sequence.
		XX
		SQ   Sequence 237 BP; 41 A; 77 C; 67 G; 52 T; 0 other;
		AA03518  Length: 237  Check: 4514  ..

			   1  aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			  61  tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			 121  ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			 181  tgagttgatt gaatg'''
			seqs = seqs.strip().split('\n\n')
##			print '/'*30
##			print 'len = ', len(seqs)
##			print seqs[0]
##			print '/'*30
##			print seqs[1]
##			print '/'*30
			self.seq1.load(seqs[0]+'\n'+seqs[1], format= 'gcg')			
			self.seq2.load(seqs[2]+'\n'+seqs[3], format= 'gcg')
			

		def load_genbank(seqs):
			'''
		LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
		DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
					rRNA and 5.8S rRNA genes, partial sequence.
		ACCESSION   U03518
		BASE COUNT	   41 a	 77 c	 67 g	 52 t
		ORIGIN	  
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc

		LOCUS	   AAU03518	  237 bp	DNA			 PLN	   04-FEB-1995
		DEFINITION  Aspergillus awamori internal transcribed spacer 1 (ITS1) and 18S
					rRNA and 5.8S rRNA genes, partial sequence.
		ACCESSION   U03518
		BASE COUNT	   41 a	 77 c	 67 g	 52 t
		ORIGIN	  
				1 aacctgcgga aggatcatta ccgagtgcgg gtcctttggg cccaacctcc catccgtgtc
			   61 tattgtaccc tgttgcttcg gcgggcccgc cgcttgtcgg ccgccggggg ggcgcctctg
			  121 ccccccgggc ccgtgcccgc cggagacccc aacacgaaca ctgtctgaaa gcgtgcagtc
			  181 tgagttgatt gaatgcaatc agttaaaact ttcaacaatg gatctcttgg ttccggc
			  '''
			seqs = seqs.strip().split('\n\n')
			self.seq1.load(seqs[0], format= 'genBank')			
			self.seq2.load(seqs[1], format= 'genBank')


		def load_ig(seqs):
			'''
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		'''
			seqs = seqs.strip().split('\n\n')
			
##			print seqs
##			print 'len = ', len(seqs)
			self.seq1.load(seqs[0], format= 'ig')			
			self.seq2.load(seqs[1], format= 'ig')

		def load_10aas(seqs):
			'''
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		
		; comment
		; comment
		U03518
		AACCTGCGGAAGGATCATTACCGAGTGCGGGTCCTTTGGGCCCAACCTCCCATCCGTGTCTATTGTACCC
		TGTTGCTTCGGCGGGCCCGCCGCTTGTCGGCCGCCGG
		'''
			seqs = seqs.strip().split('\n\n')
			
##			print seqs
##			print 'len = ', len(seqs)
			self.seq1.load(seqs[0], format= '10aas')			
			self.seq2.load(seqs[1], format= '10aas')

		def raiseSeqParsingError(seqsCount, divider):
			raise ValueError, '''Error: failed to parse input sequence pair
			into 2 sequences (you got %i !!). Something wrong in the .load()
			method of class c_2naSeqs; could be wrong value of the argument
			'divider' (you used "%s")''' % (seqsCount, divider)

		def removeExtraCodons():
			''' Compare the 2 seqs and if one of them is longer, delete
			from behind
			'''
			cd1 = self.seq1.codons
			cd2 = self.seq2.codons
			if len(cd1) > len(cd2):
				self.seq1.codons = cd1[:len(cd2)]
			elif len(cd1) < len(cd2):
				self.seq2.codons = cd2[:len(cd1)]
				
		def removeStopCodons():
			''' This remove the stopCodons from BOTH sequences.
				s1 = x   x stop x x x
				s2 = x stop x   x x x
				Then, both codon 1 and 2 will be removed from
				BOTH seqs

				This may not be what we want, but let it be for now. 
			'''
			s1 = self.seq1
			s2 = self.seq2
			sca1 = s1.stopCodonsAt()
			sca2 = s2.stopCodonsAt()
			for x in sca2:
				if x not in sca1: sca1.append(x)
			sca1.sort()
			sca1.reverse()
			for i in sca1:
				del self.seq1.codons[i]
				del self.seq2.codons[i]
##			print self.seq1.checkAttrs()
##			print self.seq2.checkAttrs()

		if format in ['fasta', 'mega', 'embl', 'gcg', 'genbank', 'ig',
					  'plain', '10aas']:
			eval('load_%s(seqs)' %format)
			removeStopCodons()
			removeExtraCodons()
		else:
			raise ValueError, '''The argument "format" in the function
			"load()" must be one of : ['fasta', 'mega', 'embl', 'gcg',
			'genbank', 'ig', 'plain', '10aas']'''	  

   
	

	def loadData(self, seqs=['', '']):
		
		if seqs[0]:
			self.__dict__['seq1']= c_naSeq(seqs[0],self.geneCodeType)
		if seqs[1]:
			self.__dict__['seq2']= c_naSeq(seqs[1],self.geneCodeType)
		return self

	def loadFile(self, file):
		dataText='\n'.join([ x.strip() for x in open(file,'r').readlines() \
						if x.strip()])
		seqPair = dataText.replace('\n','').split('seq.1')[1].split('seq.2')
		self.loadData(seqPair)
		return self
					
	def checkAttrs(self, br='\n'):
		p= self.paraDict().copy()	#{'S':S, 'V':V, 'L':L, 'P':P, 'Q':Q,
							  # 'a':a, 'b':b, 'c':c, 'A':A, 'B':B,
							  # 'K':K, 'vA':vA, 'vB':vB, 'vK':vK }
		KaKs = self.getKaKs() # {'KaByKs':..,'Ka':..., 'Ks':..,
							  #  'vKa':.., 'vKs':..}
		t=[br+ '-------------------------------------']
		t.append('geneCodeType  : %s' % self.geneCodeType)
		t.append('treatDegen3As: %s' % self.treatDegen3As)
		t.append('argCorrection : %s' % self.argCorrection)
		t.append('seq1()	 : %s' % self.seq1())
		t.append('seq2()	 : %s' % self.seq2())
		t.append('seq1.codons: %s' % self.seq1.codons)
		t.append('seq2.codons: %s' % self.seq2.codons)
		t.append('seq1 size  : %i' % len(self.seq1()))
		t.append('seq1 codons: %i' % len(self.seq1.codons))
		t.append('seq2 size  : %i' % len(self.seq2()))
		t.append('seq2 codons: %i' % len(self.seq2.codons))
		t.append('seq1.aas   : %s' % self.seq1.aas)
		t.append('seq2.aas   : %s' % self.seq2.aas)
		t.append('seq1.stopCodonsAt:%s' %self.seq1.stopCodonsAt())
		t.append('seq2.stopCodonsAt:%s' %self.seq2.stopCodonsAt())
		t.append(br+ 'degenLists():')
		t.append(' [%s,' % self.degenLists()[0])
		t.append('  %s]' % self.degenLists()[1])
		t.append(br+ 'transList():%s%s' % (br, self.transList()))
		t.append(br+ 'displayAlign():')
		t.extend(self.displayAlign())
		t.append('-------------------------------------' )
		for k in ['P', 'Q', 'a', 'b', 'c', 'A', 'B', 'K',
				  'vA', 'vB', 'vK']:
			p[k] = [(type(x)==str) and '"%s"' %x or '%1.5f' %x for x in p[k]]
			p[k] = str(p[k]).replace("'", "")
			
		t.append('S : ' + str(p['S'])) 
		t.append('V : ' + str(p['V'])) 
		t.append('L : ' + str(p['L'])) 
		t.append('P : ' + str(p['P'])) 
		t.append('Q : ' + str(p['Q']) ) 
		t.append('a : ' + str(p['a']) ) 
		t.append('b : ' + str(p['b']) ) 
		t.append('c : ' + str(p['c']) ) 
		t.append('A : ' + str(p['A']) ) 
		t.append('B : ' + str(p['B']) ) 
		t.append('K : ' + str(p['K']) ) 
		t.append('vA: ' + str(p['vA']) ) 
		t.append('vB: ' + str(p['vB']) ) 
		t.append('vK: ' + str(p['vK']) )
		t.append('-------------------------------------')		
		t.append('Ka	:' + str(KaKs['Ka']) )
		t.append('Ks	:' + str(KaKs['Ks']) )
		t.append('vKa   :' + str(KaKs['vKa']))
		t.append('vKs   :' + str(KaKs['vKs']))
		t.append('KaByKs:' + str(KaKs['KaByKs']) )
		t.append('-------------------------------------' )
		return br.join(t)
		


class c_KaKsCalculator(c_2naSeqs):
	def __init__(self,
				 seq1='',
				 seq2='',
				 geneCodeType  = 'universal',
				 treatDegen3As = 2,
				 argCorrection = 1):

		c_2naSeqs.__init__(self, seq1, seq2,
						   geneCodeType,
						   treatDegen3As,
						   argCorrection)
	def __call__(self, formula={}):
		if not formula:
			formula=Li85
			
		para = self.paraDict()
		KaKs = self.getKaKs(formula=formula)
		for k in ['Ka', 'Ks', 'vKa', 'vKs', 'KaByKs']:
			para[k] = KaKs[k]
		for k in ['Ka', 'Ks', 'vKa', 'vKs']:
			para['fml_'+k] = formula[k]
		return para

#
# This is a hacked function for calling KaKsTools for outside.
# 
# @param fasta   nucleotide sequence for the pair
# @param gct     gene code type, universal [default] or mitochondrial
# @param d3a     treat degenerate 3 A as 2 [default] or 3
# @param ag      do arg correction [1,default] or not [0]
# @param formula li93 [default] or li85
# @param call    method call from another class [1] or by command line [0]
#
def calculate(fasta,gct="universal",d3a=2,ag=1,formula="li93",call=1):
	
	kc = c_KaKsCalculator(geneCodeType=gct,treatDegen3As=d3a,argCorrection=ag)	
	
	#######################################
	# get rid of gaps in sequences
	# if modified seq is less than 80%
	# of the original, will not calculate
	#######################################
	
	def rmlb(astr):
		if astr[-2:] == "\r\n":
			astr = astr[:-2]
		elif astr[-1] == "\n":
			astr = astr[:-1]
		return astr
	
	if fasta.find("-") != -1:
		#print "   gap found..."	
		seqs = fasta.split("\n")
		seq1 = seqs[1]
		seq2 = seqs[3]		
		
		L1 = len(seq1)
		L2 = len(seq2)
		seq2m = seq1m = ""
		for i in range(len(seq1)):
			if seq1[i] == "-":
				seq2m += "-"
			else:
				seq2m += seq2[i]
		for i in range(len(seq2)):
			if seq2[i] == "-":
				seq1m += "-"
			else:
				seq1m += seq1[i]
		
		seq1 = string.joinfields(seq1m.split("-"),"")
		seq2 = string.joinfields(seq2m.split("-"),"")
		
		if(float(len(seq1))/float(L1) < 0.8 or
		   float(len(seq2))/float(L2) < 0.8):
			return ["gap","gap","gap","gap"]
		
		fasta = ">seq1\n%s\n>seq2\n%s\n" % (seq1,seq2)
	
	kc.load(fasta, format='fasta')
		
	if formula == "li85":
		formula = Li85
	else:
		formula = Li93
	
	re = kc(formula=formula)
	if call:
		return [re["Ka"],re["Ks"],re["vKa"],re["vKs"]]
	else:
		print " Ka:",re["Ka"]
		print " Ks:",re["Ks"]
		print "vKa:",re["vKa"]
		print "vKs:",re["vKs"]
	


#=============================================================
#
# Execution 
#
#=============================================================

if __name__=='__main__':
		
	fasta = """>SEQ1
ATGATCACATTTGCCAACGCACAACAAGTACCCACATTCGGTTCCGGA
CTAACATGGCCACTTATCGTCACTGCCGCCCACTTAGCCAGCTCCGCA
GTAACGCTACCTGATAGCCCACTCGTCCCCTTATTAATTGCCGTTATT
ACAATCACCCTTGAATGCGTTGTCATGAATAAATATTCTCGATTCATC
GCCACTCCATTCAATCTCCAGCGTTTAATCATATTCCCAGTCTCGTTC
CCCCTAATCAGCCGTCCTCACATCTGCTTAATATTATATACAAAATGC
GTAATAGTACTA
>SEQ2
ATCATCACATTTACAACCGCACAACAAGTACCCACATTCGGTTCCGGA
CTAACATGCCCACTTATCGTCACCACCGCCCACCTAGCAAGCTCCCCA
GTAACCCTATCTGAAAGTCCACTAGTCCTCTTATTAATCGCCGTCATA
ACAATCACCCTTGAATCCGTTGCCATTAACAAATATTCTCGATTCATC
GTCACTCCAATCAACCTCCAGCGTCTATTACTATTTCCGGTCTCGATA
CCCCTTATCAGCCGTCCTCACCTCTGCTTAATCTTATATACAAATTGC
GTATTAGTATTA
	"""
	
	aTime = time.time()
	calculate(fasta,call=0)