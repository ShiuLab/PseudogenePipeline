# Pseudogene pipeline
Scripts and wrapper file for runing the Shiu Lab's pseduogene pipeline. 

# Overview
1. Requirements 
  * python (2.7 or later, not 3; https://www.python.org/downloads/)
  * python scripts in the _pipeline_scripts folder
  * RepeatMaster (http://www.repeatmasker.org/)
  * tfasty (part of the FASTA package, version 36 or later; http://faculty.virginia.edu/wrpearson/fasta/)
    Note: The included version of tfasty was complied using "make -f ../make/Makefile.linux_sse2 all"

2. Useage

  The pipline is run using:  
  * _wrapper_scripts/CombinedPseudoWrapper.py [parameter_file]

  An example parameter_file can found in in the _example_files folder

  Additionally, a test case using chromosome 4 from A. thaliana can be run from
  the _testcase folder using the following command from within _testcase. You will
  need to supply the path to your local FASTA install and the name of the tfasty
  program you are using in test_parameter_file.txt.
  * python ../_wrapper_scripts/CombinedPseudoWrapper.py test_parameter_file.txt

  The test case takes about 15-30 minutes to run. Expected results of this run
  can be found in the _expected_results subfolder

3. Ouput

  The output of the pipeline is seperated into the following subfolder:
  
  * _intermediate: Intermediate files used to generate final results. These may be removed following a successful run, however if you want a list of pseudogenes generated prior to high confidence filtering and/or RepatMasker filtering they will be here
  * _logs: log files or the current rung
  * _results: Output files for the final list of pseudogenes following high confidence and and RepeatMasker filtering
    1. hiConf.RMfilt/hiConf.RMfilt.cdnm - position information for pseudogenes
    2. fa.hiConf.RMfilt/fa.hiConf.RMfilt.cdnm - sequence information for pseudogenes
    3. hiConf.RMfilt.cdnm.gff - gff file with pseudogene annotations
      Note: cdnm versions of the output use simplified pseudogene names.
