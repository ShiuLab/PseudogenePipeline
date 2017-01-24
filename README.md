# Pseudogene pipeline
Scripts, programs and wrapper file needed to run the pseduogene pipeline. 
This version is designed to work with tfasty36

# Useage
The pipline is run using _wrapper_scripts/CombinedPseudoWrapper.py. 
  * An example parameter file and BLOSUM matrix can be found in _example_files
  * All the required Python scripts can be foudn in _pipeline_functions
  * A precomplied version of tfasty36 is included.
    * Note: The included version of tfasty was complied using "make -f ../make/Makefile.linux_sse2 all"
  * The pipeline also requires RepeatMasker (http://www.repeatmasker.org/), which is not included.
