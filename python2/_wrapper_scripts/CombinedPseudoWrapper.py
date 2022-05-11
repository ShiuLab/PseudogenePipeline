###############################################################################
#
# 05/10/22 Shinhan - Convert code to Python 3 and release 2.0.0.
# 05/10/22 Shinhan - Create release v1.0.0 based on Nick's modification
# 09/11/13 Nick - Started to combine with other elements of the Pseudogene
#          pipeline. Final Product will incporporate the PseudoWrap script
#          Gaurav's Overlap Pipline, and the PostProcessin Wrapper. Also
#          added script to parse .gff files with wrapper and adjust gene/proteins
#          names in the case of truncation due to spaces.
#
# 02/16/12 Shinhan - Bugs... When file is empty, there is no warning. Also,
#          fast3 version change so the program name is changed.
#
# 11/21/11 Shinhan - Modify the pipeline script so it can better be used for
#          standalone applications. Make p_codes, o_codes, f_dir, and blosum
#          as required parameters. Also, activate the blast filtering part of
#          the codes.
#
# 10/14/10 Shinhan - Pseudogene pipeline wrapper. Assuming the BLAST run is
#          finished. Require a parameter file to run.
#
#

# FUNCTIONS
def GFFto4col(gff,feature,f_index,c_index,start,stop,name,truncate):
    """
    Creates a 4 column file from of selected features from a GFF file and return the name of the 4col file
    
    Inputs
        gff      := name of GFF file
        feature  := type of feature you are looking for
        f_index  := column position of features in GFF file
        c_index  := column position of chromosome in GFF file
        start    := column position of start site in GFF file
        stop     := position of stop site in GFF file
        name     := position of information from which the sequence name can be derived
        truncate := determines whether or not to truncate gene-names because of spaces
        NOTE: One is substracted from each index b/c python begin indexing at 0
    Outputs
        returns "gene_4col" the name of the 4 column file which is written to the current folder
    """

    gff_source = open(gff,'r')
    gene_4col = gff+".gene4col"
    gff_output = open(gene_4col,'w')
    for line in gff_source:
        if not line.startswith("#"):
            split_line = line.strip().split("\t")
            if split_line[f_index-1] == feature:
                info_dict = {}
                if truncate == "true":
                    sequence_information = split_line[name-1]
                    split_sequence_information = sequence_information.strip().split(";")
                    split_sequence_information = filter(None,split_sequence_information) # remove empty spaces
                    for item in split_sequence_information:
                        split_item = item.strip().split("=")
                        info_dict[split_item[0]] = split_item[1]
                    gene_name = info_dict["ID"]
                    gene_name = gene_name.strip().split(" ")[0]
                    chr_name = split_line[c_index-1]
                    chr_name = chr_name.strip().split(" ")[0]
                    newline = gene_name + "\t" + chr_name + "\t" + split_line[start-1] + "\t" + split_line[stop-1] + "\n"
                else:
                    sequence_information = split_line[name-1]
                    split_sequence_information = sequence_information.strip().split(";")
                    split_sequence_information = filter(None,split_sequence_information) # remove empty spaces
                    for item in split_sequence_information:
                        split_item = item.strip().split("=")
                        info_dict[split_item[0]] = split_item[1]
                    gene_name = info_dict["ID"]
                    newline = gene_name + "\t" + split_line[c_index-1] + "\t" + split_line[start-1] + "\t" + split_line[stop-1] + "\n"
                gff_output.write(newline)
    gff_source.close()
    gff_output.close()
    return gene_4col

def FormatTime(seconds):
    """
    Return a formatted string of the current time from the number of seconds that have passed since starting the program
    Input
         seconds := number of seconds (use time.time() - start_imter)
    Output
         current_time := formated string of time since program start
    """

    m,s = divmod(seconds,60)
    h,m = divmod(m,60)
    current_time = "%d:%02d:%02d" % (h, m, s)
    return current_time

def percentile(N, percent, key=lambda x:x):
    """
    From: http://code.activestate.com/recipes/511478/ (r1)

    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

# IMPORTS
import sys, os, time, math 

# MAIN
def main():

    #########
    ### Setup
    ########

    # Prepare timer  
    start_timer = time.time()
    print "Begin Execution  [" + FormatTime(time.time() - start_timer) + "]"

    #########
    ### Read and Check Parameter Files
    #########

    ###
    # Help Function 
    ###

    # Print help if asked for or if no argument given
    if len(sys.argv) != 2 or "--help" in sys.argv:
        print "Usage: pseudo_wrap.py parameter_file"
        print " Parameter file includes the following required values:"
        print "   b_out     Blast output (protein vs. genome)"
        print "   p_seq     Query protein sequence file"
        print "   g_seq     Genome/DNA sequence"
        print "   b_data    The type of tblastn data: tabular, blastall, or blastplus"
        print "   b_filter  Filter blast output (y) or not (n)"
        print "   p_codes   Path to pseudo pipeline scripts"
        print "   o_codes   Path to blast parser script"
        print "   f_prog    Fasta 3.x executable for tfasty"
        print "   blosum    Path to BLOSUM matrix"
        print "   gff       Maker format GFF file"
        print "   repCut    RepeatMasker cutoff parameter"
        print "   repDiv    RepeatMasker divergence parameter"
        print " And the following optional values"
        print "   ev_t      Evalue threshold in -log(E) format, default 5"
        print "   id_t      Identity threshold in %, default 40"
        print "   ml_t      Matching region length threshold in aa, default 30"
        print "   ml_p      Matching region proprotion to query threshold, default 0.05"
        print "   il_t      95th percentile intron length (skips calculating from gff file" 
        print "   "
        print " Parameter file example - /_example_files/parameter_file"
        sys.exit(0)

    ###
    # Read and check parameters
    ###
    print "Check parameter file...  [" + FormatTime(time.time() - start_timer) + "]"
    inp = open(sys.argv[1])
    inl = inp.readlines()
    P   = {}				# Paramters Dictionary
    for i in inl:
        L = i.strip().split("=")
        if len(L) != 2:
            print "Paramter file form problem in: " + i  + " [" + FormatTime(time.time() - start_timer) + "]"
            #print "Parameter file format problem:",[L]
            #print "Quit!"
            sys.exit(0)

        P[L[0]] = L[1]

    # Check if required values are missing
    required = ["b_out","p_seq","g_seq","p_codes","b_data","b_filter","o_codes","f_dir","f_prog","blosum","gff","repCut","repDiv"]
    p_missing = 0
    for i in required:
        if i not in P:
            print "Missing: " + i + " [" + FormatTime(time.time() - start_timer) + "]"
            # print " missing:",i
            p_missing = 1

    if p_missing:
        sys.exit(0)

    # Check if blast output and sequence files exist
    f_missing = 0
    r_files = [P["b_out"],P["p_seq"],P["g_seq"]]
    for i in r_files:
        if not os.path.isfile(i):
            print "File does not exist: " + i + " [" + FormatTime(time.time() - start_timer) + "]"
            #print " file does not exist:",i
            f_missing = 1

    if f_missing:
        sys.exit(0)

    if not P["b_data"] == "tabular":

        # If not not tabular, parse the raw tblastn data
        print "Process tblastn to tabular format...  [" + FormatTime(time.time() - start_timer) + "]"
        old_tblastn = P["b_out"]
        if P["b_data"] == "blast":
            os.system("python %s/ParseBlast.py -f parse_align2 -blast %s" % (P["p_codes"],old_tblastn))
        # For BLAST/blastall
        elif P["b_data"] == "blastall":
        # For BLAST+/tblastn
            os.system("python %s/ParseBlast_ModPanchy.py -f parse_align2_blastall -blast %s" % (P["p_codes"],old_tblastn)) 
        elif P["b_data"] == "blastplus":
            os.system("python %s/ParseBlast_ModPanchy.py -f parse_align2_blastPlus -blast %s" % (P["p_codes"],old_tblastn))
        P["b_out"] = old_tblastn + ".mod"

    # Check that names match between genome seq, protein seq, and blast output 
    mismatch_status = "none"
    
    os.system("grep '>' %s > %s.names" % (P["p_seq"],P["p_seq"]))
    os.system("grep '>' %s > %s.names" % (P["g_seq"],P["g_seq"]))
    genome_names_source = open(P["g_seq"]+".names",'r')
    protein_names_source = open(P["p_seq"]+".names",'r')
    genome_names_lines = genome_names_source.readlines()
    protein_names_lines = protein_names_source.readlines() 
    genome_names_list = [line.strip().strip(">") for line in genome_names_lines]
    protein_names_list = [line.strip().strip(">") for line in protein_names_lines]
    
    blast_output_source = open(P["b_out"],'r')
    blast_output_lines = blast_output_source.readlines()
    blast_output_proteinNames = [line.strip().split("\t")[0] for line in blast_output_lines]
    blast_output_contigNames = [line.strip().split("\t")[1] for line in blast_output_lines]
    
    genome_names_set = set(genome_names_list)
    protein_names_set = set(protein_names_list)

    blast_output_proteinNames_set = set(blast_output_proteinNames)
    blast_output_contigNames_set = set (blast_output_contigNames)
    protein_coverage = float(len(blast_output_proteinNames_set.intersection(protein_names_set)))/float(len(blast_output_proteinNames_set))
    contig_converage = float(len(blast_output_contigNames_set.intersection(genome_names_set)))/float(len(blast_output_contigNames_set))

    original_protein = {}
    original_contig = {}
    if protein_coverage < 0.99 or contig_converage < 0.99: # If there are missing names in the blast ouput
        print "WARNING: Blast output is not fully covered by given genome/protein sequence [" + FormatTime(time.time() - start_timer) + "]"
        print "Protein coverage: " + str(protein_coverage)
        print "Contig coverage: " + str(contig_converage)
        protein_spaces = sum([name.count(" ") for name in protein_names_list]) 
        contig_spaces = sum([name.count(" ") for name in genome_names_list])
        protein_spaces_coverage = float(protein_spaces)/float(len(protein_names_list))
        contig_spaces_coverage = float(contig_spaces)/float(len(genome_names_list))
        if protein_spaces_coverage > 0.01 or contig_spaces_coverage > 0.01: # If the coverage mismatch is caused by spaces in names
            print "The coverage problem appears to be due to the use of space characters in protein and contig names: BLAST truncates protein and contig names containg space."
            print "Temporary copies of the genome and protein sequence file will be made using truncated names. The original names will be restored in the final ouput GFF file."
            for protein_name in protein_names_list:
                original_protein[protein_name.strip().split(" ")[0]] = protein_name
            for contig_name in genome_names_list:
                original_contig[contig_name.strip().split(" ")[0]] = contig_name
            os.system("python %s/TruncateGeneNames.py %s" % (P["p_codes"],P["p_seq"]))
            os.system("python %s/TruncateGeneNames.py %s" % (P["p_codes"],P["g_seq"]))
            mismatch_status = "truncated"
            P["p_seq"] = P["p_seq"]+".truncated"
            P["g_seq"] = P["g_seq"]+".truncated"
        else: # Otherwise....
            print "The coverage problem does not appear to be due to any known issues. Pseudogene mapping will continue, but all or some pseudo-gene canidates may be omitted."
            mismatch_status = "unknown"

    # Get Intron Length Distribution
    print "Determining intron length [" + FormatTime(time.time() - start_timer) + "]"
    if "il_t" in P.keys():
        print "Using provided 95th percentile intron length"
    else:
        print "No length given. Calculating intron length distribution"
        os.system("python %s/GFFUtil.py -f checkfields -gff %s -col 2 > GFFFeatures.out" % (P["p_codes"], P["gff"]))
        feature_source = open("GFFFeatures.out",'r')
        feature_lines = feature_source.readlines()
        features = [line.strip().split(" ")[0] for line in feature_lines if "10k" not in line]
        if "intron" in features: # If there are features in the provided gff file
            print "Using introns present in GFF file"
            os.system("python %s/GFFUtil.py -f getfeat -gff %s -feat intron -col 2" % (P["p_codes"], P["gff"]))
            intron_file = P["gff"] + "_col2_intron"   
            intron_source = open(intron_file,'r')
            intron_sizes = [abs(int(line.strip().split("\t")[3]) - int(line.strip().split("\t")[4])) for line in intron_source if not line.startswith("#")] 
        else:
            print "No introns present in the GFF file, calculating intron size using gene and CDS information"
            os.system("python %s/Seq_removeGFFaltSplicing.py %s %s.noAlt" % (P["p_codes"], P["gff"], P["gff"]))
            os.system("python %s/Seq_GFFtoIntron4col_v2.py %s %s.introns" % (P["p_codes"], P["gff"], P["gff"]))
            intron_file = P["gff"] + ".introns"
            intron_source = open(intron_file,'r')
            intron_sizes = [abs(int(line.strip().split("\t")[2]) - int(line.strip().split("\t")[3])) for line in intron_source if not line.startswith("#")]
        intron_sizes.sort() # Sort the intron_sizes list
        P["il_t"] = int(math.ceil(percentile(intron_sizes,.95))) 
        os.system("rm %s" % (intron_file)) 
    print "95th percentile intron legnth is " + str(P["il_t"]) + " [" + FormatTime(time.time() - start_timer) + "]"

    p_codes = P["p_codes"]
    o_codes = P["o_codes"]
    f_dir   = P["f_dir"]
    f_prog  = P["f_prog"]
    blosum  = P["blosum"]

    # Check if other values are set, if not default is used.
    # Set default values
    default = {"ev_t"    :"5",
                                         "id_t"    :"40",
                                         "ml_t"    :"30",
                                         "ml_p"    :"0.05"}
    for i in default:
        if i not in P:
            P[i] = default[i]
            print " default: %s=%s" % (i,P[i])

    #########
    ### Pseudogene Identification
    #########

    print "Begin Candidate Pseudogene Identifiction... [" + FormatTime(time.time() - start_timer) + "]\n"

    ###
    # Filter BLAST output
    ###


    f_flag = 0
    if "b_filter" in P and P["b_filter"] == "y":
        f_flag = 1
        try:
            #print "Filter %s..." % P["b_out"]
            #print " E:%s I:%s L:%s P:%s" % (P["ev_t"],P["id_t"] ,P["ml_t"] ,P["ml_p"])
            print "Filter " + P["b_out"] + "... [" + FormatTime(time.time() - start_timer) + "]"
            print " E:" + P["ev_t"] + " I:" + P["id_t"] + " L:" + P["ml_t"] + " P:" + P["ml_t"]
        except KeyError:
            #print "Missing filtering parameter(s)! Quit!"
            print "Missing filtering parameter(s)! [" + FormatTime(time.time() - start_timer) + "]"
            sys.exit(0)
        os.system("python %s/ParseBlast.py -f get_qualified4 -blast " % o_codes +\
                                          "%s -fasta %s -E %s -I " % (P["b_out"],P["p_seq"],P["ev_t"])  +\
                                          "%s -L %s -P %s -Q 1 "   % (P["id_t"] ,P["ml_t"] ,P["ml_p"])  +\
                                          "> log_step1")

        ml_p = int(float(P["ml_p"])*100)
        os.system("mv %s_E%sI%sL%sP%sQ1.qlines %s_parsed" % \
                                                (P["b_out"],P["ev_t"],P["id_t"],P["ml_t"],ml_p,P["b_out"]))
        P["b_out"] = "%s_parsed" % P["b_out"]

        # Check if there is problem with the file

    ###
    # Get pseudoexon and phase 1 pseudogene
    ###

    print "Get pseudoexons... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/script_step2e.py %s %s > log_step2_pe" % (p_codes,P["b_out"],P["il_t"]))
    pe = "%s_G%s.PE" % (P["b_out"],P["il_t"])
    print "Pseudoexon file:" + pe + " [" + FormatTime(time.time() - start_timer) + "]"

    print "Get phase 1 pseudogene... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/script_step3bDave.py %s %s > log_step3_phase1_ps" % (p_codes,pe,P["il_t"]))
    ps1 = "%s_I%s.PS1" % (pe,P["il_t"])
    print "Phase1 file:" + ps1 + " [" + FormatTime(time.time() - start_timer) + "]"

    ###
    # Prepare ps1 sequnce file and pair file
    ###

    print "Get pair file and subject coordinates... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/script_step3.5.py %s" % (p_codes,ps1))
    print "Pair file:" + ps1 + ".pairs [" + FormatTime(time.time() - start_timer) + "]"
    print "Coordinates file:" + ps1 + ".subj_coord [" + FormatTime(time.time() - start_timer) + "]"

    # Get ps1 sequences
    print "Get phase 1 pseudogene sequences... [" + FormatTime(time.time() - start_timer) + "]"
    ps1coord = ps1 + ".subj_coord"
    os.system("python %s/FastaManager.py -f get_stretch4 -fasta " % o_codes + \
                  "%s -coords %s > log_step4_phase1_ps_seq" % (P["g_seq"],ps1coord))
    print "Phase 1 Sequence file:" + ps1 + ".subj_coord.fa [" + FormatTime(time.time() - start_timer) + "]"

    ###
    # Run SW
    ###

    print "Find stop and framshifts... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/BlastUtility.py -f batch_sw -p %s -g " % (o_codes,f_prog)+\
                  "%s.pairs -i %s -j %s.fa -bdir " % (ps1,P["p_seq"],ps1coord)    + \
                  "%s -d 1"         % f_dir)
    print "Smith-Waterman outputs:" + ps1 + "_pairs.sw* [" + FormatTime(time.time() - start_timer) + "]"

    os.system("python %s/script_step6.py %s %s.pairs.sw.out" % (p_codes,blosum,ps1))

    print "Disable_count file:" + ps1 + ".pairs.sw.out.disable_count [" + FormatTime(time.time() - start_timer) + "]"

    ###
    # End of section
    ###
    print "Finished Candidate Pseudogene Identifiction... [" + FormatTime(time.time() - start_timer) + "]\n"

    #########
    ### Find Overlaps
    #########

    print "Begin Overlap Pipeline... [" + FormatTime(time.time() - start_timer) + "]"

    # Set file prefix
    base = P["b_out"]

    # Process disable_count file
    print "Process Disable_counts file... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/1_3_annotatePseudogenes.py %s.pairs.sw.out.disable_count %s.4col" %  (p_codes,ps1,base))
    print "Processed Disable_count file:" + base + ".4col [" + FormatTime(time.time() - start_timer) + "]"

    #Process Maker format GFF into a gene 4 column file
    #Note: The numbers in the following string are the standard index position for the feature, chromosome, start poistion, stop postion, and seuqence inforamtion in GFF/GFF3 files 
    if mismatch_status == "truncated":
        genes4col = GFFto4col(P["gff"],"gene",3,1,4,5,9,"true")
    else:
        genes4col = GFFto4col(P["gff"],"gene",3,1,4,5,9,"false")
          
    # Find Overlaps
    print "Find Overlaps... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/v4_blastmap_to_gff_part2_smallchange.py %s %s.4col n"  % (p_codes,genes4col,base))
    print "Overlap file:" + base + ".4col.onlyoverlap [" + FormatTime(time.time() - start_timer) + "]"

    ###
    # End of section
    ###

    print "Finished Overlap Pipeline... [" + FormatTime(time.time() - start_timer) + "]\n"

    #########
    ### Pseudogene Post-Processing
    #########

    print "Begin Post-Processing Pipeline... [" + FormatTime(time.time() - start_timer) + "]"

    # The following section is adapted from seperate script, so some of the variable names need to be changed
    fasta = P['g_seq'] 
    prot = P['p_seq']
    gff = P['gff']
    intronLen = P['il_t']
    scriptDir = P['p_codes']
    repCut = P['repCut']
    repDiv = P['repDiv']

    # Filter out pseudogene predictions overlapping with genes
    #Note: You will always need to install/load the Repeat Masker ahead of time
    print "Filter Pseudogenes Overlaping with Genes... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/1_3_addTruePseudogenes2.py %s.4col %s.4col.onlyoverlap %s.4col.true" % (scriptDir, base, base, base))
    os.system("python %s/1_3_pseuOrigToTrue.py %s.pairs.sw.out.disable_count %s.4col.true %s.fullyFiltered.disable_count" % (scriptDir, ps1, base, base))
    #to not filter use these lines instead
    #os.system("cp %s.4col %s.4col.true" % (base, base))
    #os.system("cp %s.tblastn_parsed_G500.PE_I500.PS1.pairs.sw.out.disable_count %s.fullyFiltered.disable_count" % (base, base))
    print "Filitered File: " + base + ".fullyFiltered.disable_count  [" + FormatTime(time.time() - start_timer) + "]"

    # Repeat masking prep
    os.system("python %s/FastaManagerGaurav.py -f get_stretch4 -fasta %s -coords %s.4col.true " % (scriptDir, fasta, base))
    os.system("python %s/1_3_v2_MakeFastaCodeName.py %s.4col.true.fa %s.4col.true.fa.mod %s.4col.true.fa.ref" % (scriptDir, base, base, base))

    # Repeat Masking
    print "Begin Repeat Masking... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("grep -v '#' %s.4col.true.fa.mod > %s.4col.true.fa.mod.mod" % (base, base))
    os.system("RepeatMasker -cutoff 200 -par 8 -xsmall -species viridiplantae -gff %s.4col.true.fa.mod.mod" % (base))
    os.system("python %s/parse_RepeatMasker_gff.py %s.4col.true.fa.mod.mod.out %s %s" % (scriptDir, base, repCut, repDiv))
    os.system("python %s/NameChangers/1_3_removeCodeNamesFromRM4col.py %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col %s.4col.true.fa.ref %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames RM4col 0" % (scriptDir, base, repCut, repDiv, base, base, repCut, repDiv)) #I need this just so I don't have to modify the next script (I know it's bad programming, but it was a time saver)
    os.system("grep -v Simple_repeat %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames | grep -v Low_complexity | grep -v Satellite > %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames.filt" % (base, repCut, repDiv, base, repCut, repDiv))
    os.system("python %s/1_3_filterOutPseuRMs.py %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames.filt %s.fullyFiltered.disable_count %s.4col.true %s.4col.true.fa %s.fullyFiltered.disable_count.RMfilt %s.4col.true.RMfilt %s.4col.true.fa.RMfilt" % (scriptDir, base, repCut, repDiv, base, base, base, base, base, base))
    print "Repeat Masked files end with '.RMfilt'  [" + FormatTime(time.time() - start_timer) + "]"

    # Cleanup
    #os.system("mkdir Temp")
    #os.system("mv %s.tblastn_* Temp" % (base))
    #os.system("mv %s.4col.true.fa.mod.mod* Temp" % (base))
    #os.system("mv Temp/*count .")

    # Do high confidence filter
    print "Run High-Confidence Filter... [" + FormatTime(time.time() - start_timer) + "]"

    # Get Genome and Sequence Sizes
    os.system("python %s/FastaManagerDave.py -f get_sizes -fasta %s" % (scriptDir, fasta))
    os.system("python %s/FastaManagerDave.py -f get_sizes -fasta %s" % (scriptDir, prot))

    # Run High-Confidence Filter
    os.system("python %s/1_3_v2_verifyPseudogenes.py %s.fullyFiltered.disable_count.RMfilt  %s.4col.true.RMfilt %s.4col.true.fa.RMfilt %s.size %s.size %s.fullyFiltered.disable_count.RMfilt.hiConf %s.4col.true.RMfilt.hiConf %s.4col.true.fa.RMfilt.hiConf %s" % (scriptDir, base, base, base, fasta, prot, base, base, base, intronLen))
    print "High-Confidence files end with '.hiConf'  [" + FormatTime(time.time() - start_timer) + "]"

    # Bring back to code names %s.4col.true.RMfilt %s.4col.true.fa.mod.mod.RMfilt
    print "Correct gene-names... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/TabularManager.py -f RemCodeName -inp %s.4col.true.RMfilt.hiConf -out %s.4col.true.RMfilt.hiConf.cdnm -col 1 -ref %s.4col.true.fa.ref -FR R" % (scriptDir, base, base, base))
    os.system("python %s/FastaManagerDave.py -f replace_names -fasta %s.4col.true.fa.RMfilt.hiConf -name %s.4col.true.fa.ref -out %s.4col.true.fa.RMfilt.hiConf.cdnm" % (scriptDir, base, base, base))
    print "Name corrected files end with '.cdnm'  [" + FormatTime(time.time() - start_timer) + "]"

    print "Write Pseudogene GFF... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("python %s/1_3_getPseudoGFF.py %s.4col.true.RMfilt.hiConf.cdnm %s.fullyFiltered.disable_count.RMfilt.hiConf %s.4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.gff" % (scriptDir, base, base, base, base))
    print "GFF file:" + base + ".4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.gff  [" + FormatTime(time.time() - start_timer) + "]"

    if mismatch_status == "truncated":
        print "Rewriting Pseudogene GFF with original names... [" + FormatTime(time.time() - start_timer) + "]"
        gff_source = open(base + ".4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.gff",'r')
        gff_original_out = open(base + ".4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.original.gff",'w')
        for line in gff_source:
            split_line = line.strip().split("\t")
            chromosome = split_line[0]
            sequence_information = split_line[8]
            original_chromosome = original_contig[chromosome]
            split_sequence_information = sequence_information.strip().split(";")
            derivation = split_sequence_information[1]
            split_derivation = derivation.strip().split("=")
            genename = split_derivation[1]
            original_genename = original_protein[genename]
            original_derivation = split_derivation[0] + "=" + original_genename
            original_sequence_information = split_sequence_information[0]+";"+original_derivation+";"+split_sequence_information[2]
            original_line = original_chromosome + "\t" + "\t".join(split_line[1:8]) + "\t" + original_sequence_information + "\n"
            gff_original_out.write(original_line)
        gff_source.close()
        gff_original_out.close()
        print "GFF file with non-truncated names:" + base + ".4col.true.RMfilt.hiConf.cdnm.original.gff  [" + FormatTime(time.time() - start_timer) + "]"
        
    print "Removing temporary genome and protein sequences files... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("rm %s.*" % (P["p_seq"])) 
    os.system("rm %s.*" % (P["g_seq"]))

    print "Moving intermediate files to _intermediate... [" + FormatTime(time.time() - start_timer) + "]"
    os.system("mkdir _intermediate")
    os.system("mv %s* _intermediate" % (P["b_out"]))

    print "Move hiConf, codename and gff files to _results"
    os.system("mkdir _results")
    os.system("mv _intermediate/*.hiConf _results")
    os.system("mv _intermediate/*.hiConf.cdnm _results")
    os.system("mv _intermediate/*.hiConf.cdnm.gff _results")

    print "Move log files to _logs"
    os.system("mkdir _logs")
    os.system("mv log* _logs")

    ###
    # End of section & Pipline
    ###

    print "Finished Post-Processing Pipeline... [" + FormatTime(time.time() - start_timer) + "]"
    print "End Execution  [" + FormatTime(time.time() - start_timer) + "]"

# If not loaded as a module
if __name__ == '__main__':
    main()
