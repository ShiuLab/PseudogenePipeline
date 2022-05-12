###############################################################################
#
# 05/10/22 Shinhan - Convert code to Python 3 and release 2.0.0. 
#          To do:
#           1. Call scripts by importing them as modules
#           2. Explicitly assign variables to names of input/output files
#           3. Super long lines need to be reined in. 
# 05/10/22 Shinhan - Create release v1.0.0 based on the existing version.
#          https://github.com/ShiuLab/PseudogenePipeline/releases/tag/v1.0.0
# 09/11/13 Nick - Started to combine with other elements of the Pseudogene
#          pipeline. Final Product will incporporate the PseudoWrap script
#          Gaurav's Overlap Pipline, and the PostProcessin Wrapper. Also
#          added script to parse .gff files with wrapper and adjust gene and
#          proteins names in the case of truncation due to spaces.
# 02/16/12 Shinhan - Bugs... When file is empty, there is no warning. Also,
#          fast3 version change so the program name is changed.
# 11/21/11 Shinhan - Modify the pipeline script so it can better be used for
#          standalone applications. Make p_codes, p_codes, f_dir, and blosum
#          as required parameters. Also, activate the blast filtering part of
#          the codes.
# 10/14/10 Shinhan - Pseudogene pipeline wrapper. Assuming the BLAST run is
#          finished. Require a parameter file to run.

import sys, os, time, math, subprocess, ParseBlast

def help():
    '''Help function, need to be replaced with argParse'''
    print("Usage: pseudo_wrap.py parameter_file")
    print(" Parameter file includes the following required values:")
    print("   b_out     Blast output (protein vs. genome)")
    print("   p_seq     Query protein sequence file")
    print("   g_seq     Genome/DNA sequence")
    print("   b_format  tblastn output tabular (outfmt=6) or paiwise (outfmt=0)")
    print("   b_filter  Filter blast output (y) or not (n)")
    print("   p_codes   Path to pipeline scripts")
    print("   f_dir     Path to fasta, >=v36")
    print("   f_prog    executable file name for tfasty")
    print("   blosum    Path to BLOSUM matrix")
    print("   gff       Maker format GFF file")
    print("   r_dir     RepeatMasker executable directory")
    print("   r_species RepeatMasker species")
    print("   r_cut     RepeatMasker cutoff parameter")
    print("   r_div     RepeatMasker divergence parameter")
    print(" And the following optional values")
    print("   ev_t      Evalue threshold in -log(E) format, default 5")
    print("   id_t      Identity threshold in %, default 40")
    print("   ml_t      Matching region length threshold in aa, default 30")
    print("   ml_p      Matching region-to-query proportion threshold, default") 
    print("             0.05")
    print("   il_t      95th percentile intron length, if not provided, will")
    print("             determined based on GFF or sequenece info.") 
    print("   overlap   Rid of pseudogenes overlapping with genes (1) or not")
    print("             (0, default)") 
    print("   ")
    print(" Parameter file example - /_example_files/parameter_file")
    sys.exit(0)

def GFFto4col(gff,feature,f_index,c_index,start,stop,name,truncate):
    """Creates a 4 column file from of selected features from a GFF file and 
    return the name of the 4col file
    Args
      gff      : name of GFF file
      feature  : type of feature you are looking for
      f_index  : column position of features in GFF file
      c_index  : column position of chromosome in GFF file
      start    : column position of start site in GFF file
      stop     : position of stop site in GFF file
      name     : position of information from which the seq name can be derived
      truncate :  whether or not to truncate gene-names because of spaces
      NOTE: One is substracted from each index b/c python begin indexing at 0
    Return:
      gene_4col: the name of the 4 column file written to the current folder
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
                    seq_info = split_line[name-1]
                    split_seq_info = seq_info.strip().split(";")
                    split_seq_info = [_f for _f in split_seq_info if _f] 
                    # remove empty spaces
                    for item in split_seq_info:
                        split_item = item.strip().split("=")
                        info_dict[split_item[0]] = split_item[1]
                    gene_name = info_dict["ID"]
                    gene_name = gene_name.strip().split(" ")[0]
                    chr_name = split_line[c_index-1]
                    chr_name = chr_name.strip().split(" ")[0]
                    newline = gene_name + "\t" + chr_name + "\t" + \
                        split_line[start-1] + "\t" + split_line[stop-1] + "\n"
                else:
                    seq_info = split_line[name-1]
                    split_seq_info = seq_info.strip().split(";")
                    # remove empty spaces
                    split_seq_info = [_f for _f in split_seq_info if _f] 
                    for item in split_seq_info:
                        split_item = item.strip().split("=")
                        info_dict[split_item[0]] = split_item[1]
                    gene_name = info_dict["ID"]
                    newline = gene_name + "\t" + split_line[c_index-1] + "\t" +\
                         split_line[start-1] + "\t" + split_line[stop-1] + "\n"
                gff_output.write(newline)
    gff_source.close()
    gff_output.close()
    return gene_4col

def FormatTime(seconds):
    """Return a formatted string of the current time from the number of seconds 
    that have passed since starting the program
    Input
         seconds := number of seconds (use time.time() - start_imter)
    Output
         current_time := formated string of time since program start
    """

    m,s = divmod(seconds,60)
    h,m = divmod(m,60)
    current_time = "%d:%02d:%02d" % (h, m, s)
    return current_time

def percentile(N, percentile, key=lambda x:x):
    """Find the percentile of a list of values.
    From: http://code.activestate.com/recipes/511478/ (r1)
    Args:
        N (list): values for determining percentile. N MUST BE already sorted.
        percentile (float): The value in N will be obtained based on the 
            percentile value specified here. 0.0 to 1.0.
        key (lambda): function to compute value from each element of N. Optional.
    Return:
        d0+d1 (float): value at the specified percentile.
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

def check_parameters(parameter_file):
    """"Check whether the parameter file is formated properly"""

    print("Check parameter file...")
    inp = open(parameter_file)
    inl = inp.readlines()
    P   = {} # Paramters Dictionary

    # Check if this is delimited with "=" with 2 elements each line
    p_misformat = 0
    c = 0
    for i in inl:
        # Comment line, ignore
        if i[0] == "#": continue

        L = i.strip().split("=")
        if len(L) != 2:
            print(f"  line {c} parameter format problem: {L}")
            p_misformat = 1
        else:     
            P[L[0]] = L[1]
        c += 1

    # Check if required parameter values are missing
    required = ["b_out","p_seq","g_seq","p_codes","b_format","b_filter",
                "f_dir","f_prog","blosum","gff","r_cut","r_div"]
    p_missing = 0
    for i in required:
        if i not in P:
            print(" missing:",i)
            p_missing = 1

    # Check if blast output, sequence, and GFF files exist
    f_missing = 0
    r_files = [P["b_out"], P["p_seq"], P["g_seq"], P["gff"]]
    for i in r_files:
        if not os.path.isfile(i):
            print("  required file does not exist:",i)
            f_missing = 1

    # Convert non-tabular blast+ output to tabular
    b_unkformat = 0
    if not P["b_format"] == "tabular":
        # If not not tabular, parse the raw tblastn data
        print("\nProcess tblastn to tabular format... ")
        script_dir = P["p_codes"]
        blast_out  = P["b_out"]
        if P["b_format"] == "pairwise":
            blast_parser = ParseBlast.parser()
            blast_parser.parse_align2_blastPlus(blast_out)
            P["b_out"] = blast_out + ".mod"
        else:
            print("Error: Unknown blast format -", P["b_format"])
            b_unkformat = 1

    # Set the overlap flag to 0 if not passed and check if unk flag is passed
    o_unk = 0
    if "overlap" not in P:
        P["overlap"] = "0"
    elif P["overlap"] not in ["0", "1"]:
        o_unk = 1
        print("ERR: unknown overlap flag,", P["overlap"])

    if p_misformat or p_missing or f_missing or b_unkformat or o_unk:
        sys.exit(0)

    return P

def check_seq_names(p_seq, g_seq, b_out):
    '''Check that names match between genome seq, protein seq, and blast output
    Args:
        p_seq (str): protein sequence file name
        g_seq (str): genome sequence file name
        b_out (str): tabular blast output file name
    Return:
        mismatch_status (str): a status flag
    ''' 
    mismatch_status = "none"
    
    os.system("grep '>' %s > %s.names" % (p_seq, p_seq))
    os.system("grep '>' %s > %s.names" % (g_seq, g_seq))
    g_names_source = open(g_seq+".names",'r')
    p_names_source = open(p_seq+".names",'r')
    g_names_lines = g_names_source.readlines()
    p_names_lines = p_names_source.readlines() 
    g_names_list = [line.strip().strip(">") for line in g_names_lines]
    p_names_list = [line.strip().strip(">") for line in p_names_lines]
    
    b_out_source = open(b_out,'r')
    b_out_lines = b_out_source.readlines()
    b_out_p_names = [line.strip().split("\t")[0] for line in b_out_lines]
    b_out_c_names = [line.strip().split("\t")[1] for line in b_out_lines]
    
    g_names_set = set(g_names_list)
    p_names_set = set(p_names_list)

    b_out_p_names_set = set(b_out_p_names)
    b_out_c_names_set = set (b_out_c_names)
    b_out_p_names_set_inters = b_out_p_names_set.intersection(p_names_set)
    b_out_c_names_set_inters = b_out_c_names_set.intersection(g_names_set)
    p_cov = float(len(b_out_p_names_set_inters))/float(len(b_out_p_names_set))
    c_cov = float(len(b_out_c_names_set_inters))/float(len(b_out_c_names_set))

    original_protein = {}
    original_contig = {}

    # If there are missing names in the blast ouput
    ERR_cov1 = "The coverage problem appears to be due to the use of space "+\
               "characters in protein and contig names: BLAST truncates "+\
               "protein and contig names containg space.\n"+\
               "Temporary copies of the genome and protein sequence file "+\
               "will be made using truncated names. The original names will "+\
               "be restored in the final ouput GFF file."
    ERR_cov2 = "The coverage problem does not appear to be due to any known "+\
               "issues. Pseudogene mapping will continue, but all or some "+\
               "pseudo-gene canidates may be omitted."
    if p_cov < 0.99 or c_cov < 0.99: 
        print("WARNING: Blast output not fully covered by genome/protein seq")
        print("Protein coverage: " + str(p_cov))
        print("Contig coverage: " + str(c_cov))
        protein_spaces = sum([name.count(" ") for name in p_names_list]) 
        contig_spaces = sum([name.count(" ") for name in g_names_list])
        protein_spaces_coverage = float(protein_spaces)/float(len(p_names_list))
        contig_spaces_coverage = float(contig_spaces)/float(len(g_names_list))
        if protein_spaces_coverage > 0.01 or contig_spaces_coverage > 0.01: # If the coverage mismatch is caused by spaces in names
            print(ERR_cov1)
            for protein_name in p_names_list:
                original_protein[protein_name.strip().split(" ")[0]] = protein_name
            for contig_name in g_names_list:
                original_contig[contig_name.strip().split(" ")[0]] = contig_name
            os.system(
                "python %s/TruncateGeneNames.py %s" % (P["p_codes"],P["p_seq"]))
            os.system(
                "python %s/TruncateGeneNames.py %s" % (P["p_codes"],P["g_seq"]))
            mismatch_status = "truncated"
            P["p_seq"] = P["p_seq"]+".truncated"
            P["g_seq"] = P["g_seq"]+".truncated"
        else: # Otherwise....
            print(ERR_cov2)
            mismatch_status = "unknown"

    return mismatch_status

def get_intron_len(P):
    '''Get Intron Length Distribution'''
    il_t = 0
    if "il_t" in P:
        print("  using provided 95th percentile intron length")
        il_t = P['il_t']
    else:
        print("  no length given, calculating intron length distribution")
        os.system("python %s/GFFUtil.py -f checkfields -gff %s -col 2 > GFFFeatures.out" % (P["p_codes"], P["gff"]))
        feat_src = open("GFFFeatures.out",'r')
        feat_L = feat_src.readlines()
        feats = [L.strip().split(" ")[0] for L in feat_L if "10k" not in L]
        if "intron" in features: # If there are features in the provided gff file
            print("  using introns present in GFF file")
            os.system("python %s/GFFUtil.py -f getfeat -gff %s -feat intron -col 2" % (P["p_codes"], P["gff"]))
            intron_file = P["gff"] + "_col2_intron"   
            intron_source = open(intron_file,'r')
            intron_sizes = [abs(int(L.strip().split("\t")[3]) - int(L.strip().split("\t")[4])) for L in intron_source if not L.startswith("#")] 
        else:
            print("  no intron info in GFF, get intron size from gene/CDS info")
            os.system("python %s/Seq_removeGFFaltSplicing.py %s %s.noAlt" % \
                                            (P["p_codes"], P["gff"], P["gff"]))
            os.system("python %s/Seq_GFFtoIntron4col_v2.py %s %s.introns" % \
                                            (P["p_codes"], P["gff"], P["gff"]))
            intron_file = P["gff"] + ".introns"
            intron_source = open(intron_file,'r')
            intron_sizes = [abs(int(line.strip().split("\t")[2]) - int(line.strip().split("\t")[3])) for line in intron_source if not line.startswith("#")]
        intron_sizes.sort() # Sort the intron_sizes list
        il_t = int(math.ceil(percentile(intron_sizes,.95))) 
        os.system("rm %s" % (intron_file)) 
    
    print("  95th percentile intron legnth is", il_t)

    return il_t

def main():
    '''Main function of the pipeline'''

    #########
    ## STEP 0
    print("\n#### Step 0: read, check, and set parameters\n")

    # Print help if asked for or if no argument given
    if len(sys.argv) != 2 or "--help" in sys.argv:
        help()

    # Start checking parameters
    param_file = sys.argv[1]
    P = check_parameters(param_file)
    
    print("Check names among genome/protein seq, and blast output")
    mismatch_status = check_seq_names(P["p_seq"], P["g_seq"], P["b_out"]) 

    print("Determining intron length")
    il_t = get_intron_len(P)
    P['il_t'] = il_t

    p_codes = P["p_codes"]
    f_dir   = P["f_dir"]
    f_prog  = P["f_prog"]
    blosum  = P["blosum"]

    # Check if other values are set, if not default is used.
    print("Set default values")
    default = {"ev_t"    :"5",
                "id_t"    :"40",
                "ml_t"    :"30",
                "ml_p"    :"0.05"}
    for i in default:
        if i not in P:
            P[i] = default[i]
            print("  default: %s=%s" % (i,P[i]))

    #########
    ## STEP 1
    print("\n#### Step 1. Filtering BLAST output...\n")

    f_flag = 0
    if "b_filter" in P and P["b_filter"] == "y":
        f_flag = 1
        try:
            print("Filter %s..." % P["b_out"])
            print(" E:%s I:%s L:%s P:%s" % \
                                    (P["ev_t"],P["id_t"] ,P["ml_t"] ,P["ml_p"]))
        except KeyError:
            # This should not happen
            print("ERR: missing filtering parameter(s)! Quit!")
            sys.exit(0)

        os.system("python %s/ParseBlast.py -f get_qualified4 -blast" % p_codes+\
            " %s -fasta %s -E %s -I" % (P["b_out"],P["p_seq"],P["ev_t"])  +\
            " %s -L %s -P %s -Q 1"   % (P["id_t"] ,P["ml_t"] ,P["ml_p"])  +\
            " > log_step1")

        ml_p = int(float(P["ml_p"])*100)
        
        # Rename ParseBlast output
        os.system("mv %s_E%sI%sL%sP%sQ1.qlines %s_parsed" % \
                    (P["b_out"],P["ev_t"],P["id_t"],P["ml_t"],ml_p,P["b_out"]))
        P["b_out"] = "%s_parsed" % P["b_out"]

    #########
    ## STEP 2
    print("\n#### Step 2. Get pseudoexon\n")
    os.system("python %s/script_step2e.py %s %s > log_step2_pe" % \
                                                (p_codes,P["b_out"], P["il_t"]))
    pe = "%s_G%s.PE" % (P["b_out"],P["il_t"])
    print("Pseudoexon file:", pe)

    #########
    ## STEP 3
    print("\n#### Step 3. Get phase 1 pseudogene\n")

    print("Get initial phase 1 pseudogene")
    os.system("python %s/script_step3bDave.py %s %s > log_step3_phase1_ps" % \
                                                (p_codes,pe, P["il_t"]))
    ps1 = "%s_I%s.PS1" % (pe,P["il_t"])
    print("  phase1 file:", ps1)

    print("Get pair file and subject coordinates")
    os.system("python %s/script_step3.5.py %s" % (p_codes,ps1))
    print("  pair file:", ps1 + ".pairs")
    print("  coordinates file:" + ps1 + ".subj_coord")

    print("Get phase 1 pseudogene sequences")
    ps1coord = ps1 + ".subj_coord"
    os.system("python %s/FastaManager.py -f get_stretch4 -fasta " % p_codes + \
            "%s -coords %s > log_step4_phase1_ps_seq" % (P["g_seq"],ps1coord))
    print("  phase 1 seq:", ps1 + ".subj_coord.fa")

    #########
    ## STEP 4
    print("\n#### Step 4. Get Smith-Waterman alignments with fasta\n")

    print("Find stop and framshifts")
    os.system("python %s/BlastUtility.py -f batch_sw -p %s" % (p_codes,f_prog)+\
            " -g %s.pairs -i %s -j %s.fa -bdir" % (ps1,P["p_seq"],ps1coord)   +\
            " %s -d 1"         % f_dir)
    print("  smith-Waterman outputs:", ps1 + "_pairs.sw*")

    os.system("python %s/script_step6.py %s %s.pairs.sw.out" % \
                                                        (p_codes,blosum,ps1))
    print("  disable_count file:", ps1 + ".pairs.sw.out.disable_count")

    #########
    ## Step 5 
    print("\n#### Step 5. Convert coordinates\n")

    # Set file prefix
    base = P["b_out"]

    print("Convert S-W coordinates to genome-based ones")
    os.system("python %s/1_3_annotatePseudogenes.py %s.pairs.sw.out.disable_count %s.4col" %  (p_codes,ps1,base))
    print("  processed Disable_count file:", base + ".4col")

    #Process Maker format GFF into a gene 4 column file
    #Note: The numbers in the following string are the standard index position 
    #  for the feature, chromosome, start poistion, stop postion, and seuqence 
    #  inforamtion in GFF/GFF3 files 
    print("Process Maker format GFF if necessary")
    if mismatch_status == "truncated":
        genes4col = GFFto4col(P["gff"],"gene",3,1,4,5,9,"true")
    else:
        genes4col = GFFto4col(P["gff"],"gene",3,1,4,5,9,"false")

    #########
    ## Step 6
    print("\n#### Step 6. Determine overlap with annotated genes\n")

    # The following section is adapted from seperate script, so some of the 
    # variable names follow those in the old script
    fasta = P['g_seq'] 
    prot = P['p_seq']
    gff = P['gff']
    intronLen = P['il_t']
    scriptDir = P['p_codes']
    r_cut = P['r_cut']
    r_div = P['r_div']
    r_dir   = P['r_dir']        # 5/11,22 added this option
    r_species = P['r_species']  # 5/11,22 added this option
    overlap = P["overlap"]      # 5/11/22 added this

    print("Find Overlaps between predicted pseudogenes and genes")
    os.system("python %s/v4_blastmap_to_gff_part2_smallchange.py %s %s.4col n" %\
         (p_codes,genes4col,base))
    print("  overlap file:", base + ".4col.onlyoverlap")

    if overlap == "1":
        print("Filter out pseudogene predictions overlapping with genes")
        os.system("python %s/1_3_addTruePseudogenes2.py %s.4col %s.4col.onlyoverlap %s.4col.true" % (scriptDir, base, base, base))
        os.system("python %s/1_3_pseuOrigToTrue.py %s.pairs.sw.out.disable_count %s.4col.true %s.fullyFiltered.disable_count" % (scriptDir, ps1, base, base))
    else:
        print("Do not filter pseudogene predictions overlapping with genes")
        os.system("cp %s.4col %s.4col.true" % (base, base))
        os.system("cp %s.pairs.sw.out.disable_count %s.fullyFiltered.disable_count" % (ps1, base))
    
    print("  output file:", base + ".fullyFiltered.disable_count")

    #########
    ## Step 7
    print("\n#### Step 7. Repeat masking")
    
    print("\nRepeat masking prep")
    #Note: You will always need to install/load the Repeat Masker ahead of time
    os.system("python %s/FastaManager.py -f get_stretch5 -fasta %s -coords %s.4col.true " % (scriptDir, fasta, base))
    os.system("python %s/1_3_v2_MakeFastaCodeName.py %s.4col.true.fa %s.4col.true.fa.mod %s.4col.true.fa.ref" % (scriptDir, base, base, base))

    print("\nRun RepeatMasker")
    os.system("grep -v '#' %s.4col.true.fa.mod > %s.4col.true.fa.mod.mod" % (base, base))
    os.system("%s/RepeatMasker -cutoff 200 -par 8 -xsmall -species %s -gff %s.4col.true.fa.mod.mod" % (r_dir, r_species, base))

    print("\nParse RepeatMasker output")
    os.system("python %s/parse_RepeatMasker_gff.py %s.4col.true.fa.mod.mod.out %s %s" % (scriptDir, base, r_cut, r_div))

    # SHS: Not sure whose comment this is...
    # I need this just so I don't have to modify the next script (I know it's 
    # bad programming, but it was a time saver)
    os.system("python %s/1_3_removeCodeNamesFromRM4col.py " % scriptDir+\
              "%s.4col.true.fa.mod.mod.out." % base +\
              "%sCutoff%s.0divg.4col " % (r_cut, r_div) +\
              "%s.4col.true.fa.ref " % base+\
              "%s.4col.true.fa.mod.mod.out.%sCutoff" % (base, r_cut)+\
              "%s.0divg.4col.longNames RM4col 0" % (r_div)) 

    # os.system("python %s/NameChangers/1_3_removeCodeNamesFromRM4col.py 
    # %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col %s.4col.true.fa.ref 
    # %s.4col.true.fa.mod.mod.out.%sCutoff%s.0divg.4col.longNames RM4col 0" 
    # % (scriptDir, base, repCut, repDiv, base, base, repCut, repDiv)) #I need this just so I don't have to modify the next script (I know it's bad programming, but it was a time saver)


    print("\nGet simple repeats")
    os.system("grep -v Simple_repeat %s.4col.true.fa.mod.mod.out." % base +\
              "%sCutoff%s.0divg.4col.longNames "  % (r_cut, r_div) +\
              "| grep -v Low_complexity | grep -v Satellite > " +\
              "%s.4col.true.fa.mod.mod.out." % base +\
              "%sCutoff%s.0divg.4col.longNames.filt" % (r_cut, r_div))

    print("\nFilter out repeats")
    os.system(
        "python %s/1_3_filterOutPseuRMs.py " % scriptDir +\
        "%s.4col.true.fa.mod.mod.out." % base +\
        "%sCutoff%s.0divg.4col.longNames.filt " % (r_cut, r_div) +\
        "%s.fullyFiltered.disable_count " % base +\
        "%s.4col.true %s.4col.true.fa " % (base, base)+\
        "%s.fullyFiltered.disable_count.RMfilt " % base +\
        "%s.4col.true.RMfilt %s.4col.true.fa.RMfilt" % (base, base))
    print("  repeat Masked files end with '.RMfilt'")

    # Cleanup
    #os.system("mkdir Temp")
    #os.system("mv %s.tblastn_* Temp" % (base))
    #os.system("mv %s.4col.true.fa.mod.mod* Temp" % (base))
    #os.system("mv Temp/*count .")

    #########
    ## Step 8
    print("\n#### Step 8. Run high-confidence filter")
    # Do high confidence filter

    print("\nGet genome and protein sequence sizes")
    os.system("python %s/FastaManager.py -f get_sizes -fasta %s" % \
                                                        (scriptDir, fasta))
    os.system("python %s/FastaManager.py -f get_sizes -fasta %s" % \
                                                        (scriptDir, prot))

    print("\nRun High-Confidence Filter")
    os.system("python %s/1_3_v2_verifyPseudogenes.py " % scriptDir +\
                "%s.fullyFiltered.disable_count.RMfilt " % base +\
                "%s.4col.true.RMfilt %s.4col.true.fa.RMfilt " % (base, base) +\
                "%s.size %s.size " % (fasta, prot) +\
                "%s.fullyFiltered.disable_count.RMfilt.hiConf " % base +\
                "%s.4col.true.RMfilt.hiConf " % base +\
                "%s.4col.true.fa.RMfilt.hiConf %s" % (base, intronLen))

    print("  high-Confidence files: '.hiConf'")

    # Bring back to code names %s.4col.true.RMfilt %s.4col.true.fa.mod.mod.RMfilt
    print("\nReplace gene-names")
    os.system("python %s/TabularManager.py -f RemCodeName -inp %s.4col.true.RMfilt.hiConf -out %s.4col.true.RMfilt.hiConf.cdnm -col 1 -ref %s.4col.true.fa.ref -FR R" % (scriptDir, base, base, base))
    os.system("python %s/FastaManager.py -f replace_names -fasta %s.4col.true.fa.RMfilt.hiConf -name %s.4col.true.fa.ref -out %s.4col.true.fa.RMfilt.hiConf.cdnm" % (scriptDir, base, base, base))
    print("  name corrected files end with '.cdnm'")

    print("\nWrite Pseudogene GFF")
    os.system("python %s/1_3_getPseudoGFF.py %s.4col.true.RMfilt.hiConf.cdnm %s.fullyFiltered.disable_count.RMfilt.hiConf %s.4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.gff" % (scriptDir, base, base, base, base))
    print("  GFF files:", base + ".4col.true.fa.ref", 
                          base + ".4col.true.RMfilt.hiConf.cdnm.gff")

    # 05/12/22: SHS: The output file names below does not make sense - has a
    #  space and the string formating is specified but no string is provided.
    #  So I am commenting this out.
    if mismatch_status == "truncated":
        print("ERR: Problem with this part of the code that remain unresolved")
        """
        print("\nRewriting Pseudogene GFF with original names")
        gff_source = open(base + ".4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.gff",'r')
        gff_original_out = open(base + ".4col.true.fa.ref %s.4col.true.RMfilt.hiConf.cdnm.original.gff",'w')
        for line in gff_source:
            split_line = line.strip().split("\t")
            chromosome = split_line[0]
            seq_info = split_line[8]
            original_chromosome = original_contig[chromosome]
            split_seq_info = seq_info.strip().split(";")
            derivation = split_seq_info[1]
            split_derivation = derivation.strip().split("=")
            genename = split_derivation[1]
            original_genename = original_protein[genename]
            original_derivation = split_derivation[0] + "=" + original_genename
            original_seq_info = split_seq_info[0]+";"+original_derivation+";"+split_seq_info[2]
            original_line = original_chromosome + "\t" + "\t".join(split_line[1:8]) + "\t" + original_seq_info + "\n"
            gff_original_out.write(original_line)
        gff_source.close()
        gff_original_out.close()
        print("  GFF file with non-truncated names:", base + \
                                ".4col.true.RMfilt.hiConf.cdnm.original.gff")
        """

    #########
    ## Step 9
    print("\n#### Step 9. Clean up")

    print("\nRemoving temporary genome and protein sequences files")
    os.system("rm %s.*" % (P["p_seq"])) 
    os.system("rm %s.*" % (P["g_seq"]))
    os.system("rm temp")
    os.system("rm -rf RM_*") # repeat masker tmp folder

    dir_int = "_intermediate"
    dir_res = "_results"
    dir_log = "_logs"

    print("Moving intermediate files to _intermediate")
    if not os.path.isdir(dir_int): os.system(f"mkdir {dir_int}")
    os.system("mv %s* %s" % (P["b_out"], dir_int))
    os.system("mv test_gff.gene4col %s" % dir_int)

    print("Move hiConf, codename and gff files to _results")
    if not os.path.isdir(dir_res): os.system(f"mkdir {dir_res}")
    os.system(f"mv {dir_int}/*.hiConf {dir_res}")
    os.system(f"mv {dir_int}/*.hiConf.cdnm {dir_res}")
    os.system(f"mv {dir_int}/*.hiConf.cdnm.gff {dir_res}")

    print("Move log files to _logs")
    if not os.path.isdir(dir_log): os.system(f"mkdir {dir_log}")
    os.system(f"mv log* {dir_log}")

# If not loaded as a module
if __name__ == '__main__':
    main()

    print("\nDone!!!\n")
