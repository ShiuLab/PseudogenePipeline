# IMPORTS
import sys

# FUNCTIONS
def GFFto4col(gff,feature,f_index,c_index,start,stop,name,truncate):
    # Creates a 4 column file from of selected features from a GFF file and return the name of the 4col file
    #
    # Inputs
    #    gff      := name of GFF file
    #    feature  := type of feature you are looking for
    #    f_index  := column position of features in GFF file
    #    c_index  := column position of chromosome in GFF file
    #    start    := column position of start site in GFF file
    #    stop     := position of stop site in GFF file
    #    name     := position of information from which the sequence name can be derived
    #    truncate := determines whether or not to truncate gene-names because of spaces
    #    NOTE: One is substracted from each index b/c python begin indexing at 0
    # Outputs
    #    returns "gene_4col" the name of the 4 column file which is written to the current folder

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

# MAIN
gff = sys.argv[1]
feature = sys.argv[2]
f_index = int(sys.argv[3])
c_index = int(sys.argv[4])
start = int(sys.argv[5])
stop = int(sys.argv[6])
name = int(sys.argv[7])
truncate = sys.argv[8]
GFFto4col(gff,feature,f_index,c_index,start,stop,name,truncate)
