# For filtering BLAST output to focus on just 25 genes
#

blast = open('test_tblastn_tabular.out').readlines()
seq   = open('test_protein_sequences_25.fa').readlines()

seq_dict = {}
for i in seq:
    if i[0] == ">":
        seq_name = i.strip()[1:]
        seq_dict[seq_name] = 1

found = 0
with open('test_tblastn_tabular_25.out', 'w') as f:
    for match in blast:
        query = match.split('\t')[0]
        if query in seq_dict:
            f.write(match)
            found += 1

print("Found:", found)


