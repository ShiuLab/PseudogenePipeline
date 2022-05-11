
import sys

print "Read protein sequence file..."
inp = open(sys.argv[1])
inl = inp.readlines()
p   = {}
for i in inl:
	if i[0] == ">":
		g = i[1:-1].split(".")
		if g[0] not in p:
			p[g[0]] = [g[1]]
		else:
			p[g[0]].append(g[1])

print "Read pair file..."
inp = open(sys.argv[2])      # osv5_ps_gene.pairs
oup = open("osv5_ps_prot.pairs","w")
inl = inp.readlines()
miss = []
for i in inl:
	L = i[:-1].split("\t")
	if L[1] in p:
		for j in p[L[1]]:
			oup.write("%s\t%s.%s\n" % (L[0],L[1],j))
	else:
		if L[1] not in miss:
			miss.append(L[1])

print "The following genes are not in the prot seq file:"
for i in miss:
	print "",i
		
print "Done!"
