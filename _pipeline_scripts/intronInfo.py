#####################################
# last update 02/26/2013 by J. Mass #
# version = 0.01                    #
#####################################
import imp,sys, getopt, re
GFFFastaTools = imp.load_source('GFFFastaTools', './GFFFastaTools.py')
from GFFFastaTools import GFFParser,vprint
VERB=False
def usage():
    print ("""
    #######################################################
    # intronInfo v0.1                                     #
    # in: GFF3 file                                       #
    # out:4 col table: GENE0|INTRON0\tCHR\tSTART\tEND
    #                  ##    ##                           #
    #                  GENEn|INTRONm\tCHR\START\tEND
    #######################################################
    -g, --gff=gff3_FILE      (tested w/ .gff3 from phytozome)
    -o, --output=FILE        output file [default: gff3_FILE.intron]
    -h, --help               prints this
    """)
    sys.exit(2)
##########################################################
class Gene(object):
    def __init__(self,id):
        self._id=id
        self._locus=None
        self._geneStart=None
        self._geneStop=None
        self._chromosome=None
        self._exons = []
        self._exonStart = {}
        self._exonStop = {}
        self._introns = []
        self._intronStart = {}
        self._intronStop= {}
    def addExon(self,exonId, exonStart, exonStop):
        self._exons.append(exonId)
        self._exonStart[exonId]=exonStart
        self._exonStop[exonId]=exonStop
    def calcIntron(self):
        
        if self._geneStart == self._exonStart[self._exons[len(self._exons)-1]]:
            tmp=self._exons
            tmp.reverse()
            self._exons=tmp        
        for i in range(0,len(self._exons)-1):
            vprint(VERB, "Exon START",i,self._exons[i],self._exonStart[self._exons[i]])
            vprint(VERB,"Exon STOP",i,self._exons[i], self._exonStop[self._exons[i]])
            vprint(VERB,"Intron from", i,self._exons[i],self._exonStop[self._exons[i]] )
            vprint(VERB,"Intron to", i+1,self._exons[i],self._exonStart[self._exons[i+1]] )
            self._introns.append("intron."+str(i))
            self._intronStart["intron."+str(i)]=self._exonStop[self._exons[i]]+1
            self._intronStop["intron."+str(i)]=self._exonStart[self._exons[i+1]]-1
            #todo
def getIntronInfo(gff):
    res = None
    genes = {}
    longest={}
    lociSeen = set()
    for el in GFFParser().parse(gff):
        if (el.type=="mRNA"):
            #take only longest splice variant into account
            if (el.attrib_dct["longest"]==str(1)):
                longest[el.attrib_dct["ID"]]=el.attrib_dct["Parent"]
                lociSeen.add(el.attrib_dct["Parent"])
                if el.attrib_dct["ID"] not in genes:
                    #vprint(VERB,el.attrib_dct["Parent"]+"|"+el.attrib_dct["ID"])
                    genes[el.attrib_dct["ID"]]=Gene(id=el.attrib_dct["ID"])
                    genes[el.attrib_dct["ID"]]._chromosome=el.seqid
                    genes[el.attrib_dct["ID"]]._geneStart=el.start
                    genes[el.attrib_dct["ID"]]._geneStop=el.end
                    genes[el.attrib_dct["ID"]]._locus=el.attrib_dct["Parent"]
                elif el.attrib_dct["ID"] in genes:
                    print(el.attrib_dct["ID"],"###\tdouble entry?\n")
    seen =set()
    for el in GFFParser().parse(gff):
        if (el.type=="exon" and el.attrib_dct["Parent"] in longest):
            genes[el.attrib_dct["Parent"]].addExon(el.attrib_dct["ID"], el.start,el.end)
    #check whether all loci are covered
    for el in GFFParser().parse(gff):
        if (el.type=="gene"):
            if (el.attrib_dct["ID"] not in lociSeen):
                print(el.attrib_dct["ID"],"No longest splice form marked in gff\n")
    
    gk=sorted(genes.keys())
    for g in gk:
        genes[g].calcIntron()
        #name= genes[g]._locus+"|"+g
        #chrom=genes[g]._chromosome
        #start=genes[g]._geneStart
        #stop=genes[g]._geneStop
         
    return(genes)
##############################################################

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "g:o:he", ["gff=", "output=","help","incl_exons"])         
except getopt.GetoptError as err:
    print (str(err))
    usage()

gff=None
outfile=None
includeExons=False
for o, a in opts:
    if o in ("-g", "--gff"):
        gff = a
    elif o in ("-h", "--help"):
        usage()
    elif o in ("-e", "--incl_exons"):
        includeExons=True
    elif o in ("-o", "--output"):
        outfile=a
    else:
        assert False, "unhandled option"

if not gff:
    usage()
else:
    if not outfile:
        outfile=gff+".introns"
    outfile = open(outfile,"w")
    res = getIntronInfo(gff)
    genes = sorted(res.keys())
    noIntrons = []
    for g in genes:
        if (len(res[g]._introns) ==0):
            noIntrons.append(g)
        else:
            tmp="## "+res[g]._locus+" "+res[g]._chromosome+" "+str(res[g]._geneStart)+" "+str(res[g]._geneStop)+" ##\n"
            for i in res[g]._introns:
                name= res[g]._locus+"|"+str(i)
                chrom=res[g]._chromosome
                start=res[g]._intronStart[i]
                stop=res[g]._intronStop[i]
                tmp += "\t".join([name,chrom, str(start), str(stop)])
                tmp +="\n"
            if includeExons==True:
                for i in res[g]._exons:
                    name= res[g]._locus+"|"+str(i)
                    chrom=res[g]._chromosome
                    start=res[g]._exonStart[i]
                    stop=res[g]._exonStop[i]
                    #print(g,genes[g]._geneStart, genes[g]._geneStop)
                    tmp += "\t".join([name,chrom, str(start), str(stop)])
                    tmp +="\n"
            outfile.write(tmp)
            vprint(VERB,tmp)
    for n in noIntrons:
        tmp="## [noIntrons] "+res[n]._locus+" "+res[g]._chromosome+" "+str(res[g]._geneStart)+" "+str(res[g]._geneStop)+" ##\n"
        outfile.write(tmp)
    outfile.close()
############################################







