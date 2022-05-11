import re, string
#####################################
# last update 02/26/2013 by J. Mass #
# version = 0.04                    #
#####################################


#""" print verbosely via vprint if 'VERBOSE' is true"""
def vprint(verbose,*args):
    """ print verbosely VERBOSE is True """
    if verbose:
        print(*args)
        pass
    else:
        pass
#------------------------------------


class FastaHelper(object):
    def remove_newlines(self, string):
        yield string.replace("\n","")
    
    def insert_newlines(self, string, every=80):
        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
    
    def complement(self, string):
        tr = str.maketrans('AGTCagtc','TCAGtcag')
        res = string.translate(tr)
        return(res)
    def reverseComplement(self, string):
        tr = str.maketrans('AGTCagtc','TCAGtcag')
        res = string.translate(tr)
        return res[::-1]

class FastaParser(object):     
    def read_fasta(self, fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta' 
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"' returns ("idpart1", "ATGTGA")
        """
        name = ""
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            if delim:
                vprint("delim set to "+delim)
                line = line.split(delim)[asID]
            yield (line.rstrip(), joinedSeq.rstrip())    
        fasta.close()
        
            
        
class SeqTranslator(object):
    RNAmap = {
           "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
           "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
           "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
           "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
           "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
           "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
           "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
           "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
           "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
           "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
           "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
           "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
           "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
           "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }
    DNAmap = {
           "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
           "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
           "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
           "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
           "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
           "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
           "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
           "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
           "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
           "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
           "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
           "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
           "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
           "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
           "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
           "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }
    def triplets(self, string, frameshift = 0):
        #chop string into 3 char slices 
        for i in range(0, int(len(string)/3)):
            #print(string[i*3+frameshift:i*3+3+frameshift])
            yield string[i*3+frameshift:i*3+3+frameshift]
    def triplets1(self, string, frameshift = 0):
        if frameshift == 3:
            frameshift=0
        #chop string into 3 char slices 
        for i in range(0, int(len(string)/3)):
            #print(string[i*3+frameshift:i*3+3+frameshift])
            if frameshift==0:
                yield string[i*3+frameshift:i*3+3+frameshift]
            elif frameshift==1:
                yield string[1+ i*3+frameshift:1+ i*3+3+frameshift]
            elif frameshift ==2:
                yield string[2+ i*3+frameshift:2+ i*3+3+frameshift]
            else:
                print("\t??")
    def dna2prot(self, string, frameshift=0):
        res = ""
        for a in self.triplets(string, frameshift):
            try:
                res+=self.DNAmap[a]
            except KeyError as e:
                print(e, "Don't panic.")
                res+=""
        return res


class LargeFastaParser(object):
    def prepare_fasta(self, fasta):
        fasta = open(fasta, "r")
        fasta.seek(0)
        headerIndex = {}
        headerIndexNum = []
        l = fasta.readline()
        ind = -1
        while l:
            if l.startswith(">"):
                ind +=1
                pos = fasta.tell()
                l = l.split(" ")
                l = l[0]
                headerIndex[l.rstrip()[1:]] = pos
                headerIndexNum.append(pos)
            l = fasta.readline()
        pos = fasta.tell()
        ind+=1
        headerIndex["EOF"] = pos
        headerIndexNum.append(pos)
        fasta.close()
        return headerIndex
    
    def cacheSequence(self, fasta, offset):
        fasta = open(fasta, "r")
        fasta.seek(offset)
        tmp = fasta.readline()
        res = tmp
        while (not tmp.startswith(">")):
            tmp = fasta.readline()
            res = res+tmp
            if not tmp:
                break
        res = res.strip()
        res = res.replace("\n","")
        fasta.close()
        return res
    
    def getSequenceByCoordinatesOnLandmark(self, landmarkSeq, start, end, strand):
        start = int(start)-1
        end = int(end)
        if (strand == "+" or strand == "."):
            return(landmarkSeq[start:end])
        elif strand == ("-"):
            res = landmarkSeq[start:end]
            tr = str.maketrans('AGTCagtc','TCAGtcag')
            res = res.translate(tr)
            return res[::-1]
        else:
            print(strand)
            print("WARNING: strand was of unknown orientation\n Will return as if it was '+'\n")
            return landmarkSeq[start:end]
        
        return landmarkSeq[start-1:end-1]

class GFFParser():
    def __init__(self):
        _ids_name = {}
        _mrna_exons = {}
        _exons_mrna = {}

    @property
    def ids_name(self):
        return self._ids_name
    
    @property
    def mrna_exons(self):
        return self._mrna_exons
    
    def parse(self, gff):
        gff = open(gff,'r')
        for l in gff:
            if l.startswith("#"):
                pass
            else:
                l = l.rstrip()
                l = l.split("\t")
                tmp_seqid, tmp_source, tmp_type, tmp_start, tmp_end, tmp_score, tmp_strand, tmp_phase, tmp_attributes = l
                gffobj = gffObject(seqid=tmp_seqid, source=tmp_source, type=tmp_type, start=int(tmp_start), end=int(tmp_end), score=tmp_score, strand = tmp_strand, phase = tmp_phase, attributes =  tmp_attributes)
                yield(gffobj)
        gff.close()
    def parentsDct(self, gff):
        dct = {}
        roots = set()
        mrna = {}
        for a in self.parse(gff):
            if not "Parent" in a.attrib_dct:
                roots.add(a.attrib_dct["ID"])
                dct[a.attrib_dct["ID"]] = a.attrib_dct["ID"]
            else:
                if a.type == "mRNA":
                    mrna[a.attrib_dct["ID"]]= a.attrib_dct["Name"]
                    mrna[a.attrib_dct["Name"]]= a.attrib_dct["Name"]
                    #print(a.attrib_dct["ID"], mrna[a.attrib_dct["ID"]])
                dct[a.attrib_dct["ID"]] = a.attrib_dct["Parent"]#todo mult
        #collapse
        for k in dct:
            while (dct[k] not in roots):
                if dct[k] in mrna:
                    dct[k] = mrna[dct[k]]
                    break
                dct[k] = dct[dct[k]]
        return dct, mrna

class gffObject(object):
    """
    "seqid" (gff column 1) landmark for coordinate system
    "source" (gff column 2) source db/program etc
    "type" (gff column 3) term from the Sequence Ontology
    "start"(gff column 4) relative to the landmark seqid
    "end" (gff column 5) 1-based integer coordinates
    "score" (gff column 6) float
    "strand" (gff column 7) +/i/./? for positive/negative/none/unknown strand (relative to the landmark)
    "phase" (gff column 8) 0/1/2 for "CDS" features, start rel. to reading frame 
    "attributes" (gff column 9)  list of feature attributes 
    in the format tag=value.
    Multiple tag=value pairs are separated by semicolons.
    ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, 
    Dbxref, Ontology_term, Is_circular
    Parent: groups exons into transcripts, transcripts into genes etc.
        A feature may have multiple parents.  
    Target: Indicates the target of a nucleotide-to-nucleotide 
        or protein-to-nucleotide alignment. 
        The format of the value is "target_id start end [strand]", 
        where strand is optional and may be "+" or "-". 
    Gap: The alignment of the feature to the target if the two 
        are not collinear (e.g. contain gaps). 
        The alignment format is taken from the CIGAR format described
        in the Exonerate documentation. 
        (http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl). ("THE GAP ATTRIBUTE") 
    Parent, the Alias, Note, DBxref and Ontology_term attributes can have multiple values.
    """
    def __init__(self, seqid, source, type, start, end, score, strand, phase, attributes):
        self._seqid = seqid
        self._source = source
        self._type = type
        self._start = start
        self._end = end
        self._score = score
        self._strand = strand
        self._phase = phase
        self._attributes = attributes
        self._attrib_dct = {}
        
        pattern  = re.compile(r"""(\w+)=([\w:.,]+)""")
        tmp = pattern.findall(self._attributes)
        for a in tmp:
            key, val = (a)
            #todo ...
            self._attrib_dct[key]=val
    @property
    def seqid(self):
        return self._seqid
    @property
    def source(self):
        return self._source
    @property
    def type(self):
        return self._type
    @property
    def start(self):
        return self._start
    @property
    def end(self):
        return self._end
    @property
    def score(self):
        return self._score
    @property
    def strand(self):
        return self._strand
    @property
    def phase(self):
        return self._phase
    @property
    def attributes(self):
        return self._attributes
    @property
    def attrib_dct(self):
        return self._attrib_dct
    #todo...
    def toGFF3line(self):
        pass
    def getParents(self):
        pass
    
#---------------------------------------------------------#
class PhylipCompatibilizer: 
    @property
    def baseAlphabet(self):
        return self._baseAlphabet
    @property
    def knownIDs(self):
        return self._knownIDs
    @baseAlphabet.setter
    def baseAlphabet(self, value):
        self._knownIDs = value
    #Base62  should be plenty for shortening identifiers
    def __init__(self, baseAlphabet = string.digits + string.ascii_letters, seed = 0):
        #todo: seed from number
        self._baseAlphabet = baseAlphabet
        self._baseDict = dict((c, i) for i, c in enumerate(self._baseAlphabet))
        self._knownIDs = set()
        
    def decode(self,string):
        revBase = self._baseDict
        length = len(revBase)
        res = 0
        for i, c in enumerate(string[::-1]):
            print (i, c)
            res += (length ** i) * revBase[c]
        return res

    def encode(self,integer):
        base = self._baseAlphabet
        length = len(base)
        res = ''
        if integer in self._knownIDs:
            print("not globally unique!")
            integer = self._knownIDs+1
            self._knownIDs.add(integer)
        while integer != 0:
            res = base[integer % length] + res
            integer //= length
        if self._globalTF==True and res in self._knownIDs:
            print("not globally unique!")
            res= self.encode((self._knownIDs)+1)
        if len(res)>6:
            print("!!! Length > 6 !!!")
        
        return res