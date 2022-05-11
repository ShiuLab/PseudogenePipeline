#This script is designed to take GOs online file, gene_ontology_ext.obo, and
#make a tab delimited file from it.
#Created by David E. Hufnagel on June 27, 2012

import sys

inp = open(sys.argv[1])      #input GO file
out = open(sys.argv[2], "w") #output tabular GO file





class GOcat:
    def __init__(self, ID, name, namespace, defin, theRest):\
        #ID
        if ID != "":
            self.ID = ID
        else:
            print "\n***empty ID***\n"
            
        #name
        if name != "":
            self.name = name
        else:
            print "\n***empty name***\n"

        #namespace
        if namespace != "":
            self.namespace = namespace
        else:
            print "\n***empty ID***\n"

        #defin
        if defin != "":
            self.defin = defin
        else:
            print "\n***empty ID***\n"

        #theRest
        if type(theRest) == list:
            self.theRest = theRest
        else:
            print "\n***type error 1***\n"

    def __str__(self):
        return "id: %s    name: %s    def: %s    theRest: %s\n" % \
               (self.ID, self.name, self.defin, ";".join(self.theRest))

##    def AddToTheRest(self, more):
##        if type(more) == list:
##            for item in more:
##                self.theRest.append(item)
##        else:
##            print "\n***type error 2***\n"





#writes the users command line prompt on the first line of the output file.
out.write("#python %s\n" % (" ".join(sys.argv)))

active = False  #Whether to take in information for a GO term
started = False #when to start making and keeping GOcats (on the second line)
GOlist = []
for line in inp:
    if line.startswith("[Term]"):
        if started == True:
            GO = GOcat(ID, name, namespace, defin, theRest)
            GOlist.append(GO)
        started = True
        active = True

        ID = ""
        name = ""
        namespace = ""
        defin = ""
        theRest = []
    elif line.startswith("\n"):
        active = False
    else:
        if not line.startswith("[Typedef]"):
            lineLst = line.split(":")
            if active == True:
                if line.startswith("id:"):
                    ID = ":".join(lineLst[1:]).strip()
                elif line.startswith("name:"):
                    name = lineLst[1].strip().replace(" ", "_")
                elif line.startswith("namespace:"):
                    namespace = lineLst[1].strip().replace(" ", "_")
                elif line.startswith("def:"):
                    defin = ":".join(lineLst[1:]).strip().replace(" ", "_")
                else:
                    theRest.append(line.replace(" ", "_").strip())
        
for GO in GOlist:
    newLine = "%s\t%s\t%s\t%s\t%s\n" % (GO.ID, GO.name, GO.namespace, GO.defin, ";".join(GO.theRest))
    out.write(newLine)        








inp.close()
out.close()
