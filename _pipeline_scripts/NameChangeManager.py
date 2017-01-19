#This script is designed to do all things related to changing names within a
#file.  Example: in col 1 AT3G34029.1 --> AT3G34029
#Created by David E. Hufnagel on July 9, 2012

import sys

#functions
class NameChanger:
    #standard name change in a fasta file
    def FastaNameChanger(self):
        for line in inp:
            if not line.startswith("#"):
                if line.startswith(">"):
                    oldName = line.strip(">")
                    if spe == "At":   #AT3G34029.1 --> AT3G34029
                        newName = self.AtStandardMod(oldName)
                    elif spe == "Al": #fgenesh2_kg.1__2__AT1G02190.2 --> fgenesh2_kg.1_2_AT1G02190.2
                        newName = self.AlStandardMod(oldName)
                    elif spe == "Rr": #RrC1_p3|1.000 --> RrC1_p3
                        newName = self.RrStandardMod(oldName)
                    else:
                        print "\n***unknown species: need more functionality***\n"
                    newLine = ">%s\n" % (newName)
                    out.write(newLine)
                else:
                    out.write(line)        

    #standard name change in a .disable_count file
    def DisableCountNameChanger(self, ind):
        for line in inp:
            if not line.startswith("#p"):
                if line.startswith("#"):
                    lineLst = line.split(" ")
                    oldName = lineLst[ind].strip("#").strip()
                    if spe == "At":   #AT3G34029.1 --> AT3G34029
                        newName = self.AtStandardMod(oldName)
                    elif spe == "Al": #fgenesh2_kg.1__2__AT1G02190.2 --> fgenesh2_kg.1_2_AT1G02190.2
                        newName = self.AlStandardMod(oldName)
                    elif spe == "Rr": #RrC1_p3|1.000 --> RrC1_p3
                        newName = self.RrStandardMod(oldName)
                    else:
                        print "\n***unknown species: need more functionality***\n"
                    lineLst.pop(ind)
                    lineLst.insert(ind, newName)
                    newLine = "#%s\n" % (" ".join(lineLst).strip())
                    out.write(newLine)
                else:
                    out.write(line)

    #concatenates all words in name into one _ delimited line
    def FastaSpaceKiller1(self):
        for line in inp:
            if not line.startswith("#"):
                if line.startswith(">"):
                    newName = "_".join(line.strip(">").split(" ")).strip()
                    newLine = ">%s\n" % (newName)
                    out.write(newLine)
                else:
                    out.write(line)

    #truncates everything in the name after the first space
    def FastaSpaceKiller2(self):
        for line in inp:
            if not line.startswith("#"):
                if line.startswith(">"):
                    newName = line.strip(">").split(" ")[0].strip()
                    newLine = ">%s\n" % (newName)
                    out.write(newLine)
                else:
                    out.write(line)

    #removes all empty lines from the file
    def KillSpaceLines(self):
        for line in inp:
            if not line.startswith("#"):
                if not line == "\n":
                    out.write(line)


    #replaces specific items in the fasta name with _s
    def FastaSingleItemKiller(self, item):
        #SPECIAL CASE:
        if item == "pipe":
            item = "|"
        
        #Do the actual processing
        for line in inp:
            if not line.startswith("#"):
                if line.startswith(">"):
                    oldName = line.strip(">")
                    newName = oldName.replace(item, "_").strip()
                    newLine = ">%s\n" % (newName)
                    out.write(newLine)
                else:
                    out.write(line)

    #replaces specific items in the fasta name with _s
    def TabSingleItemKiller(self, item):
        #SPECIAL CASE:
        if item == "pipe":
            item = "|"
        
        #Do the actual processing
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.strip("\n").split("\t")
                oldName = lineLst[ind]
                newName = oldName.replace(item, "_")
                lineLst.pop(ind)
                lineLst.insert(ind, newName)
                newLine = "\t".join(lineLst) + "\n"
                out.write(newLine)

    #standard name changes in a 4col file
    def FourColNameChanger(self, ind):
        for line in inp:
            if not line.startswith("#"):
                lineLst = line.split("\t")
                oldName = lineLst[ind]
                if spe == "At":   #AT3G34029.1 --> AT3G34029
                    newName = self.AtStandardMod(oldName)
                elif spe == "Al": #fgenesh2_kg.1__2__AT1G02190.2 --> fgenesh2_kg.1_2_AT1G02190.2
                    newName = self.AlStandardMod(oldName)
                elif spe == "Rr": #RrC1_p3|1.000 --> RrC1_p3
                    newName = self.RrStandardMod(oldName)
                else:
                    print "\n***unknown species: need more functionality***\n"
                lineLst.pop(ind)
                lineLst.insert(ind, newName)
                newLine = "\t".join(lineLst).strip() + "\n"
                out.write(newLine)

    def AtStandardMod(self, oldName):
        newName = oldName.split(".")[0].strip()
        return newName

    def AlStandardMod(self, oldName):
        newName = ""
        isLastUnderscore = False
        for char in oldName:
            if char == "_":
                if isLastUnderscore == False:
                    newName += char
                isLastUnderscore = True
            else:
                newName += char
                isLastUnderscore = False
        
        newName = newName.strip()
        return newName

    def RrStandardMod(self, oldName):
        newName = oldName.split("|")[0].strip()
        #if oldName.startswith("Rr"):
        #    newName = "_".join(oldName.split("_")[:-1]).strip()
        #else:
        #    newName = oldName
        return newName

    #Writes the users command line prompt on the first line of the output file.
    def WriteCallLine(self):
        out.write("#python %s\n" % (" ".join(sys.argv)))
        
    def Help(self):
        print "\nDescription:"
        print "This script is designed to do all things related to changing names within a"
        print "file.  Example: in col 1 AT3G34029.1 --> AT3G34029"
        print "Created by David E. Hufnagel on July 9, 2012"
        print "\nItems starting with '*' are capabilities not yet included"
        print "\nParameters:"
        print "    -f function:"
        print "        FastaNC - standard name change in a fasta file. REQ: inp, out,"
        print "            spe"
        print "        4colNC - standard name change in a 4col file. REQ: inp, out,"
        print "            spe, col"
        print "        DisNC - standard name change in a disable_count file. REQ: inp,"
        print "            out, spe, col"
        print "        FastaKillSpace1 - concatenates all words in name into one _"
        print "            delimited line. REQ: inp, out"
        print "        FastaKillSpace2 - truncates everything in the name after the"
        print "            first space.  REQ: inp, out"
        print "        KillSpaceLines - removes all empty lines from the file"
        print "            REQ: inp, out"
        print "        FastaKillSingleItem - replaces a single item in fasta names with"
        print "            _s REQ: inp, out, item (for item=| use item=pipe)"
        print "            out, item"
        print "        TabKillSingleItem - replaces a single item in fasta names with"
        print "            _s REQ: inp, out, col, item (for item=| use item=pipe)"
        print "    -inp - the input file to be modified."
        print "    -out - the modified output file."
        print "    -spe - two column species ID.  Available IDs:"
        print "         At = Arabidopsis thaliana"
        print "         Al = Arabidopsis lyrata"
        print "         *Br = Brassica rapa"
        print "         Rr = Raphinus raphanistrum"
        print "    -col - the column to be changed (number form, NOT AN INDEX)"
        print "    -item - for FastaKillSingleItem, it is the item to be replaced with _"
        print
        sys.exit(0)
        

 
NC = NameChanger()
inp=""; out=""; spe=""; col=0; ind = -1; item=""

#input handling
isEmpty = True #will be set to false if there are any parameters
for i in range(1,len(sys.argv),2):
    if sys.argv[i] == "-f":
        function = sys.argv[i+1]
    elif sys.argv[i] == "-inp":
        inp = open(sys.argv[i+1])
    elif sys.argv[i] == "-out":
        out = open(sys.argv[i+1], "w")
    elif sys.argv[i] == "-spe":
        spe = sys.argv[i+1]
    elif sys.argv[i] == "-col":
        col = int(sys.argv[i+1])
        ind = col-1
    elif sys.argv[i] == "-item":
        item = sys.argv[i+1]
    elif sys.argv[i] == "-h":
        NC.Help()
    else:
        print "UNKNOWN FLAG:",sys.argv[i]
        print "add -h to get help."
        sys.exit(0)
    isEmpty = False

#no parameters will bring up the help text
if isEmpty == True:
    NC.Help()

#similar name lists for calling functions
fastaNames = ["fasta", "Fasta", "FASTA"]
yesNames = ["Yes", "yes", "YES", "y", "Y"]
noNames = ["No", "NO", "no", "N", "n"]

#calling funcitons
NC.WriteCallLine()
if function == "FastaNC" or function == "fastaNC":
    if inp == "" or out == "" or spe == "":
        print "\n***missing parameters REQ: inp, out, spe***\n"
        NC.Help()
    NC.FastaNameChanger()

elif function == "4colNC" or function == "4ColNC":
    if inp == "" or out == "" or spe == "" or col == 0:
        print "\n***missing parameters REQ: inp, out, spe, col***\n"
        NC.Help()
    NC.FourColNameChanger(ind)

elif function == "DisNC" or function == "disNC":
    if inp == "" or out == "" or spe == "" or col == 0:
        print "\n***missing parameters REQ: inp, out, spe, col***\n"
        NC.Help()
    NC.DisableCountNameChanger(ind)

elif function == "FastaKillSpace1" or function == "fastaKillSpace1":
    if inp == "" or out == "":
        print "\n***missing parameters REQ: inp, out***\n"
        NC.Help()
    NC.FastaSpaceKiller1()

elif function == "FastaKillSpace2" or function == "fastaKillSpace2":
    if inp == "" or out == "":
        print "\n***missing parameters REQ: inp, out***\n"
        NC.Help()
    NC.FastaSpaceKiller2()
elif function == "KillSpaceLines" or function == "killSpaceLines":
    if inp == "" or out == "":
        print "\n***missing parameters REQ: inp, out***\n"
        NC.Help()
    NC.KillSpaceLines() 
elif function == "FastaKillSingleItem" or function == "fastaKillSingleItem":
    if inp == "" or out == "" or item == "":
        print "\n***missing parameters REQ: inp, out, item***\n"
        NC.Help()
    NC.FastaSingleItemKiller(item)

elif function == "TabKillSingleItem" or function == "tabKillSingleItem":
    if inp == "" or out == "" or item == "" or col ==0:
        print "\n***missing parameters REQ: inp, out, item***\n"
        NC.Help()
    NC.TabSingleItemKiller(item)    

else:
    print "\n***unknown function: need more functionality***\n"
