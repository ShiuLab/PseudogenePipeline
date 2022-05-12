#!/usr/local/bin/python

##
#
# 08/21,02
#  Implement compare_lists()
# 09/23,02
#  Problem with DOS style files. In get_sequences, "\r\n" is explicitly detected
#  but should use a generate function to convert file first.
# 01/02,03
#  The function get_group_seq don't output sequences in the order of the input
#  list. Make it so.
# 05/27,03
#  This is very significant so I am writing this down. A new code seqment for
#  fasta_to_dict is written, now it proceed at BLAZING SPEED.
# 08/24/2012 and beyond
#  Added some functions (David E. Hufnagel)
##

import os, sys, FileUtility, string, Translation

class fasta_manager:

    def __init__(self):
        pass

    #
    # @param coords  [seq_id][L][R]
    #
    def mask(self,fasta,coords):
        
        print("Read fasta into dict...")
        F = self.fasta_to_dict(fasta,dflag=0)
        
        print("Read coords into dict...")
        C   = {}
        inp = open(coords)
        inl = inp.readline()
        while inl != "":
            L = inl.split("\t")
            if L[0] in C:
                C[L[0]].append([int(L[1]),int(L[2])])
            else:
                C[L[0]] = [[int(L[1]),int(L[2])]]
            inl = inp.readline()
        
        print("Concatenate seq...")
        for i in C:
            if i in F:
                print("",i,"%i features" % len(C[i]))
                c = 0
                for j in C[i]:
                    if c % 100 == 0:
                        print(" %i x100" % (c/100))
                    c += 1
                    F[i] = "%s%s%s" % \
                            (F[i][:j[0]-1],"N"*(j[1]-j[0]+1),F[i][j[1]:])
            else:
                print("SEQ ABSENT:",i)
        
        print("Write sequences...")
        self.dict_to_fasta(F,fasta+".mask")
        
        
    def dict_to_fasta(self,fdict,fname):
        
        oup = open(fname,"w")
        for i in fdict:
            oup.write(">%s\n" % i)
            s = fdict[i]
            c = 0
            while c < len(s):
                oup.write("%s\n" % s[c:c+80])
                c += 80
        oup.close()

    ##
    # If multiple coords are passed, only the last one of the coords for each
    # seqeunce will be stored. Implement multiple get later...
    #
    # @param coords  this can be a string with coords separated by "," or a file
    #                with coordinates specified.
    # @param format  0 (default). 1 to get nt sequences based on peptide coords.
    #                The sequence file passed HAS TO BE nt and synchronized with
    #                the peptide sequences where coords are derived from.
    # @param verbose for fasta_to_dict
    # @param wc      ouput sequence name with coordinate or not
    # @param call    non-command line call [1]. Return a dict instead of
    #                file output.
    # @param add     add domain name to the sequence name. The coords passed
    #                have to be [gene][L][R][dom]
    ##
    def get_stretch(self,fasta,coords,format=0,wc=0,verbose=0,call=0,add=0):
        
        print("Fasta  :",fasta)
        print("Coord  :",coords)
        print("Format :",format)
        print("W/coord:",wc)
        print("verbose:",verbose)
        print("Call   :",call)
        print("Add_dom:",add)
        
        if format:
            print("CAUTION: you try to get nt sequences using pep coords")
        
        print("Read sequence into dict")
        fdict = self.fasta_to_dict(fasta,0,verbose)
        fkeys = list(fdict.keys())
        print("%i sequences..." % len(fkeys))
        
        #print fdict    
        cdict = {}
        only1 = 0  # all sequences have the same coords
        
        print("Read coords...")
        # internal call, dict will be passed
        if call == 1:
            cdict = coords
        else:
            # check if coords is a file
            try:
                inp = open(coords,"r")
                inl = inp.readline()
                while inl != "":
                    L = self.rmlb(inl).split("\t")
                    if L[0] in cdict:
                        cdict[L[0]].append(L[1:])                   
                    else:
                        cdict[L[0]] = [L[1:]]
                    inl = inp.readline()
            except IOError:
                print("Only 1 pair of coords...")
                only1 = 1
                coords = coords.split(",")
                for i in fdict:
                    cdict[i] = 0
                
        # for each sequence
        oup = open(fasta+".seg.fa","w")
        count_g = 0
        count_t = 0
        count_a = 0
        ckeys = list(cdict.keys())
        ckeys.sort()
        for i in ckeys:
            #print i
            if i not in fdict:
                print("E1: seq absent,",i)
                count_a += 1
                continue
            seq = fdict[i]
            slen = len(seq)
            out3prime = 0
            
            if only1 and len(coords) != 2:
                # take multiple stretches
                c = 0
                seg = ""
                while c < len(coords):
                    try:
                        if int(coords[c]) < 1 or int(coords[c+1]) < 1:
                            print(" coords smaller than 1, quit!\n")
                            sys.exit(0)
                    except ValueError:
                        print(coords)
                        sys.exit(0)
                    
                    coords[c]   = int(coords[c])
                    coords[c+1] = int(coords[c+1])
                    if format:
                        coords[c]   = coords[c]*3-2
                        coords[c+1] = coords[c+1]*3

                    seg = seg + seq[coords[c]-1:coords[c+1]]                
                    if coords[c+1] > slen:
                        out3prime = 1           
                    c += 2
                    
                if wc:
                    oup.write(">%s_%s\n%s\n"%(i,string.joinfields(coords),seg))
                else:
                    oup.write(">%s\n%s\n" % (i,seg))
                count_t += 1
            else:
                if only1:
                    clist = [coords]
                else:
                    clist = cdict[i]
                count_c = 1
                
                for j in clist:
                    seg = ""
                    
                    # check the first pair, if it is in reverse ori and there
                    # are more than one pair of coord, will stop right here.
                    if len(j) > 2 and int(j[0]) > int(j[1]):
                        print("Coord format problem... Quit!")
                        sys.exit(0)
                    
                    # iterate coordinate sets, some may be multiple
                    for k in range(0,len(j)-1,2):
                        cL = int(j[k])
                        cR = int(j[k+1])
                        if int(j[k]) > int(j[k+1]):
                            cR = int(j[k])
                            cL = int(j[k+1])
                        if format:
                            cL = cL*3-2
                            cR = cR*3                       
                        seg += seq[cL-1:cR]
                                            
                    # add domain name, dom should be last token
                    if add:
                        idx = "%s|%s" % (i,j[-1]) 
                    elif len(clist) == 1:
                        idx = i
                    else:
                        idx = "%s_seg%i" % (i,count_c)
                    if wc:
                        idx += "_%s" % (string.joinfields(j,"|"))

                    # now, the coord may be in reverse. only check the first
                    # pair to see if this is the case. if so, the reverse
                    # complement is generated. Obviously won't apply to peptide
                    # sequences
                    if int(j[0]) > int(j[1]):
                        print(j)
                        seg = trans.rc(seg)
                    
                    oup.write(">%s\n%s\n" % (idx,seg))
                    count_c += 1
                    
                    if int(j[1]) > slen:
                        out3prime = 1
                    count_t += 1
            
            count_g += 1
            if out3prime:
                pass
                #print "%s:out of 3' end, take as much as it can" % i
        
        print("Total %i stretches in %i genes, %i id not in seq file" % \
                (count_t,count_g,count_a))
                          
        
    ##
    # This is written to get stretch from a LARGE sequence file with a single
    # sequence, such as the chromosome sequence. The other method takes forever
    # because the fasta is concatenated into a single string. A SERIOUS negative
    # for this function is that the coords CANNOT overlap.
    #
    # @param fasta  the implementation only work for one id one seq in the file
    # @param coords separated by ",". Sequence obtained is inclusive of the
    #               coords OR a file name can be passed. The file contains
    #               a pair of coords separated by "," in each line
    # @param call   whether return the sequence as dict [1] or written [0]
    # @param isfile coords in a file
    ##
    def get_stretch2(self,fasta,coords,call=0,isfile=0):
        
        print("Fasta   :",fasta)
        print("Coord   :",coords)
        print("isfile  :",isfile)

        # check if coord is a file
        if isfile:
            inp = open(coords,"r")
            inl = inp.readline()
            
            # first read to a dict to sort the coords
            order = {}
            while inl != "":
                L = self.rmlb(inl).split(",")
                # assume coord are not overlapping
                order[int(L[0])] = int(L[1])
                inl = inp.readline()
            
            okeys = list(order.keys())
            okeys.sort()
            clist = []
            for i in okeys:
                clist.extend([i,order[i]])
            inp.close()
        else:
            clist = coords.split(",")
            clist = [int(clist[0]),int(clist[1])]
        
        inp = open(fasta,"r")
        idx = self.rmlb(inp.readline())[1:]
        
        # rid of anything after the 1st space
        if idx.find(" ") != -1:
            idx = idx[:idx.find(" ")]
        
        inl = inp.readline()
        c     = 0
        get   = 0
        odict = {}
        ostr  = ""
        end   = 0
        read  = 1
        
        saved = 0
        while inl != "":
            inl = self.rmlb(inl)
            #print [inl]
            
            # only increment if this line has been read
            if read:
                c += len(inl)
            else:
                read = 1
            
            if int(clist[0]) < int(clist[1]):
                cL = int(clist[0])
                cR = int(clist[1])
            else:
                cR = int(clist[0])
                cL = int(clist[1])
            
            if c > cL:              
                if get == 0:
                    #print "%s-%s" % (clist[0],clist[1]),cL,cR
                    #print " start:",
                    if c < cR:
                        ostr += "%s\n" % inl[cL-c+len(inl)-1:]
                        get = 1
                    else:
                        #print "end1"
                        ostr += "%s\n" % inl[cL-c+len(inl)-1:\
                                             cR-c+len(inl)]
                        end = 1
                else:
                    if c < cR:
                        #print ">%s\n" % inl
                        #print ".",
                        ostr += "%s\n" % inl
                    else:
                        #print "end2"
                        ostr += "%s\n" % inl[:cR-c+len(inl)]
                        end = 1
                                        
                # store seq into dict and reset everything
                if end:
                    #print " reset"
                    #print " store->",[ostr]
                    
                    if saved % 10000 == 0:
                        print(" %i x 10k" % (saved/10000))
                    saved += 1
                    
                    # deal with coords in reverse
                    if int(clist[0]) > int(clist[1]):
                        ostr = trans.rc(ostr)
                        print("rc:",[ostr])
                    
                    odict["%s|%s" % (clist[0],clist[1])] = ostr                 
                    # reset everything
                    end    = 0
                    get    = 0
                    ostr   = ""
                    clist = clist[2:]
                    if clist == []:
                        break
                    # just in case that the next entry is also in this line.
                    if c > cL:
                        #print ">>>no_read"
                        read = 0
            # read line or not                      
            if read:        
                #print " new_line"
                inl = inp.readline()        
        
        if call:
            return odict
        else:
            if isfile:
                oup = open(coords+".seg.fa","w")
            else:
                oup = open(fasta +".seg.fa","w")
            
            for i in odict:
                oup.write(">%s_%s\n%s" % (idx,i,odict[i]))
                
                  

    # for multiple chr, derived from get_stretch3.
    # 2/24/12, coordinates used to be [chr][L][R] or can be multiple coordinates
    #          separated by ",". now allow only [seq][L1,R1,...]
    #
    # @param fasta  can have multiple sequences in one file
    # @param coords [seq][L1,R1,L2,R2...]. If L > R, take reverse complement.
    # @param seqid  Use 4th column of the coords file as sequence ID if available.
    #               (1) or not (0, default)
    def get_stretch4(self,fasta,coords,seqid):
        print("Sequence to dict...")
        seq = self.fasta_to_dict(fasta,0)
        # fasta_to_dict got rid of "\n" already
    
        #print seq.keys()
    
        c = 0 # count total
        m = 0 # count not in fasta
        if coords.find(",") == -1:
            print("Read coordinates...")
            inp = open(coords)
            oup = open(coords+".fa","w")
            oup2= open(coords+".missing","w")
            inl = inp.readline()
            while inl != "":    # Go through each coord
                if c % 1000 == 0:
                    print(" %i k" % (c/1000))
                c  += 1
                L  = inl.strip().split("\t")
                seqName = L[0]  # Sequnece name
                if seqName in seq:
                    if len(L) >= 2:
                        # L = [name, L, R], some may have 4th col which is IDs to be given.
                        if L[1].find(",") == -1:
                            # Deal with reverse ori
                            ori= 1; cL = int(L[1]); cR = int(L[2])
                            if cL > cR:
                                ori = -1
                            # Get sequence
                            if ori == -1:
                                S = seq[seqName][cR-1:cL]
                            else:
                                S = seq[seqName][cL-1:cR]
    
                            if S == "":
                                print("ERR COORD: %s,[%i,%i]" % (seqName,cL,cR))
                            else:
                                if ori == -1:
                                    S = trans.rc(S)
                                
                                # If there is 4th column, use them as sequence IDs.
                                if len(L) == 4 and seqid:
                                    oup.write(">%s\n%s\n" % (L[3],S))
                                else:
                                    oup.write(">%s|%i-%i\n%s\n" %   (seqName,cL,cR,S))
                        # name <\t> "L1,R1,L2,R2..." <\t> whatever
                        elif L[1].find(",") != -1:
                            coordList = L[1].split(",")
                            S = ""
                            # Set orientation, only consider the first pair
                            if int(coordList[0]) < int(coordList[1]):
                                ori = 1
                            else:
                                ori = -1
                            for j in range(0,len(coordList),2):
                                cL = int(coordList[j]); cR = int(coordList[j+1])
                                if ori == -1:
                                    cL = int(coordList[j+1])
                                    cR = int(coordList[j])
                                    S = seq[seqName][cL-1:cR] + S
                                else:
                                    S += seq[seqName][cL-1:cR]
                            oup.write(">%s|%s\n%s\n" % (seqName,"-".join(coordList),S))
                        else:
                            print("Unknown cooord format:",L)
                            print("Quit!")
                            sys.exit(0)
                else:
                    m += 1
                    oup2.write(inl)
                inl = inp.readline()
            print("Total coords:",c)
            print("Not in seq  :",m)
            oup2.close()

        # coordinates are passed
        else:
            coords = coords.split(",")
            print("Coords:",coords)
            oup = open("%s_%s.fa" % (fasta,"-".join(coords)),"w")
            C = []
            for i in coords:
                C.append(int(i))
            for i in seq:
                s = ""
                for j in range(0,len(C),2):
                    s += seq[i][C[j]-1:C[j+1]]
                oup.write(">%s\n%s\n" % (i,s))
                    
        oup.close()
              

    # From Gaurav
    # for multiple chr, derived from get_stretch3
    #
    # @param fasta  can have multiple sequences in one file
    # @param coords [chr][L][R], L can be bigger than R, take RC in this case
    def get_stretch5(self,fasta,coords):
        
        print("Sequence to dict...")
        seq = self.fasta_to_dict(fasta,0)
        # fasta_to_dict got rid of "\n" already
        
        print("Read coordinates...")
        inp = open(coords)
        oup = open(coords+".fa","w")
        #oup2 = open(coords+".fa.table","w")
        inl = inp.readline()
        c = 0
        while inl != "":
            if inl.startswith('#'):
                    oup.write(inl)
            else:
                if c % 1000 == 0:
                    print(" %i k" % (c/1000))
                c  += 1
                L  = self.rmlb(inl).split("\t")
                #Changed format is
                #Pname/Any ID -- Cname --Start--End
                #~~~~Gaurav
                
                pname=L[0]
                ch = L[1]
                left=L[2]
                right=L[3]
                '''
                ch=L[3]
                left=L[4]
                right=L[5]
                '''
                S  = ""
                #print L
                # deal with reverse ori
                ori = 1
                if int(L[2]) > int(L[3]):
                        ori = -1
                        #ori2='-1'
                #else:
                        #ori2 ='+1'
                
                #print [seq[ch]],len(seq[ch])
                for j in range(2,len(L),2):             
                    cL = int(L[j])
                    cR = int(L[j+1])
                    #print cL, cR
                    #sys.exit()
                    
                    if ori == -1:
                        cL = int(L[j+1])
                        cR = int(L[j])
                        S = seq[ch][cL-1:cR] + S
                    else:
                        S += seq[ch][cL-1:cR]
                    #print cL,cR,S
                if S == "":
                    print("ERR COORD: [%i,%i,%s]" % (cL,cR,pname))
                else:
                    if ori == -1:
                        S = trans.rc(S)
                    #oup.write(">%s|%s|%s\n%s\n" %(pname,ch,string.joinfields(L[2:],"-"),S))
                    oup.write(">%s\n%s\n" %(pname,S))
                    #oup2.write('%s\t%s\n'%\
                    #      ('\t'.join(L[0:]),S))
                    #oup2.write('%s\t%s\t%s\t%s\t%s\t%s\n'%\
                    #      (L[0],L[1],ori2,L[2],L[3],S))
            inl = inp.readline()
        
          
        
    #
    # Well... get_strech2 is problematic when overlapping sequences are needed.
    # Now this method is much smarter... and much shorter too.
    # 07/19,05. New lines are not formed in some entries. Fixed?
    # 04/27,06. Reverse complement coordinates are not reflected. Fixed.
    #           Modify this for multiple sequences
    #
    # @para fasta  file with single sequence
    # @para coords coord file, tab-delim, the last two token have to be coords.
    #
    def get_stretch3(self,fasta,coords):
        
        print("Covert fasta into a string...")
        inp = open(fasta)
        inl = inp.readlines()
        idx = self.rmlb(inl[0])[1:]
        seq = string.joinfields(inl[1:],"")

        # get rid of line breaks, sometimes \n and \r\n are mixed.
        if seq.find("\r\n") != -1:
            seq = string.joinfields(seq.split("\r\n"),"")
        if seq.find("\n") != -1:
            seq = string.joinfields(seq.split("\n"),"")
        print(len(seq))
        
        print("Read coords and output seq...")
        inp = open(coords)
        oup = open(coords+".fa","w")
        inl = inp.readline()
        c = 0
        while inl != "":
            print([inl])
            if c % 10000 == 0:
                print(" %i x 10k" % (c/10000))
            c += 1
            L = inl.split("\t")
            
            # deal with reverse ori
            cL = int(L[-2])
            cR = int(L[-1])
            ori = 1
            if int(L[-2]) > int(L[-1]):
                cL = int(L[-1])
                cR = int(L[-2])
                ori = -1
            
            S = seq[cL-1:cR]
            if S == "":
                print("ERR COORD: [%i,%i]" % (cL,cR))
            else:
                if ori == 1:
                    oup.write(">%s|%i-%i\n%s\n" % (idx,cL,cR,S))
                else:
                    S = trans.rc(S)
                    oup.write(">%s|%i-%i\n%s\n" % (idx,cR,cL,S))
                    #print ">%s|%i-%i %s" % (idx,cR,cL,S[:10])
                        
            inl = inp.readline()
        
        
        

    ##
    # coord format: from id<space>L-R to id_((L-1)*3+1)_R
    ## 
    def convert_header(self,fasta):

        inp = open(fasta,"r")
        oup = open(fasta+".coord.fa","w")
        inl = inp.readline()
        while inl != "":
            if inl[0] == ">":
                if inl.find(" ") == -1 or inl.find("-") == -1:
                    print("Wrong descriptor format, Quit!")
                    sys.exit(0)
                    
                L = inl.split(" ")
                left = int(L[1][:L[1].find("-")])
                righ = int(L[1][L[1].find("-")+1:-1])

                left = (left-1)*3+1
                righ = righ*3
                oup.write("%s_%i_%i\n" % (L[0],left,righ))

            else:
                oup.write(inl)
            inl = inp.readline()

        
    # rid of coords
    #
    # @style  coord connected to id via underscore: 0 [default], or space [1]
    def delete_coord(self,fasta,style):

        print("NOT SURE THIS IS WORKING YET")

        inp = open(fasta,"r")
        oup = open(fasta+".mod.fa","w")

        inl = inp.readline()
        while inl != "":
            if inl[0] == ">":
                if style == 0:
                    if inl.find("_") == -1:
                        oup.write(inl)
                elif style == 1:
                    if inl.find(" ") == -1:
                        oup.write(inl)


            inl = inp.readline()


    ##
    # Return a dict with seq id as key, seq as value. Will convert things to
    # UNIX format without "\r".
    #
    # 11/15 Major rewrite of the method
    # 11/17 Deal with redundant id, before they are simply replaced.
    #
    # @param fasta
    # @param desc    include desc[1] or not[0,default]. The description is
    #                separated from seq id with a space. If include, it will
    #                be the first element of the list, the second will be seq
    # @param verbose display line count [1], or not [0, default]
    # @param dflag   delimit by space or not, default no
    ##
    def fasta_to_dict(self,fasta,dflag=0,verbose=0):
    
        # each idx should only occur once.
        inp    = open(fasta,"r")
        inl = inp.readline()
        fdict  = {} # idx as key, seq as value
        c = 0
        N = 0 
        while inl != "":
            inl = self.rmlb(inl)
            next = 0
            idx  = ""
            desc = ""
            
            if inl[0] == ">":
                if verbose and c%10000 == 0:
                    print(" %i x 10k" % (c/10000))
                c += 1
                # rid of anything after space if asked
                if dflag and inl.find(" ") != -1:
                    desc = idx[idx.find(" ")+1:]
                    idx = idx[:idx.find(" ")]
                else:
                    idx = inl[1:]
                # count lines and store seq into a list
                slist = []
                inl = inp.readline()
                while inl[0] != ">":
                    slist.append(self.rmlb(inl))
                    inl = inp.readline()
                    if inl == "":
                        break
                seq = "".join(slist)
                
                if idx in fdict:
                    if verbose:
                        print("Redundant_id:",idx, end=' ')
                    if dflag:
                        if len(fdict[idx][1]) < len(seq):
                            fdict[idx] = [desc,seq]
                            if verbose:
                                print("longer")
                        else:
                            if verbose:
                                print("shorter")
                    else:
                        if len(fdict[idx]) < len(seq):
                            fdict[idx] = seq
                            if verbose:
                                print("longer")
                        else:
                            if verbose:
                                print("shorter")
                else:
                    N += 1
                    if dflag:
                        fdict[idx] = [desc,seq]
                    else:
                        fdict[idx] = seq
                next = 1        
            
            # so no extra line is read, because of the innder while
            if not next:
                inl = inp.readline()
        inp.close()
        if verbose:
            print("Total %i sequences, %i non-redundant" % (c,N))
        
        return fdict
        


    ##
    # Generate fasta file based on species passed, assuming desc like:
    #  >XX_id
    #
    # @param fasta
    ##
    def get_sp(self,fasta,species):

        print("Get sequences for a particular species:")
        print(" Fasta  :",fasta)
        print(" Species:",species)
        
        inp = open(fasta,"r")
        oup = open(fasta+".%s.fa" % species,"w")
        inl = inp.readline()
        slist = []

        w = 1 # write or not
        c = 0
        while inl != "":
            if inl[0] == ">":
                if inl[1:3] == species:
                    w = 1
                    c += 1
                else:
                    w = 0

            if w:
                oup.write(inl)       
            inl = inp.readline()

        print("Total %i sequences" % c)

    #
    # New rename function. But this function is not so good. It's slow and
    # only output seq if new name is found.. Should throw this away.
    #
    def rename2(self,fasta,name,ignore=""):
        
        print("Read fasta file...")
        fdict = manager.fasta_to_dict(fasta,0)
                
        print("Rename sequence...")
        inp = open(name,"r")
        inl = inp.readline()
        oup = open(fasta+"_rename.fa","w")
        absent = []
        ndict  = {}
        countN = 0
        countS = 0
        
        while inl != "":
            countN += 1
            inl = self.rmlb(inl).split("\t")
            if inl[1] in fdict:
                countS += 1
                oup.write(">%s\n%s\n" % (inl[0],fdict[inl[1]]))

                if inl[1] in ndict:
                    ndict[inl[1]]+= 1
                else:
                    ndict[inl[1]] = 1
            else:
                absent.append(inl[1])
            
            inl = inp.readline()

        
        oup = open(fasta + "_rename.log","w")
        oup.write("In name file but not in sequence file:\n ")
        oup.write("%s" % string.joinfields(absent,"\n "))
        
        oup.write("\nSequence id redundant:\n ")
        for i in ndict:
            if ndict[i] > 1:
                oup.write(" %s: %i\n" % (i,ndict[i]))
        
        print("Fasta: %i entries, %i with new name" % (len(list(fdict.keys())),countS))
        print("Name : %i entries, %i unique, %i not in fasta" % \
                (countN,len(list(ndict.keys())),len(absent)))
        
        
        

    ##
    # Rename id within fasta file
    #
    # @param fasta
    # @param name  name file in [new][old] format. Tab-delimited.
    # @param ignore Ignore whatever is after the passed char. Default nothing.
    #
    ##
    def rename(self,fasta,name,ignore=""):

        # put names to a dict, old as key, new as value
        inp    = open(name,"r")
        inl = inp.readline()
        ndict  = {}
        countR = 0
        found = 0
        
        print("Read name file...")
        oup_log = open(fasta+"_rename.log","w")
        oup_log.write("Redundant names:\n")
        while inl != "":
            inl = self.rmlb(inl)
            L = inl.split("\t")
            if L[1].lower() in ndict:
                oup_log.write(" %s\n" % L[1])
                found = 1
                countR = countR+1
            else:
                ndict[L[1].lower()] = L[0]
            inl = inp.readline()
        if not found:
            print(" none")

        # scan fasta descriptor line
        print("Read fasta file and change names...")
        inp    = open(fasta,"r")
        oup    = open(fasta+"_rename.fa","w")
        inl = inp.readline()
        countF = 0  # found
        countN = 0  # not found
        oup_log.write("\nFasta name not found:\n")
        while inl != "":
            if inl[0] == ">":
                inl = self.rmlb(inl)
                if ignore == "" or inl.find(ignore) == -1:  
                    if inl[1:].lower() in ndict:
                        oup.write(">%s\n" % ndict[inl[1:].lower()])
                        countF = countF+1
                    # if no new name, just output old ones.
                    else:
                        oup.write(inl+"\n")
                        oup_log.write(" %s\n" % inl[1:])
                        countN = countN+1
                else:
                    if inl[1:inl.find(ignore)].lower() in ndict:
                        oup.write(">%s%s\n" % \
                                 (ndict[inl[1:inl.find(ignore)].lower()],\
                                  inl[inl.find(ignore):]))
                        countF = countF+1
                    # if no new name, just output old ones.
                    else:
                        oup.write(inl+"\n")
                        oup_log.write(" %s\n" % inl[1:])
                        countN = countN+1
                    
            else:
                oup.write(inl)
            inl = inp.readline()
            
        print("Found %s, not found %s, %i redundant" % (countF,countN,countR))
        

    ##
    # Convert fasta seq id back to its original ones
    #
    # @param fasta      the fasta desc line can have just ">name\n" or something
    #          like ">name description". The script looks for space as
    #          delimiter in this case, so name should NOT contain any
    #          of it.
    # @param desc_flag  should see string after space as desc [1] or not [0]
    # @param name_file  in the format [new_name][old_name]
    ##
    def change_names(self,fasta,name_file,desc_flag,delim=" "):

        print("Change ID in Fasta file")
        print(" Fasta:",fasta)
        print(" Name :",name_file)
        
        # read names into a dict
        ndict = f_util.file_to_dict(name_file,5)
        print(ndict)
        sys.exit()

        #print ndict
        
        # compare names against the fasta file
        inp = open(fasta,"r")
        oup = open(fasta+".mod.fa","w")
        inl = inp.readline()
        countF = 0
        countA = 0
        while inl != "":
            if inl[0] == ">":
                countA = countA + 1

                old  = inl[1:-1]
                desc = ""
                if desc_flag and inl.find(" ") != -1:
                    old  = old[:old.find(" ")]
                    desc = inl[inl.find(" ")+1:-1]
                    print(old,">>",desc)
                
                if old in ndict:
                    if desc != "":
                        oup.write(">%s %s\n" % (ndict[old],desc))
                    else:
                        oup.write(">%s\n" % ndict[old])
                    countF = countF+1
                else:
                    oup.write(inl)
                    print("Unknwon id:",inl[1:-1])
            else:
                oup.write(inl)       

            inl = inp.readline()

        print("Total %i sequences, %i with new names" % (countA,countF))
        
        

    ##
    # Convert names into indices in a fasta file
    #
    # @param start  the starting index
    ##
    def index_names(self,fasta,start):

        print("Index name for:")
        print(" Fasta:",fasta)
        print(" Start index:",start)
        inp    = open(fasta,"r")
        oup1   = open(fasta+".fa","w")
        oup2   = open(fasta+".name","w")
        inl = inp.readline()
        c      = start
        ndict  = {}
        while inl != "":
            if inl[0] == ">":
                seq_id = inl[1:-1]
                if seq_id not in ndict:
                    ndict[seq_id] = c
                    oup1.write(">%i\n" % c)
                    oup2.write("%s\t%i\n" % (seq_id,c))
                    c = c+1
                else:
                    print("Redundant:",seq_id)
            else:
                oup1.write(inl)

            inl = inp.readline()

        print(" Ending index:",c)
        
        return ndict

    ##
    # 
    # @param fasta
    # @param return_dict boolean flag for deciding if file output [0, default]
    #                    should be generated or a dict should be returned [1].
    # @param x           determine also the masked sequences. Only works for 
    #                    text output but not direct call.
    ##
    def get_sizes(self,fasta,return_dict=0,x=0):

        inp = open(fasta,"r")
        inl= inp.readline()
        oup = open(fasta+".size","w")

        S  = 0
        X  = 0
        ID = ""
        sdict = {}
        T  = 0 # total characters
        while inl != "":
            if inl[-2:] == "\r\n":
                inl = inl[:-2]
            elif inl[-1] == "\n":
                inl = inl[:-1]          
            
            # just in case there is empty lines between sequences
            if inl == "" or inl.startswith("#"):
                pass
            elif inl[0] == ">":
                # output size, reset variables
                if S != 0:
                    if return_dict:
                        sdict[ID] = S
                    else:
                        if x:
                            pMasked = float(X)/float(S)
                            oup.write("%s\t%i\t%i\t%f\n" %(ID,S,X,pMasked))
                        else:
                            oup.write("%s\t%i\n" %(ID,S))
                        T += S
                ID = inl[1:]
                S  = 0
                X  = 0
            else:
                # Don't count stop
                if inl[-1] == "*":
                    inl = inl[:-1]
                S = S + len(inl)
                X = X + len(inl.split("x")) - 1
                
            inl= inp.readline()
        
        # write the last entry
        T += S
        print("Total length:",T)
        if return_dict:
            sdict[ID] = S
            return sdict
        else:
            if x:
                pMasked = float(X)/float(S)
                oup.write("%s\t%i\t%i\t%f\n" %(ID,S,X,pMasked))
            else:
                oup.write("%s\t%i\n" %(ID,S))
            

    ##
    # @param fasta  the fasta sequence file. Should include a sequence called
    #               OUT if the include_o is set to 1.
    # @param group  the tab-delimited file with group specification
    # @param create_d create dir for group. Default no [0]
    # @param includ_o include outgroup sequence or not. Default 0, no outgroup
    #               included.
    ##
    def get_group_seq(self,fasta,group,create_d,includ_o):
        
        # tdict = {group:{gene:[rest]}}
        tdict = f_util.file_to_dict(group,7)
        fdict = self.fasta_to_dict(fasta)
        
        if includ_o and "OUT" not in fdict:
            print("No sequence called OUT but you want to include outgroup.")
            print("QUIT!")
            sys.exit(0)
        
        print("Get sequences for group:\n")
        keys = list(tdict.keys())
        keys.sort()
        # iterate group
        for i in keys:
            print(" ",i)
            fname = i           
            if(create_d):
                os.system("mkdir %s" % i)
                fname = "./%s/%s" % (i,i)
                    
            oup = open(fname+".fa","w")
            for j in tdict[i]:
                if j in fdict:
                    # notice that I only take the 2nd element of fdict[j].
                    # the first contain some desc-like stuff, not added.
                    oup.write(">%s\n%s\n" % (j,fdict[j]))
                else:   
                    print("   %s missing" % j)  
                    
            if includ_o:
                oup.write(">OUT\n%s\n" % fdict["OUT"][1])
                        
            oup.close()
        
        # this is for direct method call in AlnUtility.batch_tree()
        return keys
            
    #
    # A simpler get_sequences without reading seq into dict
    #
    # @param fasta
    # @param name   single column or tab delimited file
    # @oaram tokens the tokens to regarded as ID, separated by ',', default ""
    #   
    def getseq2(self,fasta,name,tokens=""):
        
        print("Read name file:",name)
        if tokens == 0 or tokens == "":
            tokens = [0]
        else:
            tmp = tokens.split(",")
            tokens = []
            for i in tmp:
                tokens.append(int(i))
        
        inp = open(name)
        inl = inp.readline()
        ndict = {}
        countR = 0
        rlist  = []
        while inl != "":    
            L = self.rmlb(inl).split("\t")
            for j in tokens:
                if L[j] == "":
                    pass
                elif L[j] not in ndict:
                    ndict[L[j]] = 0
                else:
                    countR += 1
                    rlist.append(L[j])
            inl = inp.readline()
        print(" %i redundant"    % countR)
        countN = len(list(ndict.keys()))
        print(" %i unique names" % countN)
        
        print("\nRead fasta file:")
        oup = open(name+".fa","w")
        countA = 0
        countF = 0
        for i in fasta.split(","):
            print("",i)
            inp = open(i)
            inl = inp.readline()
            flagW  = 0
            while inl != "":
                if ">" in inl:
                    if countA % 1e4 == 0:
                        print(" %i x 10k" % (countA/1e4))
                    countA += 1
                    # number found is the same as number of unique names
                    if countF == countN:
                        break                   
                    N = self.rmlb(inl)[1:]
                    # in dict and has not been written
                    if N in ndict and ndict[N] == 0:
                        #print " ",N
                        ndict[N] = 1
                        countF += 1
                        flagW = 1
                    else:
                        flagW = 0
                if flagW:
                    oup.write(inl)

                inl = inp.readline()
                
        print(" %i found\n" % countF)
        oup.close()
    
        if countF < len(list(ndict.keys())):
            print("Missing:", end=' ')
            countM = 0
            for i in ndict:
                if ndict[i] == 0:
                    #print "",i
                    countM += 1
            print(countM)
        
            
    ##
    # All names are coverted to lower case when compared but the output file
    # retain the same case as the fasta file. If the fasta file sequence header
    # contain space, the sequence id is whatever after ">" and  before the space.
    #
    # @param segment  no[0,default], yes[1]. This is done based on the name fie
    #        where coordinates are specified like:
    #          >name [whatever] L-R
    #        The method uses space as delimiter. The last token is
    #        assumed to be the one with coord info.
    # @param type     fasta seqeunce is protein [0,default] or nucleotide [1]
    # @param match    exact match to source fasta [0, default], or just match
    #        the passed name to the exact number of char after ">" in
    #        fasta file [1]
    # @param call     0 (default): run whole module, or 1: method call from other
    #        module.
    ##
    def get_sequences(self,fasta,name_file,segment=0,type=0,match=0,call=0,
                      outname=""):
        
        print("Fasta:",fasta)
        print("Name :",name_file)
        
        # construct name list
        if not call:
            inp     = open(name_file,"r")
            inls = inp.readlines()
            names   = {}
            print("Redundant names in name file:")
            flag = 0
            for i in inls:
                # check format
                if i[-2:] == "\r\n":
                    i = i[:-2]
                elif i[-1] == "\n":
                    i = i[:-1]
                    
                # skip empty line
                if i == "":
                    continue
                
                # see if anything else is encoded in the name
                if i.find(" ") != -1:
                    llist = i.split(" ")
                    #llist[-1] == llist[1][:-1] disabled due to format checking
                else:
                    llist = [i]
                    #llist = [i[:-1]], modified due to format checking

                if llist[0].lower() not in names:
                    if len(llist) == 1:
                        names[llist[0].lower()] = [0,""]
                    else:
                        names[llist[0].lower()] = [0,llist[1:]]
                else:
                    print(" ",i)
                    flag = 1
            
            if not flag:
                print(" None...")

            print("Total %i names" % len(list(names.keys())))
        else:
            names = name_file
        
        if outname != "":
            name_file = outname
        
        countP = 0 # counting problem coord
        countO = 0 # counting ok coord
        
        # scan fasta file
        inp     = open(fasta,"r")              
        oup     = open(name_file+".fa","w")
        inl  = inp.readline()

        count = 0
        flag = 0
        nkeys = list(names.keys())
        #print nkeys
        key_found = []
        while inl != "":
            found = 0
            if inl[0] == ">":
                #print inl[:-1]
                
                # convert id to lower case
                if inl.find(" ") != -1:
                    id = inl[1:inl.find(" ")].lower()
                else:
                    id = inl[1:-1].lower()
                    if id[-1] == "\r":
                        id = id[:-1]
                
                # check if any id match the same length string after ">", slow
                gotit = 0
                nname = ""
                if match:
                    for i in nkeys:
                        #print "",i,"vs",id[:len(i)]
                        if i == id[:len(i)]:
                            nname = i
                            gotit = 1
                            #print "  found!"
                            break
                
                # include modification that won't take "redundant" sequence.
                # but this modification may impact the codes getting domains
                # need to see if this is the case.
                if (id in names and names[id][0] != 1) or gotit:
                    count = count +1
                    #print "Found:",count,id
                    
                    if segment == 0:
                        # take care of two domains in the same sequecnce
                        if names[id][0] > 0:
                            if inl[-2:] == "\r\n":
                                inl = inl[:-2]  
                            elif inl[-1] == "\n":
                                inl = inl[:-1]
                            oup.write("%s %i\n" % (inl,names[id][0]))                       
                        else:
                            oup.write(inl)

                        inl = inp.readline()
                        while inl[0] != ">":
                            oup.write(inl)
                            inl = inp.readline()                     
                            if inl == "":
                                break

                    elif segment > 0:
                        coord = names[id][1][-1]
                        L     = int(coord[:coord.find("-")])
                        R     = int(coord[coord.find("-")+1:])
                        
                        if (abs(R-L)+1)%3 != 0:
                            countP = countP+1
                            #print "",id,L,R
                        else:
                            countO = countO+1

                        # take care of two domains in the same sequecnce
                        if names[id][0] > 0:
                            if inl[-2:] == "\r\n":
                                inl = inl[:-2]  
                            elif inl[-1] == "\n":
                                inl = inl[:-1]
                            oup.write(">%s %i %s\n" % (id,names[id][0],coord))                      
                        else:
                            oup.write(">%s %s" % (id,coord))

                        inl = inp.readline()
                        seq    = ""
                        while inl[0] != ">":
                            seq = seq + inl[:-1]
                            inl = inp.readline()             
                            if inl == "":
                                break                      

                        if R > L:
                            if type:
                                seq = seq[L-1:R-3]
                            else:
                                seq = seq[L-1:R-1]
                        else:
                            if type:
                                seq = seq[R-1:L-3]
                            else:
                                seq = seq[R-1:L-1]

                        oup.write(seq+"\n")
    
                    # set names value to 1
                    if match:
                        names[nname] = [1,names[nname][1]]
                    else:
                        names[id] = [1,names[id][1]]
                    found = 1
                    

            if not found:
                inl = inp.readline()

        # output a list of sequences not found
        oup.close()
        oup = open(name_file+".log","w")
        if not call:
            for i in list(names.keys()):
                if names[i][0] == 0:
                    flag = 1
                    oup.write("%s\n" % i)

            if not flag:
                print(" None...")
            print("Found %i out of %i" % (count,len(list(names.keys()))))
            if type:
                print("%i are in frame, %i out of frame" % (countO,countP))

    ##
    # Convert fasta style file to one line format with [seq_id][seq]
    #
    # @param format this is only for direct method call, not from command line.
    #          0: default, output as a file
    #          1: output as a dict, with id as key, seq as value
    #        d delimiter, defailt the whole line will be taken.
    ##
    def fasta_to_oneline(self,fasta,format=0,d=""):
        
        #print "Convert Fasta file to one line format:"
        
        inp = open(fasta,"r")
        
        inl = inp.readline()
        count = 0
        odict = {}
        seq   = ""
        id = ""
        while inl != "":
            if inl[-2:] == "\r\n":
                inl = inl[:-2]
            elif inl[-1] == "\n":
                inl = inl[:-1]
            if inl == "":
                inl = inp.readline()
                continue
        
            found = 0
            if inl[0] == ">":
                
                if seq != "":
                    odict[id] = seq
                    seq = ""
                
                if d != "" and inl.find(d) != -1:
                    id = inl[1:inl.find(d)]
                else:
                    id = inl[1:]

                odict[id] = ""
                inl = inp.readline()
                while inl[0] != ">":
                    if inl[-2:] == "\r\n":
                        inl = inl[:-2]
                    elif inl[-1] == "\n":
                        inl = inl[:-1]

                    seq = seq + inl
                    inl = inp.readline()
                    
                    if inl == "":
                        break

                found = 1
                count = count +1        

            if not found:
                inl = inp.readline()

        # put the last seq into dict
        odict[id] = seq

        if format:
            return odict
        else:
            oup = open(fasta+".pep","w")
            okeys = list(odict.keys())
            okeys.sort()
            for i in okeys:
                oup.write("%s\t" % i)
                oup.write("%s\n"  % odict[i])
                
        print("Total %i sequences converted" % count)


    def oneline_to_fasta(self,oneline):

        inp = open(oneline,"r")
        oup = open(oneline+".fa","w")
        inl = inp.readline()
        while inl != "":
            llist = self.rmlb(inl).split("\t")

            oup.write(">"+llist[0]+"\n")

            while len(llist[1]) > 0:
                if len(llist[1]) > 80:
                    oup.write(llist[1][:80]+"\n")
                else:
                    oup.write(llist[1][:80])
                try:
                    llist[1] = llist[1][80:]
                except IndexError:
                    break
            oup.write("\n")
            inl = inp.readline()

        print("Output to %s.fa" % oneline)
        
    
    
    ##
    # Parse the description out of fasta.
    #
    # @param fasta  the fasta file with meaningful header
    # @param style  the type of fasta file:
    #       GenBank(gb) : gi|xxxxx|db|acc| desc [species]
    #       MAtDB(matdb): gene_name desc
    #       EnsEMBL(ensembl): [id] [Gene:] [Clone:] [Contig:] [Chr:] 
    #                [Basepair:] [status:]
    ##
    def parse_desc(self,fasta,style,delim):
        
        print("Fasta:",fasta)
        print("Style:",[style])
        print("Delim:",[delim])
        
        inp = open(fasta,"r")
        oup = open(fasta+".desc","w")
        
        #if style == "gb":
        #       oup.write("GI#\tdatabase\tAcc\tAcc2\tDescription\tOrganism\n")
        #elif style == "ensembl":
        #       oup.write("ID\tGene\tClone\tContig\tChr\tBasePair\tStatus\n")
        
        print("\nRead through fasta file...")
        inl = inp.readline()
        c = 0
        while inl != "":
            if inl[0] == ">":
                if c % 1000 == 0:
                    print(" %i k" % (c/1000))
                c += 1
                    
                # rid of ">" and "\n"
                inl = self.rmlb(inl)[1:]
                desc = species = ""

                # GenBank style
                if style == "gb":
                    glist  = inl[:inl.find(" ")]
                    glist  = glist.split("|")                      
                    if glist[-1] != "" and inl.find("[") == -1:
                        # PIR style header
                        if glist[2] == "pir":
                            glist[3] = glist[4]
                            desc     = inl[inl.find(" ")+1:
                                          inl.find(" - ")]
                            species  = inl[inl.find(" - ")+3:]
                        # GB screwup type
                        elif glist[2] == "gb":
                            desc     = inl[inl.find(" ")+1:]
                            species  = "unknown"
                        # No idea where this come form
                        elif glist[2] == "sp":
                            desc     = inl[inl.find(" ")+1:]   + " (" +\
                                       glist[4][:glist[4].find("_")] + ")"
                            species  = glist[4][glist[4].find("_")+1:]
                        # all other problems
                        else:
                            print("PROBLEMATIC ENTRY:")
                            print(inl)          
                    else:                                          
                        desc   = inl[inl.find(" ")+1:inl.find("[")-1]
                        species= inl[inl.find("["):]
                    
                    # this will be: "gi",gi,db,acc,???
                    oup.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (glist[1],glist[2],
                                                glist[3],glist[4],desc,species))    
                # use 1st space as delimiter style
                elif style == "1stspace":
                    L = inl.split(" ")
                    if len(L) > 1:
                        oup.write("%s\t%s\n" %
                                        (L[0],string.joinfields(L[1:],"\t")))
                    else:
                        oup.write("%s\n" % inl)
                elif style == "":
                    if delim != "":
                        L = string.joinfields(self.rmlb(inl).split(delim),"\t")
                        oup.write("%s\n" % L)
                    else:
                        oup.write("%s\n" % inl)
                # all others
                else:
                    L = self.rmlb(inl).split(" ")
                    oup.write("%s\n" % string.joinfields(L,"\t"))
                    
            inl = inp.readline()
        

    def parse_ensembl_fasta(self,fasta,sp):

        print("Parse ensembl fasta:",fasta)
        inp = open(fasta,"r")
        oup = open(fasta+".ensembl","w")
        inl = inp.readline()
        seq = ""
        c   = 0
        while inl != "":
            if inl[0] == ">":
                if seq != "":
                    oup.write("%i\t%s\n" % (len(seq),seq))
                    
                llist = inl[1:-1].split(" ")
                
                if len(llist) == 1:
                    if llist[0].find("SINFRU") != -1:
                        index = llist[0][llist[0].find("SINFRU")+6:]
                        llist.append("SINFRUG"+index)                      
                        llist.extend(["","","0","0",""])
                    elif sp == "ce":
                        try:
                            int(llist[0][-1])
                            llist.append(llist[0])
                        except ValueError:
                            llist.append(llist[0][:-1])
                        llist.extend(["","","0","0",""])
                    else:
                        print("Line token number mismatch, quit!")
                        sys.exit(0)

                elif sp == "dm":
                    tlist = []
                    tlist.append(llist[1][llist[1].find(":")+1:]) # name
                    tlist.append(llist[3][llist[3].find("[")+1:]) # gene
                    tlist.extend(["",""])
                    clist = llist[8].split(":")
                    tlist.append(clist[1][1:])          # chr
                    tlist.extend(["0",""])
                        
                    llist = tlist
                    #print llist
                        
                if sp == "ce" or sp == "dm":
                    id = c
                else:
                    id = int(llist[0][-8:])
                
                #print llist
                oup.write("%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
                         (id,                 # id
                          llist[0],         # name
                      llist[1][llist[1].find(":")+1:], # gene
                      llist[2][llist[2].find(":")+1:], # clone
                      llist[3][llist[3].find(":")+1:], # contig
                      llist[4][llist[4].find(":")+1:], # chr
                      llist[5][llist[5].find(":")+1:], # basepair
                      llist[6][llist[6].find(":")+1:]))# status
                seq = ""
                c += 1
            else:
                seq = seq + inl[:-1]
            
            inl = inp.readline()
        

    #
    # This is for getting the longest transcript or protein seq of the same
    # gene.
    #
    def ensembl_longest(self,fasta):
        # assume fasta header is in the format:
        # >pep_id pep:xxx chromosome:xxx gene:xxx transcript:xxx
        
        # put fasta into dict
        print("Read fasta into dict:")
        inp = open(fasta)
        fd  = {}
        inl = inp.readline()
        seq = gid = pid = ""
        countG = 0
        countP = 0
        while inl != "" or seq != "":
            if inl == "" and seq != "" or inl[0] == ">":
                if seq != "":
                    if gid == "" or pid == "":
                        print("Problem g or pid:",[gid,pid])
                        sys.exit(0)
                    if gid in fd:
                        countP += 1
                        if pid in fd[gid]:
                            print("This should not happen:",pid)
                        else:
                            fd[gid][pid] = seq
                    else:
                        countG += 1
                        countP += 1
                        fd[gid] = {pid:seq}
                    # reset seq
                    seq = ""
                    # just put the last entry into dict, out of the loop
                    if inl == "":
                        break
                # set gid and pid
                L = self.rmlb(inl[1:]).split(" ")
                pid = L[0]
                for j in L[1:]:
                    if j.find("gene") != -1:
                        gid = j[j.find(":")+1:]
                        break
                #print gid,pid
            else:
                seq += self.rmlb(inl)
            inl = inp.readline()
        
        print(" %i genes, %i peptides" % (countG,countP))
        
        # go through each gene and output the longeset one
        print("Get the longest one...")
        oup = open(fasta+".long","w")
        for i in fd:
            longest = []
            for j in fd[i]:
                if longest == [] or longest[1] < len(fd[i][j]):
                    longest = [j,len(fd[i][j])]
            oup.write(">%s\n%s\n" % (longest[0],fd[i][longest[0]]))
        
          
    
    #
    # Delimiter: if "\t", should send "tab"
    #
    def simplify_desc(self,fasta,style,delim,tokens):

        print("Process:%s" % fasta)
        print("Output :%s.mod.fa" % fasta)
        print("Delim  :",[delim])
        print("Tokens :",[tokens])
        
        if delim in ["\\t","tab"]:
            delim = "\t"
        
        if tokens != "":
            tmp = tokens.split(",")
            tokens = []
            for i in tmp:
                tokens.append(int(i))
        
        print("Style  :%s\n" % style)
        inp = open(fasta,"r")
        oup = open(fasta+".mod.fa","w")
        inl = inp.readline()
        count = 0
        while inl != "":
            inl = self.rmlb(inl)
            if inl != "" and inl[0] == ">":
                if count % 1e3 == 0:
                    print(" %i k" % (count/1e3))
                count = count+1
                if style == "gb":
                    if inl.find("|") == -1:
                        print(inl)
                    llist = inl.split("|")
                    inl = ">" + llist[1]
                elif style == "ensembl":
                    inl = inl[:inl.find(" ")]
                    if inl.find(":") != -1:
                        inl = ">" + inl[inl.find(":")+1:]
                elif delim != "":
                    idx = inl.find(delim)
                    # if delimiter does not exist, just don't change anything
                    if idx == -1:
                        pass
                    else:
                        # no token specified just take the 1st one
                        if tokens == "":
                            inl = inl[:idx]
                        else:
                            L = inl.split(delim)
                            tmp = []
                            for j in tokens:
                                tmp.append(L[j])
                            inl = string.joinfields(tmp,"|")
                            if inl[0] != ">":
                                inl = ">"+inl
                else:
                    inl = inl[:inl.find(" ")]
            
            if inl != "":       
                oup.write(inl+"\n")            
            
            inl = inp.readline()

        print("Total %i entries" % count)
        

    ##
    # @param fasta   the seq file to change
    # @param prefix  species code
    ##
    def add_prefix(self,fasta,prefix):
    
        print("Add prefix...")
        inp = open(fasta,"r")
        oup = open(fasta+".mod.fa","w")
        inl = inp.readline()
        while inl != "":
            if inl[0] == ">":
                inl = ">" + prefix + "_" + inl[1:]

            if inl != "\n" and inl != "\r\n":
                if inl[-2:] == "\r\n":
                    inl = inl[:-2] + "\n"
                oup.write(inl)
                        
            inl = inp.readline()
            
        print("Output file %s generated." % (fasta+".mod.fa"))
        

        #Gaurav script: Add suffix
    def add_suffix(self,fasta,suffix):
    
        print("Add suffix...")
        inp = open(fasta,"r")
        oup = open(fasta+".mod.fa","w")
        inl = inp.readline()
        while inl != "":
            if inl[0] == ">":
                inl = inl[:-1]+suffix+"\n"

            if inl != "\n" and inl != "\r\n":
                if inl[-2:] == "\r\n":
                    inl = inl[:-2] + "\n"
                oup.write(inl)
                        
            inl = inp.readline()
            
        print("Output file %s generated." % (fasta+".mod.fa"))
        
    
    ##
    # This function assumes that each line of the fasta file is ended in "\n"
    # and that no empty line is present between sequences. 
    ##
    def size_filter(self,fasta,threshold,out_name):
        
        print("Retreive sequences with length above threshold...")
        
        if out_name == "":
            out_name = fasta + "_T%i" % threshold + ".fa"
        
        inp = open(fasta,"r")
        oup = open(out_name,"w")
        
        size = 0
        inl = inp.readline()
        llist = []
        sdict = {}
        while inl != "":             
            if inl[0] == ">":
                # write llist if not empty and size is above threshold
                if llist != []:
                    for i in llist[1:]:
                        if i[-2:] == "\r\n":
                            size = size + len(i) - 2
                        elif i[-1] == "\n":
                            size = size + len(i) - 1
                        else:
                            size = size + len(i)                            
                    if size in sdict:
                        sdict[size] = sdict[size]+1
                    else:
                        sdict[size] = 1
                    
                    if size >= threshold:
                        for i in llist:
                            oup.write(i)
                
                # reset values
                llist = []
                size  = 0
                
            llist.append(inl)           
            inl = inp.readline()
        
        # deal with the last sequence
        for i in llist[1:]:
            if i[-2:] == "\r\n":
                size = size + len(i) - 2
            elif i[-1] == "\n":
                size = size + len(i) - 1
            else:
                size = size + len(i)

        if size in sdict:
            sdict[size] = sdict[size]+1
        else:
            sdict[size] = 1
        
        if size >= threshold:
            for i in llist:
                oup.write(i)
        
        print("Fasta output %s generated" % out_name)
        
        out_name = out_name + ".stat"
        
        keys = list(sdict.keys())
        keys.sort()
        oup = open(out_name,"w")
        for i in keys:
            oup.write("%i\t%s\n" % (i,sdict[i]))
            
        print("Length statistics %s generated" % out_name)
             
    
    ##
    # Compare two lists and generate a output with ".comp" 
    #
    # @param files   a string of file names separated by ","
    # @param outname name for output file
    ##
    def compare_lists(self,files,outname):
        
        lists = files.split(",")
        list1 = lists[0]
        list2 = lists[1]
        
        # put list1 into dict
        inp = open(list1,"r")
        dict1 = {}
        inl = inp.readline()
        while inl != "":
            if inl[:-1].lower() not in dict1:
                dict1[inl[:-1].lower()] = 0
            inl = inp.readline()
            
        # put list2 into dict
        inp = open(list2,"r")
        dict2 = {}
        inl = inp.readline()
        while inl != "":
            if inl[:-1].lower() not in dict2:
                dict2[inl[:-1].lower()] = 0  
            inl = inp.readline()
            
        
        if outname == "":
            outname = "TEMP.COMP"
        oup = open(outname,"w")
        
        # compare dict1 against dict2
        keys1 = list(dict1.keys())
        keys2 = list(dict2.keys())
        list1not2 = 0
        for i in keys1:
            if i in dict2:
                oup.write(i+"\t"+i+"\n")
                dict2[i] = 1
            else:
                list1not2 = list1not2 +1
                oup.write(i+"\t-\n")
        
        # write the rest of dict2 not compared
        list2not1 = 0      
        for i in keys2:
            if dict2[i] == 0:
                list2not1 = list2not1 +1
                oup.write("-\t"+i+"\n") 
        
        print("In list1 not list2:",list1not2)
        print("   list2 not list1:",list2not1)


    ##
    # The extract_cds function in ParseBlast will leave the assembly info in
    # the sequnence within the brackets. Rid of these things and replace it
    # with "X"
    ##
    def cleanup(self,fasta):
        
        inp = open(fasta,"r")                              
        oup = open(fasta+"_clean.fa","w")
        
        inl = inp.readline()
        while inl != "":
            if inl.find("(") != -1:
                oup.write(inl[:inl.find("(")]+"X"+\
                          inl[inl.find(")")+1:])
            else:
                oup.write(inl)
            inl = inp.readline()


    #
    # Get rid of redundant sequences by name
    # 
    # 11/17,04 Major rewrite. Use fasta_to_dict
    #
    def del_redun(self,fasta):

        fdict = self.fasta_to_dict(fasta,0,verbose=1)
        fkeys = list(fdict.keys())
        fkeys.sort()
        oup   = open(fasta+".mod","w")
        for i in fkeys:
            oup.write(">%s\n%s\n" % (i,fdict[i]))
        
        

        ##
    # Get rid of redundant sequences by seqeunce
    # Created by David E. Hufnagel on Mar 25, 2013
    ##
    def del_redun_by_seq(self, fasta):
                
                oup = open(fasta+".mod","w")

                
                #go through fasta and make a dict of key: seq val: name
                fdict = self.fasta_to_dict(fasta,0,verbose=0)
                fdictRev = {}
                for name in fdict:
                        seq = fdict[name]
                        self.SaveIntoDict(seq, name, fdictRev)
                                        
                #go through dict and if there are multiple names in the val concatenate them, seperated by ";" and output the info
                outDict = {}
                for seq in fdictRev:
                        namesLst = fdictRev[seq]
                        namesStr = ";".join(namesLst)
                        newLines = ">%s\n%s\n" % (namesStr, seq)
                        oup.write(newLines)
                        

    
    def check_redun(self,dir):
        
        f = os.listdir(dir)
        for i in f:
            try:    
                print(i)
                self.fasta_to_dict(i,0,verbose=1)
                print("")
            except IndexError:
                continue
        
        


    ##
    # This is written to get rid of sequences within a particular fasta file
    #
    # @param fasta  fasta file with entries to be deleted
    # @param dlist  file with a list of names to be deleted, or a string with
    #               species header. ALWAYS look for 2 charcters
    ##
    def delete(self,fasta,dlist):
        
        dflag = 0
        try:
            print("Read delete list:")
            ddict = f_util.file_to_dict(dlist,0)
            print("To be deleted:",len(list(ddict.keys())))
        except IOError:
            print(" not a file.")
            dlist = dlist.split(",")
            dflag = 1
        
        
        inp = open(fasta,"r")
        oup = open(fasta+".mod","w")
        inl = inp.readline()
        write  = 1
        countD = 0
        countW = 0
        countT = 0
        while inl != "":
            inl = self.rmlb(inl)
            if inl == "":
                inl = inp.readline()
                continue
            
            if inl[0] == ">":
                countT += 1                 
                if dflag:
                    if inl[1:3] in dlist:
                        write = 0
                        countD += 1
                    else:
                        write = 1
                        countW += 1
                else:
                    if inl[1:] in ddict:
                        ddict[inl[1:]] = -1
                        write = 0
                        countD += 1
                    else:
                        write = 1
                        countW += 1
            if write == 1:
                oup.write(inl+"\n")
            
            inl = inp.readline()
        
        if dflag == 0:
            oup = open(fasta+".log","w")
            for i in ddict:
                if ddict[i] != -1:
                    oup.write("%s\n" % i)
        
        print("Total %i, write %i, deleted %i" % (countT,countW,countD))
        print("Deletion done!")
        
    #
    # Break fasta file into several different files.
    # 
    # @param by     divide by the int passed
    # @param setdir make dir for each file
    #   
    def divide(self,fasta,by,setdir=0):
        
        fdict = self.fasta_to_dict(fasta,0)
        fkeys = list(fdict.keys())
        fkeys.sort()
        leng  = len(fkeys)/by
        
        C = 0
        for i in range(by):
            oup  = open(fasta+"_%i" % (i+1),"w")
            if C == by -1:
                keys = fkeys[leng*C:]   
            else:
                keys = fkeys[leng*C:leng*(C+1)]         
            
            for j in keys:
                oup.write(">%s\n%s\n" % (j,fdict[j]))
            
            if setdir:
                os.system("mkdir %i" % i)
                os.system("mv %s_%i %i" % (fasta,i,i))
            C = C+1
            
        

    #
    # Concatenate sequences in a file together and record the coordinates
    #
    def concat(self, fasta):
        
        fdict = self.fasta_to_dict(fasta)
        oup1  = open(fasta+".concat","w")
        oup1.write(">%s.concat\n" % fasta)
        oup2  = open(fasta+".concat.coord","w")
        
        c = 0
        for i in fdict:
            seq = self.rmlb(fdict[i][1])
            
            oup1.write("%s" % seq)
            oup2.write("%s\t%i\t%i\n" % (i,c+1,c+len(seq)))
            c = c + len(seq)
        
        
    
    #
    # Extension of concat. Use in the following situation, for example:
    #
    # Several rice contig sequences are concatenated together. Then genes are
    # predicted from such concatenated sequence file. Now we want to map the
    # genes back to the contigs.
    #
    # @param coord1 the coord output file of concat. L should be smaller than R.
    # @param coord2 a file with [seq][L][R]. The coordinates should be within
    #               the range of coord1. L does not have to be smaller than R.
    #
    def concat_locate(self, coord1, coord2):
        
        c1  = f_util.file_to_dict(coord1,3)
        c2  = f_util.file_to_dict(coord2,3)
        oup = open(coord2+".mapped","w")
        for i in c2:
            if int(c2[i][0]) < int(c2[i][1]):
                gL = int(c2[i][0])
                gR = int(c2[i][1])
            else:
                gL = int(c2[i][1])
                gR = int(c2[i][0])
            
            hit = "not_mapped"
            for j in c1:
                #print int(c1[j][0]),int(c1[j][1]),"--->",gL,gR
                if int(c1[j][0]) < gL and int(c1[j][1]) > gR:
                    hit = j
                    #print "HERE"
                    break
            oup.write("%s\t%s\n" % (i,hit))
            
        
    
    #
    # Store each fasta squence in a file into its own file
    #
    def indiv(self,fasta):
        
        inp = open(fasta)
        inl = inp.readline()
        print("Write seq...")
        while inl != "":
            if inl[0] == ">":
                seq = self.rmlb(inl[1:])
                seq = string.joinfields(seq.split(" "),"_")
                print("",seq)
                oup = open(seq+".fa","w")
                oup.write(">%s\n" % seq)
            else:
                oup.write(inl)
            inl = inp.readline()
        
        
    
    def count(self,fasta):
        inp = open(fasta)
        inl = inp.readline()
        print("Count seq...")
        c = 0
        while inl != "":
            if inl[0] == ">":
                c += 1
            inl = inp.readline()
        print("Total %i sequences" % c)
              

        #Created by David E. Hufnagel on Jan 9, 2013
        def cntUnique(self,fasta):
                inp = open(fasta)
                print("Count seq...")
                cntSet = set()
                for line in inp:
                        if line.startswith(">"):
                                name = line.strip(">").strip("\n")
                                cntSet.add(name)
                print("Total %s sequences" % (len(cntSet)))
        

    def get_longest(self,fasta,D):
        
        fa = self.fasta_to_dict(fasta)
        fk = list(fa.keys())
        gn = {} # gn = {gene_name:{ver:length}}
        countT = 0
        countG = 0
        for i in fk:
            countT += 1
            g = i.split(D)
            if len(g) != 2:
                print("Gene name format problem with this delimiter:")
                print("",g)
                print("Quit!")
                sys.exit(0)
            
            if g[0] not in gn:
                countG += 1
                gn[g[0]] = {g[1]:len(fa[i])}
            else:
                gn[g[0]][g[1]] = len(fa[i])
        print("Total %i pep belong to %i genes" % (countT,countG))
        
        oup1 = open(fasta+".longest","w")
        oup2 = open(fasta+".longest.log","w")
        for i in gn:
            # compare length of a gene_name
            maxLen = 0
            maxVer = ""
            for j in gn[i]:
                if gn[i][j] > maxLen:
                    maxLen = gn[i][j]
                    maxVer  = j
            maxID = "%s%s%s" % (i,D,maxVer)
            oup1.write(">%s\n%s\n" % (maxID,fa[maxID]))
            for j in gn[i]:
                flag = "-"
                if j == maxVer:
                    flag = "long"
                ID = "%s%s%s" % (i,D,j)
                oup2.write("%s\t%s\t%i\t%s\n" % (i,j,gn[i][j],flag))
        
        

    def rmlb(self,astr):
        if astr[-2:] == "\r\n":
            astr = astr[:-2]
        elif astr[-1] == "\n":
            astr = astr[:-1]
        return astr

    #condenses all sequence lines into one.  Created By David E. Hufnagel on Aug 24, 2012
    def CondenseSeq(self):
        #open files
        inpFd = open(fasta)
        outp = open(out_name, "w")
        
        #Write the users command line prompt on the first line of the output file.
        outp.write("#python %s\n" % (" ".join(sys.argv)))

        #do the actual processing
        tempSeq = ""
        lastName = ""
        for line in inpFd:
            if not line.startswith("#"):
                if line.startswith(">"):
                    if not tempSeq == "":
                        outp.write(lastName)
                        outp.write(tempSeq + "\n")
                    tempSeq = ""
                    lastName = line
                else:
                        tempSeq += line.strip()

    #removes all sequence lines containing the designated keyword.  Created By David E. Hufnagel on Oct 7, 2012
    def RemoveByKeyword(self):
        #open files
        inpFd = open(fasta)
        outp = open(out_name, "w")
        
        #Write the users command line prompt on the first line of the output file.
        outp.write("#python %s\n" % (" ".join(sys.argv)))
                        
        #do the actual processing
        doWrite = True
        for line in inpFd:
            if not line.startswith("#"):
                if line.startswith(">"):
                    if keyword in line:
                        doWrite = False
                    else:
                        outp.write(line)
                        doWrite = True
                else:
                    if doWrite == True:
                        outp.write(line)

    #basically my version of change_names.  Replaces a code name in a fasta file with the
    #real name (or vice versa).
    def ReplaceNames(self, reverse=False):
        #open files
        inpFd = open(fasta)
        nameFd = open(name_file)
        outp = open(out_name, "w")

        #Write the users command line prompt on the first line of the output file.
        outp.write("#python %s\n" % (" ".join(sys.argv)))

        #take the ref file into a dict
        refDict = {}
        for line in nameFd:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                if reverse == True:
                    refDict[lineLst[0]] = lineLst[1].strip()
                else:
                    refDict[lineLst[1].strip()] = lineLst[0]

        #go through the fasta and replace the names
        for line in inpFd:
            if not line.startswith("#python"):
                if line.startswith(">"):
                    newName = refDict[line[1:-1]]
                    outp.write(">%s\n" % (newName))
                else:
                    outp.write(line)

    #Calculates and prints n50
    #Created by David E. Hufnagel on Mar 30, 2013
    def N50(self):
        #make the .size file for the fasta
        #self.get_sizes(fasta)

        #go through the fasta .size file and make a list of tuples of (size, name) and keep track of the total of all sequence lengths
        size = fasta + ".size"
        sizeFD = open(size)
        sizeLst = []
        total = 0
        for line in sizeFD:
            if not line.startswith("#"):
                lineLst = line.strip().split("\t")
                sizeLst.append((int(lineLst[1]), lineLst[0]))
                total += int(lineLst[1])

        #sort the list
        sizeLst.sort()

        #go through the list from biggest seq to smallest (in reverse order)
        tempTotal = 0 #a dynamic count of the sizes thus far
        for group in sizeLst[::-1]:
            tempTotal += group[0]
            if tempTotal >= total / 2.0:
                print("n50: ", group[0])
                break

    #This script is designed to take a clustalw output fasta file and determine %ID
    #between all pairwise clusters, outputting in the format:
    #gene1  gene2  %ID
    #WARNING: CLUSTERS MUST BE IN GROUPS OF 2
    def GetPercentIDFromClustalW(self):
        inpFd = open(fasta)
        outp = open(out_name, "w")
        tempDict = {}
        name = ""
        for line in inpFd:
            if not line.startswith("#python"):
                lineLst = line.strip().split("\t")
                #1) get cluster info
                if line.startswith("###"):
                    #print "1\n"
                    #2)calculate %ID
                    if type:
                        perID = self.CalculatePerIDNuc(tempDict[lastName], tempDict[name])
                    else:
                        perID = self.CalculatePerIDProt(tempDict[lastName], tempDict[name])
                    #3) output info
                    newLine = "%s\t%s\t%s\n" % (lastName, name, perID)
                    outp.write(newLine)
                    tempDict = {}
                elif line.startswith(">"):
                    #print "2\n"
                    lastName = name
                    name = line.strip().strip(">")
                    tempDict[name] = ""
                else:
                    #print "3\n"
                    seq = line.strip()

                    tempDict[name] = tempDict[name] + seq
                                    
    #A subfunction of GetPercenIDFromClustalW.  Calculates percent identity
    #given two globally alligned amino acid sequence strings 
    def CalculatePerIDProt(self, first, second):
        zippy = list(zip(first, second))
        sameCnt = 0
        difCnt = 0
        #make the list of amino acids to deal with
        okayLst = "ACDEFGHIKLMNPQRSTVWY"
        
        #go through zippy and make a count of when the two are the same vs different
        for pair in zippy:
            #make sure the sequences are part of the understood amino acid list
            if not (pair[0] in okayLst and pair[1] in okayLst):
                if pair[0] == "-" or pair[1] == "-":
                    continue
                #handle the unknown amino acid, X
                elif pair[0] == "X" or pair[1] == "X":
                    difCnt += 1
                    continue
                else:
                    print("ERROR: UNKNOWN CHARACTER", pair)
                    sys.exit()

            #go through pairs and add to the counts
            if pair[0] == pair[1]:
                sameCnt += 1
            else:
                difCnt += 1

        #determine percent identity from counts
        total = sameCnt + difCnt
        perID = float(sameCnt) / float(total)
        
        return perID

    #A subfunction of GetPercenIDFromClustalW.  Calculates percent identity
    #given two globally alligned nucleotide sequence strings 
    def CalculatePerIDNuc(self, first, second):
        zippy = list(zip(first, second))
        sameCnt = 0
        difCnt = 0
        #make the list of amino acids to deal with
        okayLst = "ATCG"
        
        #go through zippy and make a count of when the two are the same vs different
        for pair in zippy:
            #make sure the sequences are part of the understood amino acid list
            if not (pair[0] in okayLst and pair[1] in okayLst):
                if pair[0] == "-" or pair[1] == "-":
                    continue
                #handle the unknown amino acid, X
                elif pair[0] == "N" or pair[1] == "N" or pair[0] == "M" or pair[1] == "M":
                    difCnt += 1
                    continue
                else:
                    print("ERROR: UNKNOWN CHARACTER", pair)
                    sys.exit()

            #go through pairs and add to the counts
            if pair[0] == pair[1]:
                sameCnt += 1
            else:
                difCnt += 1

        #determine percent identity from counts
        total = sameCnt + difCnt
        perID = float(sameCnt) / float(total)
        
        return perID

    # Created by David E. Hufnagel on Mar 25, 2013
    def SaveIntoDict(self, key, val, dictX):
        if key not in dictX:
            dictX[key] = [val]
        else:
            dictX[key].append(val)
                                
        

    def help(self):
        print(" -f function")
        print("    get_sequences - get specified sequences. REQUIRES: -fasta")
        print("       -name, Optional: -type, -match")
        print("    getseq2 - simpler function. REQUIRES: fasta,name. OPT:tokens")     
        print("    get_group_seq - get sequences for groups based on the group")
        print("       specification passed. REQUIRES: -fasta, -group. Optional:")
        print("       create_d, includ_o")
        print("    get_stretch - get segments of all sequences in fasta")
        print("       file, this is meant for alignment output. Need:")
        print("       -fasta, -coords. OPTIONAL: -format, -wcoord, -verbose")
        print("    get_stretch2 - get segment of a SINGLE sequence. Best for")
        print("       chromosome seq. REQUIRES: fasta, coords, OPT: isfile")
        print("    get_stretch3 - Take single seq fasta REQUIRES: fasta, coords")
        print("    get_stretch4 - Take multiple seq fasta, REQ: fasta,coords")
        print("    get_stretch5 - Take multiple seq fasta, REQ: fasta,coords")
        print("       included by Gaurav.")
        print("    fasta_to_oneline - convert fasta file to one line")
        print("       format. Requires: -fasta, coords. OPT: d")
        print("    oneline_to_fasta - convert the other way around")
        print("       REQUIRES: -oneline")
        print("    parse_desc - get the description out of fasta file")
        print("       NEED: fasta, OPT: style, D")
        print("    simplify - simplify the header by taking the id imme-")
        print("       diately after '>'.")
        print("       Requires: -fasta, Optional:-style, -D, tokens")
        print("    prefix - add whatever prefix to all sequences in")
        print("       a fasta file. Requires -fasta, -prefix")
        print("    suffix - add whatever suffix to all sequences in")
        print("       a fasta file. Requires -fasta, -suffix")      
        print("    size_filter - get sequences longer than the treshold")
        print("       passed. Requires: -fasta, -T. Optional: out_name")
        print("    compare_lists - compare two lists, requires: -files,")
        print("       optional: out")
        print("    cleanup - get rid of bracketed characters and replace")
        print("       with an '-'. Requires: -fasta")
        print("    get_sizes - get the sizes of sequences within a fasta")
        print("       file. REQUIRES: -fasta. OPTIONAL: x")
        print("    index_names - convert names into indices. Generate a")
        print("       fasta file with seq_id converted to index, and a")
        print("       '.name' file with both seq_id and indices")
        print("       REQUIRES: -fasta")
        print("    change_names - change names in the Fasta header.")
        print("       REQUIRES: -fasta, -name, -desc OPTIONAL:-delim")
        print("    rename - Change id within the fasta file")
        print("       REQUIRES: -fasta, -name")
        print("    del_redun - delete redundant sequences from Fastafile by name")
        print("       REQUIRES: -fasta")
        print("    del_redun_by_seq - delete redundant sequences from Fastafile")
        print("       by sequence. REQUIRES: -fasta")
        print("    check_redun - check fasta in a folder. NEED: targetd")
        print("    get_sp - get species for a particular species")
        print("       REQUIRES: -fasta, -sp")
        print("    parse_ensembl_fasta - parse desc, get len, one line")
        print("       seqeunce. REQUIRES: -fasta (in ensembl format)")
        print("    covert_header - coord from id<space>L-R to id_((L-1)*3+1)_R,")
        print("       NEED: fasta")
        print("    delete - rid of seqeunces, REQUIRES: -fasta,-dlist")
        print("    divide - divide fasta file into parts. NEED: fasta,by,setdir")
        print("    concat - concatenate sequences in a file, REQUIRES: -fasta")
        print("    locate - a companion function of concat. REQUIRES: c1, c2")
        print("    indiv  - Store each fasta squence in a file into its own")
        print("        REQUIRES: fasta")
        print("    mask   - masking the areas based on the passed coord.")
        print("        REQUIRES: fasta, coords")
        print("    count  - count the number of sequences, REQUIRES: fasta")
        print("    cnt_unique - count the number of unique sequences, REQUIRES: fasta")
        print("    get_longest - get the longest among alternative predictions")
        print("        NEED: fasta, D")
        print("    condense_seq - condenses all sequence lines into one.")
        print("        Need: fasta, out")
        print("    remove_by_keyword - removes all sequence lines containing")
        print("        the designated keyword. NEED: fasta, keyword, out")
        print("    replace_names - Replaces a code name in a fasta file with the real")
        print("        name (or vice versa). NEED: fasta, name, out")
        print("    replace_namesR - like replace_names, but takes the second column")
        print("        as the key for the ref dict instead of the first. NEED: fasta, name, out")
        print("    n50 - calculates and prints n50 NEED: fasta")
        print("    getPercentID - take a clustalw output fasta file and determine %ID between all")
        print("          pairwise clusters, outputting in the format: gene1  gene2  %ID")
        print("          WARNING: CLUSTERS MUST BE IN GROUPS OF 2. REQ: fasta, out, type")
        print(" -fasta : location and name of fasta file for getting or con-")
        print("   verting sequences.")
        print(" -group : group specification [group_name][seq_id]")
        print(" -oneline: sequence in the format [name][seq]")
        print(" -files : names separated by ','") 
        print(" -name  : file with lines of gene names. For change_names(), the")
        print("   file is in the format [new_name][old_name]")
        print(" -style : The fasta description file format:")
        print("     GenBank - gb")
        print(" -prefix: the prefix to be added to the sequence descriptors") 
        print(" -out   : output file name")
        print(" -T     : threshold length of sequences included in the output")
        print(" -start : the starting index, default = 1")
        print(" -type  : peptide [0,default], nucleotide [1]")
        print(" -D     : delimiter")
        print(" -species: two CAPITAL characters species abbreviation")
        print(" -dlist : a list of names to be deleted from fasta file")
        print(" -delim : delimiter between id and other info")
        print(" -match : 0: exacth match to fasta descriptor (default) or 1:")
        print("   exact match to the passed name right after '>'")
        print(" -desc  : whether the name contain description [1, default] or")
        print("          not [0]")
        print(" -coords: the coordinates of the segments to get, the number of")
        print("   coordinates should be even and none of them should be")
        print("   smaller than 1. They are separated by ','. Alternatively,")
        print("   a file with [id][L][R] can be specified (for mask).")
        print("   For get_stretch2: file contains two coord separated by ',' ")
        print("       get_stretch3: whatever...[L][R]")
        print("       get_stretch4: [gene][chr][L][R]")
        print("   in each line")
        print(" -wcoord: output sequence name with coords")
        print(" -create_d: create dir for each sequence groups, default no [0]")
        print(" -includ_o: include outgroup sequence [1], or not [0, default]")
        print("          For this to work, the fasta file passed should have a")
        print("          sequence called OUT.")
        print(" -ignore: For rename, ignore the string after the specified char")
        print("          during name matching. And add the string back later")
        print(" -format: 0 (default) for pep or nt coords for pep or nt seq file")
        print("          respectively. 1 for pep coords to get nt seq.")
        print(" -verbose: count line for fasta_to_dict [1], or not [0, default")
        print(" -by    : divide fasta file by the passed integer")
        print(" -add   : add domain info in get_strech. The last token in the")
        print("          coord file should be domain name. Default 0")
        print(" -x     : count masked region size")
        print(" -tokens: indice seperated by ',' defining the tokens to be")
        print("          taken as names in the passed name file. DEFAULT ''")
        print(" -setdir: make dir for each divided file [1], default no [0]")
        print(" -targetd : target dir")
        print(" -d     : delimiter for fasta_to_oneline, defailt ''")
        print(" -keyword: the keyword to use for parsing")
        print("")
        sys.exit(0)

        
#
# Function call
#

f_util = FileUtility.file_util()

if __name__ == '__main__':
    fasta = name_file = function = style = prefix = suffix=out_name = files = \
            oneline   = delimiter= species = dlist = delim = coords = group = \
            ignore    = by = c1 = c2 = tokens = targetD = d = keyword = "" 
    threshold = type = match = create_d = includ_o = format = verbose = \
            wcoord = segment = add =  call = x = setdir = seqid = 0
    start     = 1
    desc      = 1
    isfile    = 1
    manager   = fasta_manager()
    trans     = Translation.translate()

    for i in range(1,len(sys.argv)):
        if sys.argv[i] == "-fasta":
            fasta      = sys.argv[i+1]
        if sys.argv[i] == "-oneline":
            oneline    = sys.argv[i+1]
        if sys.argv[i] == "-name":        
            name_file  = sys.argv[i+1]
        if sys.argv[i] == "-f":  
            function   = sys.argv[i+1]
        if sys.argv[i] == "-style":      
            style      = sys.argv[i+1]
        elif sys.argv[i] == "-seqid":        
            seqid      = int(sys.argv[i+1])
        if sys.argv[i] == "-prefix":        
            prefix     = sys.argv[i+1]
        if sys.argv[i] == "-suffix":        
            suffix     = sys.argv[i+1]
        if sys.argv[i] == "-out":          
            out_name   = sys.argv[i+1]
        if sys.argv[i] == "-files":      
            files      = sys.argv[i+1]        
        if sys.argv[i] == "-D":  
            delimiter  = sys.argv[i+1]        
        if sys.argv[i] == "-T":  
            threshold  = int(sys.argv[i+1])
        if sys.argv[i] == "-start":      
            start      = int(sys.argv[i+1])
        if sys.argv[i] == "-segment":      
            segment    = int(sys.argv[i+1])
        if sys.argv[i] == "-type":        
            type       = int(sys.argv[i+1])
        if sys.argv[i] == "-match":      
            match      = int(sys.argv[i+1])
        if sys.argv[i] == "-species":      
            species    = sys.argv[i+1]
        if sys.argv[i] == "-dlist":      
            dlist      = sys.argv[i+1]
        if sys.argv[i] == "-delim":      
            delim      = sys.argv[i+1]
        if sys.argv[i] == "-desc":        
            desc       = int(sys.argv[i+1])
        if sys.argv[i] == "-coords":        
            coords     = sys.argv[i+1]
        if sys.argv[i] == "-group":     
            group      = sys.argv[i+1]
        if sys.argv[i] == "-create_d":       
            create_d   = int(sys.argv[i+1])
        if sys.argv[i] == "-includ_o":       
            includ_o   = int(sys.argv[i+1])
        if sys.argv[i] == "-ignore":        
            ignore     = sys.argv[i+1]
        if sys.argv[i] == "-format":        
            format     = int(sys.argv[i+1])         
        if sys.argv[i] == "-verbose":       
            verbose    = int(sys.argv[i+1])         
        if sys.argv[i] == "-by":        
            by         = int(sys.argv[i+1])         
        if sys.argv[i] == "-c1":        
            c1         = sys.argv[i+1]
        if sys.argv[i] == "-c2":        
            c2         = sys.argv[i+1]
        if sys.argv[i] == "-wcoord":        
            wcoord     = int(sys.argv[i+1])
        if sys.argv[i] == "-add":       
            add        = int(sys.argv[i+1])
        if sys.argv[i] == "-call":      
            call        = int(sys.argv[i+1])
        if sys.argv[i] == "-isfile":        
            isfile      = int(sys.argv[i+1])
        if sys.argv[i] == "-x":     
            x           = int(sys.argv[i+1])
        if sys.argv[i] == "-tokens":        
            tokens      = sys.argv[i+1]
        if sys.argv[i] == "-targetd":
            targetD     = sys.argv[i+1]
        if sys.argv[i] == "-d":
            d     = sys.argv[i+1]
        if sys.argv[i] == "-keyword":
            keyword     = sys.argv[i+1]

    if function == "get_sequences":
        if fasta == "" or name_file == "":
            print("\nNeed fasta file and name file\n")
            manager.help()
        manager.get_sequences(fasta,name_file,segment,type,match)
        
    elif function == "getseq2":
        if fasta == "" or name_file == "":
            print("\nNeed fasta and name.\n")
            manager.help()
        manager.getseq2(fasta,name_file,tokens)

    elif function == "get_group_seq":
        if fasta == "" or group == "":
            print("\nNeed fasta file and group specification.\n")
            manager.help()
        manager.get_group_seq(fasta,group,create_d,includ_o)

    elif function == "get_stretch":
        if fasta == "" or coords == "":
            print("\nNeed fasta file and coordinates\n")
            manager.help()
        manager.get_stretch(fasta,coords,format,wcoord,verbose,call,add)

    elif function == "get_stretch2":
        if fasta == "" or coords == "":
            print("\nNeed fasta file and coordinates\n")
            manager.help()
        manager.get_stretch2(fasta,coords,call,isfile)

    elif function == "get_stretch3":
        if fasta == "" or coords == "":
            print("\nNeed fasta file and coordinates\n")
            manager.help()
        manager.get_stretch3(fasta,coords)
    elif function == "get_stretch4":
        if fasta == "" or coords == "":
            print("\nNeed fasta file and coordinates\n")
            manager.help()
        manager.get_stretch4(fasta,coords,seqid)
    elif function == "get_stretch5":
        if fasta == "" or coords == "":
            print("\nNeed fasta file and coordinates\n")
            manager.help()
        manager.get_stretch5(fasta,coords)

    elif function == "fasta_to_oneline":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.fasta_to_oneline(fasta,d)
    elif function == "oneline_to_fasta":
        if oneline == "":
            print("\nNeed one line file\n")
            manager.help()
        manager.oneline_to_fasta(oneline)
    elif function == "simplify":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.simplify_desc(fasta,style,delimiter,tokens)

    elif function == "prefix":
        if fasta == "" or prefix == "":
            print("\nNeed fasta file and prefix\n")
            manager.help()
        manager.add_prefix(fasta,prefix)
    elif function == "suffix":
        if fasta == "" or suffix == "":
            print("\nNeed fasta file and suffix\n")
            manager.help()
        manager.add_suffix(fasta,suffix)        
    elif function == "size_filter":
        if fasta == "" and treshold == 0:
            print("\nNeed fasta files and length threshold\n")
            manager.help()
        manager.size_filter(fasta,threshold,out_name)
    elif function == "parse_desc":
        if fasta == "":
            print("\nNeed fasta file and description style info\n")
            manager.help()
        manager.parse_desc(fasta,style,delimiter)
    elif function == "compare_lists":
        if files == "" or len(files.split(",")) != 2:
            print("\nNeed files and it should be exactly 2 names\n")
            manager.help()
        manager.compare_lists(files,out_name)
    elif function == "cleanup":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.cleanup(fasta)  
    elif function == "get_sizes":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.get_sizes(fasta,0,x)
    elif function == "index_names":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.index_names(fasta,start)
    elif function == "change_names":
        if fasta == "" or name_file == "" or desc == -1:
            print("\nNeed fasta file and name file. Also need to specify desc\n")
            manager.help()
        manager.change_names(fasta,name_file,desc,delim)
    elif function == "rename":
        if fasta == "" or name_file == "":
            print("\nNeed fasta and name files\n")
            manager.help()
        manager.rename2(fasta,name_file,ignore) 
    elif function == "del_redun":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.del_redun(fasta)
    elif function == "del_redun_by_seq":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.del_redun_by_seq(fasta)
    elif function == "check_redun":
        if targetD == "":
            print("\nNeed dir\n")
            manager.help()
        manager.check_redun(targetD)
    elif function == "convert_header":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.convert_header(fasta)
    elif function == "get_sp":
        if fasta == "" or species == "":
            print("\nNeed fasta file and a species name specified\n")
            manager.help()
        manager.get_sp(fasta,species)
    elif function == "parse_ensembl_fasta":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.parse_ensembl_fasta(fasta,species)
    elif function == "delete":
        if fasta == "" or dlist == "":
            print("\nNeed fasta file and list of names\n")
            manager.help()
        manager.delete(fasta,dlist)
    elif function == "divide":
        if fasta == "" or by == "":
            print("\nNeed fasta file and the number of divisions\n")
            manager.help()
        manager.divide(fasta,by,setdir)
    elif function == "concat":
        if fasta == "":
            print("\nNeed fasta file\n")
            manager.help()
        manager.concat(fasta)
    elif function == "locate":
        if c1 == "" or c2 == "":
            print("\nNeed coordinate files\n")
            manager.help()
        manager.concat_locate(c1,c2)
    elif function == "indiv":
        if fasta == "":
            print("\nNeed fasta\n")
            manager.help()
        manager.indiv(fasta)
    elif function == "mask":
        if fasta == "" or coords == "":
            print("\nNeed fasta and coords\n")
            manager.help()
        manager.mask(fasta,coords)
    elif function == "fasta_to_dict":
        manager.fasta_to_dict(fasta)
    elif function == "count":
        if fasta == "":
            print("\nNeed fasta\n")
            manager.help()
        manager.count(fasta)
    elif function == "cnt_unique":
        if fasta == "":
            print("\nNeed fasta\n")
            manager.help()
        manager.cntUnique(fasta)
    elif function == "ensembl_longest":
        if fasta == "":
            print("\nNeed fasta\n")
            manager.help()
        manager.ensembl_longest(fasta)
    elif function == "get_longest":
        if "" in [fasta,delimiter]:
            print("\nNeed fasta and deliminter\n")
            manager.help()
        manager.get_longest(fasta,delimiter)
    elif function == "condense_seq":
        if fasta == "" or out_name == "":
            print("\nNeed fasta and/or out\n")
            manager.help()
        manager.CondenseSeq()
    elif function == "remove_by_keyword":
        if fasta == "" or keyword == "" or out_name == "":
            print("\nNeed fasta and/or keyword\n")
            manager.help()
        manager.RemoveByKeyword()
    elif function == "replace_names":
        if fasta == "" or name_file == "" or out_name == "":
            print("\nNeed fasta and/or keyword\n")
            manager.help()
        manager.ReplaceNames()
    elif function == "replace_namesR":
        if fasta == "" or name_file == "" or out_name == "":
            print("\nNeed fasta and/or keyword\n")
            manager.help()
        manager.ReplaceNames(reverse=True)
    elif function == "n50":
        if fasta == "":
            print("\nNeed fasta\n")
            manager.help()
        manager.N50()
    elif function == "getPercentID":
        if fasta == "" or out_name == "":
            print("\nNeed fasta and/or out_name\n")
            manager.help()
        manager.GetPercentIDFromClustalW()
    else:
        print("\nUnknown function...\n")
        manager.help()

