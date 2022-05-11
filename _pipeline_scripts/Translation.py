

import FastaManager,sys,os,ParseBlast,FileUtility,BlastUtility,string,time
from string import *
from bisect import bisect

class translate:

    ##
    # This is not for general things but for peptides derived from EST contigs.
    #
    # If searching for a nt_code and it is not found after a gap larger than 2
    # and the next 3 nt doesn't match the next aa, it will return an error
    # message for that particular sequence.
    #
    # This method ABSOLUTELY dependent that the correct corresponding start of
    # the sequence is given. It is not good at all in finding the start then
    # get the rest.
    #
    # @param  pep  peptide sequence, this should be derived from the nt seq
    #              either through direct translation or through extracted cds
    #              The seq comes with orientation information and coordinates.
    #              using tblastn. Either identical or very similar.
    # @param  nt   nucleotide sequence.
    # @param  exclstop [1] or not [0, default]
    #
    ##
    def bt(self,pep,cds,exclstop):

        # store sequences into dict
        manager = FastaManager.fasta_manager()

        print("Load pep file...")
        pdict = manager.fasta_to_dict(pep,0)

        print("Load cds file...")
        ndict = manager.fasta_to_dict(cds,0)
        code  = self.get_aa_code()

        oup = open(pep+".cds.fa","w")

        print("Total %i pep sequences:" % len(list(pdict.keys())))
        countS = 0
        for i in list(pdict.keys()):     # each peptide
            if countS % 1e2 == 0:
                print(" %i x 100" % (countS/1e2))
            countS = countS+1

            id = i
            aa = pdict[i]     # peptide seq
            nt = ""     # get the corresponding nucleotide seq
            L  = 0
            R  = 0
            
            # parse coordinates
            if id.find("_") == -1:
                print("No coord delimiter found, skip:",id)
                continue

            id = i.split("_")
            L  = int(id[-2])
            R  = int(id[-1])
            # reassemble ID so those with "_" in ID can be conserved
            id = i[:i.find("_%s_%s" %(id[-2],id[-1]))]
            #print "",id,L,R
            # extract nt sequence and make sure it is the right orientation
            if id in ndict:
                if L<R:
                    nt = ndict[id][L-1:R]
                else:
                    nt = ndict[id][R-1:L]
                    nt = self.complement(nt)
                    nt = self.reverse2(nt)
            else:
                print("NOT in ndict, skip:",id)
                continue

            #print ">",aa
            #print ">",nt
            if exclstop:
                seq,firstMis = self.matching(code,aa,nt)
            else:
                seq,firstMis = self.back_translate(code,aa,nt)
                
            if seq != "":
                oup.write(">%s\n%s\n" % (i,seq))
            else:
                print(" %s ->no sequence..." % i)
                    
            if firstMis:
                print(" %s ->first AA mismatch" % i)
    
    
    def matching(self,code,aa,nt):
        
        firstMis = 0
        # check if the first aa is at the beginning, skip the first if it is 
        # mismatched. Assuming from the second one, it will be ok...
        if not nt[:3] in code[aa[0]]:
            firstMis = 1
            P = 1
        else:
            P = 0
                    
        # read aa seq and match nt, generate synthetic cds
        seq = ""
        N = 0
        for i in range(P,len(aa)):
            if aa[i] not in code:
                print("Unknown code:",aa[i])
                if aa[i] == "B":
                    N += 3
                else: # such as X
                    pass
            else:
                for j in range(N,len(nt)):
                    if nt[j:j+3] in code[aa[i]]:
                            seq += nt[j:j+3]
                            N = j+3
                            break
                    elif j+4<len(nt) and nt[j+1:j+4] in code[aa[i]]:
                            seq += nt[j+1:j+4]
                            N = j+4
                            break
                    elif j+5<len(nt) and nt[j+2:j+5] in code[aa[i]]:
                            seq += nt[j+2:j+5]
                            N = j+5
                            break

            if seq[-3:] in code["*"]:
                seq = seq[:-3]
    
        return seq,firstMis 

    ##
    # This is part of the back translate function. Just make the reiteration
    # part into an independent method
    ##
    def back_translate(self,code,aa,nt):
        
        #print aa
        #print nt
        firstMis = 0
        # check if the first aa is at the beginning
        if not nt[:3] in code[aa[0]]:
            firstMis = 1
            print("",aa[0])
            print("",nt[:3])
            
        cn   = 0  # count nt

        # skip the first if it is mismatched. Assuming from the second one, it
        # will be ok...
        if firstMis: 
            ca = 1
        else:
            ca = 0  
        g    = 0  # count gap
        seq  = "" 
        last = "" # the previous nt code
        reset= 0
        # each aa in that peptide
        while ca < len(aa):

            if not aa[ca] in list(code.keys()):
                ca = ca+1
                continue
            #print "",aa[ca]
            
            if cn+3 > len(nt):
                #print nt[cn:]
                #print "End of nt sequence!"
                break

            if nt[cn:cn+3] in code[aa[ca]]:
                if last != "":
                    #print " addL>",last,"gap:%i-%i" % (cn-3-g,cn-3-1)
                    seq  = seq + last
                    last = ""
                    g    = 0
                #print " add >",nt[cn:cn+3],aa[ca]
                seq = seq + nt[cn:cn+3]
                cn = cn+3
            else:
                # give the previous one up,evaluate the current one again
                if last != "":
                    #print " omit > %s, no match" % aa[ca]
                    # reset
                    last = ""
                    cn   = cn-g
                    g    = 0
                    continue
                    
                # find recursively
                else:
                    #print " recursive find:",aa[ca]
                    found = 0
                    #print " cn :",cn
                    while cn < len(nt):
                        if nt[cn:cn+3] in code[aa[ca]]:
                            last = nt[cn:cn+3]
                            #print " last>",nt[cn:cn+3]
                            cn = cn+3
                            found = 1
                            break
                            
                        cn = cn+1
                        g  = g +1
                    #print " cn,gap:",cn,g
                    
                    if not found:
                        #print " omit> %ith aa, not found" % (ca+1)
                        cn  = cn-g
                        g   = 0
                    # if found and it is the last, write it
                    elif ca+1 == len(aa):
                        pass
                        #print " add >",nt[cn:cn+3],aa[ca]
            ca = ca+1

        return seq,firstMis      

    ##
    # A back translation utility based on synchronized pep and nt sequences
    #
    # 4/14,03, problems when sequence is in lower case. Taken care of.
    # 3/11,03, implement verbose
    ##
    def back_translate2(self,pep,nt,v=0):
        
        if v:
            print("\nBack translate:")
        # verify length first
        plen = 0
        if pep.find("-") != -1:
            plist = pep.split("-")
            for i in plist:
                plen += len(i)
        else:
            plen = len(pep)
            
        # nt sequence don't have "-", but need to rid of stop
        if nt[-3:] in ["TAA","TGA","TAG","taa","tga","tag"]:
            nt = nt[:-3]
        
        flag = 0
        if len(nt) != plen*3:
            if v:
                print("Size discrepancy: nt",len(nt),"pep",plen)
            flag = 1
            
        seq = ""
        if not flag:
            countN = 0
            codes = self.get_nt_code()
            for i in pep:
                if i == "-":
                    seq += "---"
                else:
                    seq += nt[countN:countN+3]
                    
                    # checking if it translate alright
                    if nt[countN:countN+3].upper() not in codes:
                        if v:
                            print("Code unknown:",i.upper(),\
                                  nt[countN:countN+3].upper())
                    elif codes[nt[countN:countN+3].upper()] != i.upper():
                        if v:
                            print("Discrepancy:",i.upper(),\
                                  nt[countN:countN+3].upper())                      
                    countN += 3
            
        return seq, flag
    

    
    #
    # @param id  for translating one particular sequence in a Fasta file.
    #            default is an empty string means translate the whole file.
    # 
    def translate(self,cds,id="",frame=0,discard=0,indiv=0,lenT=0):
        
        print("Read fasta into dict...")
        cdict = manager.fasta_to_oneline(cds,1)
        
        if id != "":
            print("Translate:", id)
            oup = open("%s_F%i.trans" % (id,frame),"w")
            oup.write(">%s\n" % id)
            oup.write("%s\n" % self.translate_passed(cdict[id].upper()))
        else:
            oup = open(cds+"_F%i.trans"     % frame,"w")
            oup1= open(cds+"_F%i.trans.log" % frame,"w")
            ckeys = list(cdict.keys())
            ckeys.sort()
            countT = 0
            count1 = 0
            count2 = 0
            count3 = 0
            countW = 0
            udict = {}
            print("Translate each seq...")
            for i in ckeys:
                if frame == 1:
                    cdict[i] = "NN" + cdict[i]
                elif frame == 2:
                    cdict[i] = "N" + cdict[i]
                [pep,error,unk_dict] = self.translate_passed2(cdict[i].upper())
                write = 0
                if error == 1:
                    print(" ERR-length:",i)
                    count1 += 1
                elif error == 2:
                    print(" ERR-unk codons:",i)
                    count2 += 1
                else:
                    write = 1
                
                # check if there is stop that is not the last residue
                if pep[:-1].find("*") != -1 and pep[:-1][-1] != "*":
                    print(" ERR-inner stop:",i)
                    #print pep
                    #sys.exit(0)
                    count3 += 1
                    error = 3
                    write = 0
                
                for j in unk_dict:
                    if j in udict:
                        udict[j] += 1
                    else:
                        udict[j] = 1
                
                if discard and not write:
                    #print " here"
                    pass
                else:
                    countW += 1
                    #print " there"
                    oup.write(">%s\n" % i)
                    oup.write("%s\n" % pep)
                    oup1.write("%s\terror%i\n" % (i,error))
                countT += 1
            
            print("Total %i sequences translated."      % countT)
            print("  %i has length not divisable by 3." % count1)
            print("  %i has inner stop codon(s)."       % count3)
            print("  %i has unknown codons."            % count2)
            print("Discard flag: %i"                    % discard)
            print("  %i in output file"                 % countW)
                        
            if count2 > 0:
                print(" listed:")
                for i in udict:
                    print("",i,udict[i])

    #
    # Break ORFs into subORFs that >= size threshold. All ORF should start with
    # M.
    #
    def suborf(self,orfcds,lenT,check):
        
        print("Read sequence file...")
        fdict = manager.fasta_to_dict(orfcds,0)
        print(" total %i\n" % len(list(fdict.keys())))
        
        print("Iterate throughs sequences...")
        oup = open(orfcds+".sub","w")
        countS = 0
        countT = 0
        for i in fdict:
            if countT % 1e4 == 0:
                print(" %i x 10k" % (countT/1e4))
            countT += 1
            #print i
            # chr|ori|frame|L-R
            sid     = i.split("|")
            sid[-1] = sid[-1].split("-")
            
            # there are situations where there is negative coord like:
            # Contig989|-|2|-2-192
            # in all cases, they are -2. And should be corrected by +3
            
            
            if check and len(sid[-1]) == 3:
                if sid[-1][0] == "":
                    sid[-1] = ["1",sid[-1][2]]
                else:
                    print("ERR1:",[i])
                    continue
                        
            L       = int(sid[-1][0])
            R       = int(sid[-1][1])
            cds     = fdict[i]
            
            # check size
            if check and abs(R-L)+1 != len(cds):
                if abs(R-L)+1-3 == len(cds):
                    R = R-3
                else:
                    print("ERR2:",[i],abs(R-L)+1,len(cds))
                    continue
            
            # construct new seq names
            seqname = "%s|%s|%s|%i-%i" % (sid[0],sid[1],sid[2],L,R)
                            
            pep    = self.translate_passed2(cds)[0]
            
            m = 0
            mprime = m 
            l = len(pep)
            #print "> ",m,l,lenT
            #print [pep]
            #print [cds]
            while l > lenT:
                countS += 1
                oup.write(">%s|%i\n%s\n" % (seqname,m*3,cds))
                #print ">>write"
                mprime = pep[1:].find("M")
                m += mprime + 1
                pep = pep[1:][mprime:]
                cds = cds[3:][mprime*3:]
                l = len(pep)
                #print ">>",m,l
                #print [pep]
                #print [cds]
        
        print("%i sequences, %i suborfs" % (countT,countS)) 
        print("Done!\n")            
    
    
    
    
    ###
    # 6 frame translation. Not for getting ORFs, just straight translation.
    #
    # @param seq    fasta sequence file
    ###
    def sixpack_simple(self,seq):
        print("Read fasta into dict:")
        fdict = manager.fasta_to_dict(seq,0)
        
        oup1 = open("%s.6pack_cds" % seq,"w")
        oup2 = open("%s.6pack_pep" % seq,"w")
        for i in fdict: 
            f   = fdict[i]          # forward
            fn0 = f                 # forward frame 0
            fn1 = f[1:]
            fn2 = f[2:]
            # Rid of extra 1 or 2 nt
            if len(fn0) % 3 != 0:
                fn0 = fn0[:-(len(fn0)%3)]
            if len(fn1) % 3 != 0:
                fn1 = fn1[:-(len(fn1)%3)]
            if len(fn2) % 3 != 0:
                fn2 = fn2[:-(len(fn2)%3)]
                
            r   = self.rc(f)        # reverse
            rn0 = f                 # reverse frame 0
            rn1 = f[1:]
            rn2 = f[2:]
            if len(rn0) % 3 != 0:
                rn0 = rn0[:-(len(rn0)%3)]
            if len(rn1) % 3 != 0:
                rn1 = rn1[:-(len(rn1)%3)]
            if len(rn2) % 3 != 0:
                rn2 = rn2[:-(len(rn2)%3)]
            
            # Get translation
            fp0 = self.translate_passed2(fn0,0)[0]
            fp1 = self.translate_passed2(fn1,0)[0]
            fp2 = self.translate_passed2(fn2,0)[0]
            rp0 = self.translate_passed2(rn0,0)[0]
            rp1 = self.translate_passed2(rn1,0)[0]
            rp2 = self.translate_passed2(rn2,0)[0]
            
            oup1.write(">%s_f0\n%s\n>%s_f1\n%s\n>%s_f2\n%s\n" % (i,fn0,i,fn1,i,fn2) +\
                       ">%s_r0\n%s\n>%s_r1\n%s\n>%s_r2\n%s\n" % (i,rn0,i,rn1,i,rn2))
            oup2.write(">%s_f0\n%s\n>%s_f1\n%s\n>%s_f2\n%s\n" % (i,fp0,i,fp1,i,fp2) +\
                       ">%s_r0\n%s\n>%s_r1\n%s\n>%s_r2\n%s\n" % (i,rp0,i,rp1,i,rp2))
                       
        print("Done!")
    
    #
    # Batch call six pack
    #
    def batch_6pack(self,seq,lenT,lenT2,met,inc,verbose):
    
        oup0  = open("%s_T%i-%im%i.cds"   % (seq,lenT,lenT2,met),"w")
        oup1  = open("%s_T%i-%im%i.pep"   % (seq,lenT,lenT2,met),"w")
        oup2  = open("%s_T%i-%im%i.coord" % (seq,lenT,lenT2,met),"w")
        oup2.write("SeqID\tOri\tFrame\tL\tR\n")
        
        print("Read fasta into dict:")
        fdict = manager.fasta_to_dict(seq,0)
        c = 0
        
        print("Do 6pack:")
        for i in fdict:
            print("",i)
            eflag = self.sixpack([i,fdict[i]],lenT,lenT2,met,inc,
                                  oup0,oup1,oup2,verbose)
            if eflag:
                break
        
        print("Done!")

    #
    # sixpack: 6 frame translation
    #
    # PROBLEM: met set to 2 is not working for antisense strand.
    # 05/16,05 Stop codon position is included
    # 05/17,05 lenT2 is not applied to met = 0
    #
    # @param seq   a fasta file with 1 seq
    # @param lenT  leng cutoff for ORF, 25 default 
    # @param lenT2 leng threshold for ORF, 10000 default
    # @param met   methionine as ORF start [1,default] or any [0], or get 
    #              all ORFs start with Met in any full length ORF [2] NOT
    #              WORKING!!!!!
    # @param oup0  Passed from batch_6pack, for ORF cds
    # @param oup1  Passed from batch_6pack, for ORF pep
    # @param oup2  Passed from batch_6pack, for ORF coord
    #
    def sixpack(self,seq,lenT,lenT2,met,inc,oup0="",oup1="",oup2="",verbose=1):
        
        if verbose:
            if oup1 == "":
                print("\nSequence:",seq)
            else:
                print("\nSequence:",seq[0])
            print("lenT    :",lenT)
            print("lenT2   :",lenT2)
            print("Met_flag:",met)
            print("TL_inc  :",inc)
        
        if inc % 3 != 0:
            print("ERR: TL_inc has to be multiple of 3")
            if oup1 == "":
                sys.exit(0)
            else:
                return 1
        
        ###################
        # split() start
        # seq_index,ori,ntseq,TL_seq,frame,lenT,met,antisense,antisense length
        ###################
        def split(idx,ori,ntseq,s,f,a=0,alen=0,verbose=1):          
            
            # L coord as key, [R_coord,pep,nt] as value
            #print s
            sdict = {}
            s = s.split("*")
            #print s
            c = f
            countT = 0
            
            if verbose:
                print("  split_orf:",len(s))
            for i in s:
                if verbose and countT != 0 and countT % 10000 == 0:
                    print("   %i x 10k" % (countT/10000))
                countT += 1
                #print "ORF:",i
                if i == "":
                    c += 3
                    continue
                
                EXC = ["?","X"] 
                #EXC = ["X"]
                # if orf starts or ends in ? or X, this orf will be disgarded
                if len(i) >= lenT and i[0] not in EXC and i[-1] not in EXC:
                    iL = len(i)
                    if met == 1 and i.find("M") != -1:
                        m = i.find("M")
                        mprime = m
                        x = i[m:]
                        l = len(x)
                        #print m,x,l
                        while bisect([lenT,lenT2+1],l) == 2:
                            #print " .."
                            mprime = x[1:].find("M")
                            m += mprime + 1
                            x = x[1:][mprime:]
                            l = len(x)
                        
                        # len(x) is between lenT and lenT2
                        if bisect([lenT,lenT2],l) == 1:
                            if a:
                                ntL    = c+m*3+1           # ntseq coord
                                ntR    = c+(iL+1)*3 
                                coordL = alen-(c+iL*3)+1-3 # sense strand
                                coordR = alen-(c+m*3+1)+1
                            # sense
                            else:
                                ntL    = c+m*3+1           # nt seq coord
                                ntR    = c+(iL+1)*3 
                                coordL = c+m*3+1           # sense strand
                                coordR = c+(iL+1)*3                             
                            if len(x) > 0:
                                sdict[coordL] = [coordR,x,ntseq[ntL-1:ntR]]
                            #print "IN :",x
                            #print "    ",ntL,ntR,coordL,coordR

                        else:
                            #print "OUT:",x
                            pass

                    elif met == 0:
                        # NOTICE THAT lenT2 is not applied.
                        # antisense
                        if a:
                            sdict[alen-(c+iL*3)+1-3] = [alen-(c+1)+1,i,
                                                        ntseq[c+1-1:c+(iL+1)*3]]
                        # sense
                        else:
                            sdict[c+1] = [c+(iL+1)*3,i,ntseq[c+1-1:c+(iL+1)*3]]         

                    """ PROBLEM HERE, INACTIVETED FOR NOW
                    elif met == 2:
                        m = i.find("M") # initial position for MET
                        n = m           # variable aa position for Met
                        x = i
                        while m != -1:                          
                            x = x[m:]
                            if len(x) >= lenT:
                                # antisense
                                if a:
                                    coordL = alen-(c+iL*3)+1
                                    coordR = alen-(c+n*3+1)+1
                                    # deal with ? and X
                                    if x.find("?") != -1:
                                        ...                                 
                                
                                # sense
                                else:
                                    coordL = c+n*3+1
                                    coordR = c+iL*3
                                
                                if len(x) > 0:
                                    sdict[coordL]  = [coordR,x]
                            else:
                                break
                            x = x[1:]
                            m = x.find("M")
                            n += m+1
                    """
                            
                c += len(i)*3 + 3           
            
            skeys = list(sdict.keys())
            skeys.sort()
            if verbose:
                print("  qualified:",len(sdict))
            for i in skeys:
                #print idx,ori,f,i,sdict[i][0]
                #print sdict[i][1]
                #print sdict[i][2]
                
                # Get rid of sequences that have ? or X
                if sdict[i][1].find("?") != -1 or sdict[i][1].find("X") != -1:
                    continue
                
                # 06/-6,07 for ORFs involving the 1st or the last nucleotide, 
                # the coordinates are -3 or +3 more than they should. This is
                # not so much of a problem when working with chromosomes, but
                # an issue when looking to multiple contigs. The following
                # is to check if there is discrepancy between the length of sub
                # seq and if they are the beginning or ending entries. If so,
                # coordinates are corrected accordingly.
                err = 0
                modi = 0
                if i < 0:
                    if i == -2 and sdict[i][0] == len(sdict[i][2]):
                        modi = 1
                    else:
                        " ERR:",idx,ori,f,i,sdict[i][0]
                #print sdict[i][0]-i+1, len(ntseq)
                if not modi:
                    if sdict[i][0]-i+1 > len(ntseq):
                        if sdict[i][0]-len(ntseq) == 2 and \
                        sdict[i][0]-i+1-3 == len(ntseq):
                            sdict[i][0] -= 3
                        else:
                            " ERR:",idx,ori,f,i,sdict[i][0]             
                else:
                    if sdict[i][0]-modi+1 > len(ntseq):
                        if sdict[i][0]-len(ntseq) == 2 and \
                        sdict[i][0]-modi+1-3 == len(ntseq):
                            sdict[i][0] -= 3
                        else:
                            " ERR:",idx,ori,f,i,sdict[i][0]             
                
                if not err:
                    oup0.write(">%s|%s|%i|%i-%i\n%s\n" % \
                            (idx,ori,f,i,sdict[i][0],sdict[i][2]))                      
                    oup1.write(">%s|%s|%i|%i-%i\n%s\n" % \
                            (idx,ori,f,i,sdict[i][0],sdict[i][1]))
                    oup2.write("%s\t%s\t%i\t%i\t%i\n"  % \
                            (idx,ori,f,i,sdict[i][0]))
        ###################
        # split() end
        ###################
                            
        if verbose:
            print("Read sense strand...")
        
        # output stream not specified, must be direct call
        if oup1 == "":
            inp   = open(seq,"r")
            inl   = inp.readlines()
            idx   = self.rmlb(inl[0])[1:]
            sense = string.joinfields(inl[1:],"")
            oup0  = open("%s_T%i-%im%i.cds"   % (idx,lenT,lenT2,met),"w")
            oup1  = open("%s_T%i-%im%i.pep"   % (idx,lenT,lenT2,met),"w")
            oup2  = open("%s_T%i-%im%i.coord" % (idx,lenT,lenT2,met),"w")
            oup2.write("SeqID\tOri\tFrame\tL\tR\n")
        # batch call, seq passed is a list with [idx,sequence]
        else:
            idx   = seq[0]
            sense = seq[1]

        # get rid of line breaks, sometimes \n and \r\n are mixed.
        if sense.find("\r\n") != -1:
            sense = string.joinfields(sense.split("\r\n"),"")
        if sense.find("\n") != -1:
            sense = string.joinfields(sense.split("\n"),"")
            
        if verbose:
            print("Reverse...")
        antis = self.reverse2(sense)        
        if verbose:
            print("Complement...")
        antis = self.complement(antis)
        
        # sense, frame 0,1,2
        if verbose:
            print("Translate..")
            print(" sense, frame 0")
        s0 = self.translate_passed(sense,0,inc)[0]
        split(idx,"+",sense,s0,0,0,0,verbose)

        if verbose:
            print(" sense, frame 1")
        s1 = self.translate_passed(sense,1,inc)[0]
        split(idx,"+",sense,s1,1,0,0,verbose)
        
        if verbose:
            print(" sense, frame 2")
        s2 = self.translate_passed(sense,2,inc)[0]
        split(idx,"+",sense,s2,2,0,0,verbose)
        # antisense, frame 0,1,2
        if verbose:
            print(" antisense, frame 0")
        a0 = self.translate_passed(antis,0,inc)[0]
        split(idx,"-",antis,a0,0,1,len(antis),verbose)
        
        if verbose:
            print(" antisense, frame 1")
        a1 = self.translate_passed(antis,1,inc)[0]
        split(idx,"-",antis,a1,1,1,len(antis),verbose)
        
        if verbose:
            print(" antisense, frame 2")
        a2 = self.translate_passed(antis,2,inc)[0]
        split(idx,"-",antis,a2,2,1,len(antis),verbose)
        
        #oup1.close()
        #oup2.close()

        if verbose:
            print("Done!")  

    def translate_passed(self,seq,frame=0,inc=1000):
        
        code = self.get_nt_code()
        
        if frame == 1:
            seq = seq[1:]
        elif frame == 2:
            seq = seq[2:]
        
        if len(seq) > 1000000:
            print("  nt_len: %i" % len(seq))
        
        inc = inc*1000
        
        prot_all = []
        for i in range(0,len(seq),inc):
            if len(seq) > 1000000:
                print("   %i-%i kb" % (i/1000,(i+inc)/1000))
                
            # 090911: Deal with fragments that are not multiple of 3.
            segment = seq[i:i+inc]
            #print "pre :",segment  
            blah = len(segment)%3
            segment = segment[:-blah]
            #print "post:",segment
            
            prot = ""
            for j in range(0,len(segment),3):
                prot += code.get(segment[j:j+3], "?")
            prot_all.append(prot)
        prot_all = string.joinfields(prot_all,"")
        
        """
        prot_all = ""
        for i in xrange(0,len(seq),3):
            prot_all += code.get(seq[i:i+3], "?")
        """
        
        return [prot_all,0,{}]          
    
    
    #
    # Exclude
    #
    def exclude(self,orf_coord,exc_coord):
    
        print("Read %s...\n" % exc_coord)
        inp = open(exc_coord)
        inl = inp.readline()
        # {idx:{L:R}}
        E = {}
        while inl != "":
            inl = inl.split("\t")
            i = inl[0]
            L = int(inl[1])
            R = int(inl[2])
            if i in E:
                if L in E[i]:
                    if R > E[i][L]:
                        E[i][L] = R
                else:
                    E[i][L] = R
            else:
                E[i] = {L:R}
            inl = inp.readline()
        
        print("Generated sorted exclusion coordinates...\n")
        Esorted = {}
        for i in E:
            ikeys = list(E[i].keys())
            ikeys.sort()
            Esorted[i] = ikeys
            
        print("Read %s and generate outputs..." % orf_coord)
        inp = open(orf_coord)
        oup = open(orf_coord+".qualified","w")
        inp.readline() # first line is header
        inl = inp.readline()
        countA = 0
        countQ = 0
        c = 0
        while inl != "":
            if countA % 1e4 == 0:
                print(" %i x 10k" % (countA/1e4))
            countA += 1
            S   = inl.split("\t")
            L   = int(S[3])
            R   = int(S[4])
            elist = Esorted[S[0]]
            ins = bisect(elist,L)

            #print [inl]
            #print "5':",elist[ins-1],E[S[0]][elist[ins-1]]
            #print "3':",elist[ins],E[S[0]][elist[ins]]     
            # evaluate 5', as long as it is not the first
            O5 = O3 = 1
            if ins != 0:
                if L > E[S[0]][elist[ins-1]]:
                    O5 = 0
            # evaluate 3', as long as it is not the last
            if ins != len(elist):
                if R < elist[ins]:
                    O3 = 0
                    
            if not O5 and not O3:
                oup.write(inl)
                countQ += 1
            
            #print "O5:%i,O3:%i" % (O5,O3)
            
            inl = inp.readline()
        
        print(" total %i, %i qualified\n" % (countA,countQ))
        
    
    #
    # Read a fasta file with one seq, return a list with [id,seq]
    #
    def read_1seq(self,fasta):
        inp = open(fasta,"r")
        inl = inp.readlines()
        idx = self.rmlb(inl[0])[1:]
        fasta = string.joinfields(inl[1:],"")
        
        # get rid of line breaks, sometimes \n and \r\n are mixed.
        if fasta.find("\r\n") != -1:
            sense = string.joinfields(fasta.split("\r\n"),"")
        if fasta.find("\n") != -1:
            fasta = string.joinfields(fasta.split("\n"),"")
        
        return([idx,fasta])
    
    #
    # Do reverse complement on sequences based on passed names
    # > whatever_xxxx|yyyy
    # 
    # if x > y, then it will be rc'ed.
    #
    def rc2(self,nt):
        fdict = manager.fasta_to_dict(nt,0)
        oup = open(nt+".rc_correct","w")
        
        fkeys = list(fdict.keys())
        countR = 0
        for i in fkeys:
            rc = 0 # rc or not flag
            if i.find("|") != -1 and i.find("_") != -1:
                # coordinates
                c = i[i.rfind("_")+1:].split("|")
                if int(c[0]) > int(c[1]):
                    rc = 1
                    countR += 1
            if rc:
                oup.write(">%s\n%s\n" % \
                    (i,self.complement(self.reverse2(fdict[i]))))
            else:
                oup.write(">%s\n%s\n" % (i,fdict[i]))
        
        print("%i sequences, %i rc'ed" % (len(fkeys), countR))

    #
    # Do reverse complement on sequences based on passed names. If no name
    # passed, all will be rc'd.
    #
    def batch_rc(self,nt,names):
        fdict = manager.fasta_to_dict(nt,0)
        oup = open(nt+".rc_select","w")
        
        if names != "":
            inp = open(names)
            inl = inp.readlines()
            names = {}
            for i in inl:
                names[self.rmlb(i)] = 0
                
        fkeys = list(fdict.keys())
        fkeys.sort()
        countR = 0
        for i in fkeys:
            if i in names:
                names[i] = 1
                oup.write(">%s\n%s\n" % \
                    (i,self.complement(self.reverse2(fdict[i]))))
                countR += 1
            else:
                oup.write(">%s\n%s\n" % (i,fdict[i]))
        
        for i in names:
            if names[i] == 0:
                print("Not found:",i)
        
        print("%i sequences, %i rc'ed" % (len(fkeys), countR))
    
    #
    # @param call internal[0,default] or outside[1].
    #
    def rc(self,nt,call=0):
        # external call, read file
        if call:
            fn  = nt
            nt  = self.read_1seq(nt)
            idx = nt[0]
            nt  = nt[1]
        
        nt = self.complement(self.reverse2(nt)) 
        
        if call:
            oup = open(fn+".rc","w")
            oup.write(">%s\n%s\n" % (idx,nt))   
        else:
            return nt

    def complement(self,nt):
		# 5/11/22: string.maketrans is deprecated, use str.maketrans.
        comp = nt.translate(str.maketrans("AGCTagct","TCGAtcga"))
        return comp

    def reverse2(self,a):
        b = list(a)
        b.reverse()
        a = "".join(b)
        # check if the reversed line has "\n" in the beginning, if so, add it 
        # to the end. If this is not done, some fasta file will have the last
        # sequence line fused with the next sequence id
        try:
            if a[0] == "\n":
                a = a[1:] + "\n"
        except IndexError:
            print("IndexErr:",a)
        return a

    
    ##
    # This is for direct call
    #
    # @param seq   nucleotide seqeunce
    # @param frame the frame to translate, default 0
    ##
    def translate_passed2(self,seq,frame=0):
        
        error = 0
        unk_codon = {}
        codes = self.get_nt_code()
        if frame < 0 or frame > 2:
            print("Frame should be 0-2, quit...")
            sys.exit(0)

        if len(seq)%3 != 0:
            error = 1
            
        c       = frame
        out_str = ""
        countAA = 1
        while c < len(seq):
            try:
                out_str = out_str + codes[seq[c:c+3].upper()]
            except KeyError:
                if seq[c:c+3] in unk_codon:
                    unk_codon[seq[c:c+3]] += 1
                else:
                    unk_codon[seq[c:c+3]] = 1
                if error != 1:
                    error = 2
                out_str += "X"
                #sys.exit(0)
            c = c+3
            countAA += 1

        return [out_str,error,unk_codon]
    
    #
    # Take the output of ParseBlast.parse_align4 to generate a pair file and
    # associated sequences that:
    # 1) has the same length
    # 2) has no stop, gap, or unknown codons
    #
    def tl_mindless(self,align_seq,lenT=75):
        inp = open(align_seq)
        inl = inp.readlines()
        codes = self.get_nt_code()
        o1 = open("%s.pairs" % align_seq,"w")
        o2 = open("%s.fake_cds" % align_seq,"w")
        Q = {}
        S = {}
        c = 0
        for i in range(0,len(inl),3):
            if c%1e4 == 0:
                print(" %i x10k" % (c/1e4))
            c += 1
            n = inl[i][1:].strip().split(" ")
            q = n[0] + "|" + n[2].split("|")[0]
            s = n[1] + "|" + n[2].split("|")[1]
            # Although the name may be redundant, because how the sequences
            # are being proecssed, the sequences may end up to be different.
            # So a count for when a q or s is encoutered is added.
            if q not in Q:
                Q[q] = 1
                q2 = "%s.1" % q
            else:
                Q[q] += 1
                q2 = "%s.%i" % (q,Q[q])
            if s not in S:
                S[s] = 1
                s2 = "%s.1" % s
            else:
                S[s] += 1
                s2 = "%s.%i" % (s,S[s])             
            
            qS= inl[i+1].strip().upper()
            sS= inl[i+2].strip().upper()
            
            # Get them so they are divisable by 3
            remainder = len(qS) % 3
            if remainder != 0:
                qS = qS[:-remainder]
                sS = sS[:-remainder]
            #print "len%3:",len(qS)%3,len(sS)%3
            #print qS
            #print sS
            qS2 = ""
            sS2 = ""
            for j in range(0,len(qS),3):
                qC = qS[j:j+3]
                sC = sS[j:j+3]
                #print "",qC,sC
                # Don't allow gaps
                if "-" in qC or "-" in sC:
                    #print "  gap"
                    continue
                elif qC in ["TAA","TGA","TAG"] or sC in ["TAA","TGA","TAG"]:
                    #print "  stop"
                    continue
                elif qC not in codes or sC not in codes:
                    continue
                qS2 += qC
                sS2 += sC
                #print ">>",qS2
                #print ">>",sS2
                
            if len(qS2) >= lenT:
                o1.write("%s\t%s\n" % (q2,s2))
                o2.write(">%s\n%s\n" % (q2,qS2))
                o2.write(">%s\n%s\n" % (s2,sS2))
                
        print("Done!")
                
    ##
    # Check the correspondance between cds and pep entries
    #
    # @param frame 0 (default), 1, or 2
    ##
    def validate(self,cds,pep,frame=0):

        # read files into dict
        
        cdict = manager.fasta_to_oneline(cds,1)
        pdict = manager.fasta_to_oneline(pep,1)

        countM = 0
        countS = 0
        for i in list(cdict.keys()):
            tlseq = self.translate_passed(cdict[i],frame)
            if pdict[i] == tlseq:
                countM += 1
            else:
                countS += 1
                #print " ",i
                #print tlseq
                #print pdict[i]
            
                seq1    = "tmp1.fa"
                seq2    = "tmp2.fa"
                bls_out = "tmp.out"
                outname = "validate.bl2"

                oup = open(seq1,"w")
                oup.write(">%s\n%s\n" % (i,tlseq))
                oup = open(seq2,"w")
                oup.write(">%s\n%s\n" % (i,pdict[i]))
                oup.close()   # THIS IS CRUCIAL!!

                os.system("bl2seq -p blastp -i %s -j %s -o %s -F F -e 1e-10" % \
                          (seq1,seq2,bls_out))
                os.system("cat tmp.out %s > tmp.combined" % outname)
                os.system("mv tmp.combined %s" % outname)

                # rid of the files
                os.system("rm tmp*")

        print("%i match, %i mismatch" % (countM,countS))

    #
    # This function takes the pair list, get the sequences for pairs, do a bl2seq
    # parse it, and output for verification.
    #
    # @param pairs  a file with pairs to be examined, in the format:
    #       [id1][id2][whatever]
    # @param fasta  file with the sequences
    #
    def bl2_pairs(self,fasta,pairs,blast_dir,prog):

        # clean up the list so id2 only match once, take the highest value on
        # the third column
        
        inp    = open(pairs,"r")
        inline = inp.readline()
        pdict  = {}
        while inline != "":
            llist = self.rmlb(inline).split("\t")
            if len(llist) == 3:
                score = float(llist[2])
                if llist[1] not in pdict:
                    pdict[llist[1]] = [llist[0],score]
                else:
                    if score > pdict[llist[1]][1]:
                        print(llist[1])
                        print(pdict[llist[1]])
                        print([llist[0],score])
                        pdict[llist[1]] = [llist[0],score]
                        
            else:
                if llist[1] not in pdict:
                    pdict[llist[1]] = [llist[0],1]
                else:
                    print("Redundant:",llist[1])
                                   
            inline = inp.readline()
        
        # read fasta 
        fasta = manager.fasta_to_dict(fasta,0)
        
        print("Conduct bl2seq:")
        oup = open(pairs+".mredun","w")
        for i in list(pdict.keys()):
            id2 = i
            id1 = pdict[i][0]

            print("%s,%s" % (id1,id2))
            oup_tmp = open("tmp1.fa","w")
            oup_tmp.write(">%s\n%s\n" % (id1,fasta[id1]))
            oup_tmp = open("tmp2.fa","w")
            oup_tmp.write(">%s\n%s\n" % (id2,fasta[id2]))
            oup_tmp.close()
            outname = pairs+".bl2"
            
            os.system("%s/bl2seq -p %s -i %s -j %s -o %s -F F -e 1e-10" % \
                      (blast_dir,prog,"tmp1.fa","tmp2.fa","tmp.out"))

            # blast,T,wself,format,style,query=""
            parser.parse_align("tmp.out",70,0,1,0,id1)

            os.system("cat tmp.out.log %s > tmp.combined" % outname)
            os.system("mv tmp.combined %s" % outname)

    
    # define univeral codes and return a dict with such
    def get_nt_code(self):
        # /: between codon frameshift, reported by e.g. tfasty
        # |: within codon frameshift
        code = {"TTT":"F","TCT":"S","TAT":"Y","TGT":"C",
                "TTC":"F","TCC":"S","TAC":"Y","TGC":"C",
                "TTA":"L","TCA":"S","TAA":"*","TGA":"*",
                "TTG":"L","TCG":"S","TAG":"*","TGG":"W",
                 
                "CTT":"L","CCT":"P","CAT":"H","CGT":"R",
                "CTC":"L","CCC":"P","CAC":"H","CGC":"R",
                "CTA":"L","CCA":"P","CAA":"Q","CGA":"R",
                "CTG":"L","CCG":"P","CAG":"Q","CGG":"R",

                "ATT":"I","ACT":"T","AAT":"N","AGT":"S",
                "ATC":"I","ACC":"T","AAC":"N","AGC":"S",
                "ATA":"I","ACA":"T","AAA":"K","AGA":"R",
                "ATG":"M","ACG":"T","AAG":"K","AGG":"R",

                "GTT":"V","GCT":"A","GAT":"D","GGT":"G",
                "GTC":"V","GCC":"A","GAC":"D","GGC":"G",
                "GTA":"V","GCA":"A","GAA":"E","GGA":"G",
                "GTG":"V","GCG":"A","GAG":"E","GGG":"G",
                
                "NNN":"X","---":"-","///":"Z","|||":"Z"}

        return code

    # return a dict with aa as key, nt code as value
    def get_aa_code(self):

        nt_code = self.get_nt_code()

        code = {}
        for i in list(nt_code.keys()):
            if nt_code[i] not in code:
                code[nt_code[i]] = [i]
            else:
                code[nt_code[i]].append(i)

        return code
    
    def get_x4(self,sp=""):
        
        code = {"TCT":"S",
                "TCC":"S",
                "TCA":"S",
                "TCG":"S",
                 
                "CTT":"L","CCT":"P","CGT":"R",
                "CTC":"L","CCC":"P","CGC":"R",
                "CTA":"L","CCA":"P","CGA":"R",
                "CTG":"L","CCG":"P","CGG":"R",

                "ACT":"T",
                "ACC":"T",
                "ACA":"T",
                "ACG":"T",

                "GTT":"V","GCT":"A","GGT":"G",
                "GTC":"V","GCC":"A","GGC":"G",
                "GTA":"V","GCA":"A","GGA":"G",
                "GTG":"V","GCG":"A","GGG":"G"}
                    
        if sp == "fungi":
            code["CTG"] = "S"

        return code
    ##
    # This function is for editing predicted pseudogene cDNA sequences from
    # Torrents et al. Rid of out of frame nts in cDNA. Assuming:
    #
    # 1. The first codon is in frame in the cDNA sequence.
    # 2. cDNA and pep sequence names are the same
    #
    #
    ##
    def pseudo_cds(self,pep,cdna):
        
        print("Make sure the order is correct, or this thing won't work\n")
        
        print("Read sequences...")
        print(" PEP :",pep)
        print(" cDNA:",cdna)
        pdict = manager.fasta_to_dict(pep,0)
        cdict = manager.fasta_to_dict(cdna,0)
        pcode = self.get_aa_code()
        oup   = open("%s.mod" % cdna,"w")
        
        print("Compare pep and cdna sequences...")
        c = 0
        for i in pdict:
            if c % 1000 == 0:
                print(" %i x 1000" % (c/1000))
            c += 1
            
            try:
                pseq = pdict[i]
                cseq = cdict[i]
            except KeyError:
                print("Seq missing:",i) 
            
            cdna_mod = ""
            pos  = 0
            for j in pseq:
                #print j,
                if j in pcode:
                    #print pcode[j],pos
                    code = pcode[j]
                    while pos < len(cseq):
                        # code found
                        if cseq[pos:pos+3] in code:
                            #print "  y:",pos,cseq[pos:pos+3]
                            cdna_mod += cseq[pos:pos+3]
                            pos += 3
                            break
                        # not found, increment by 1
                        else:
                            #print "  n:",pos,cseq[pos:pos+3]
                            pos += 1
                else:
                    #print ":code invalid"
                    pass
            
            oup.write(">%s\n%s\n" % (i,cdna_mod))
        
        print("Done!")
    
    #
    # Assume gene and pep has the same sequence id. Try dynamic programming
    #
    #
    def gene_2_cds(self,gene,pep):
                
        pdict = manager.fasta_to_dict(pep,0)
        gdict = manager.fasta_to_dict(gene,0)
        ntcode= self.get_nt_code()
        for i in gdict:
            print(i)
            if i not in pdict:
                print("ERR - not in pdict:",i)
            else:
                g = gdict[i]
                p = pdict[i]
                a = []
                #print len(p),p
                #print len(g),g
                
                ########
                # INITIALIZATION
                ########
                for i in range(len(g)-2):
                    a.append([0]*len(p))
                
                # Assume the 1st codon is definitely going to be a match.
                # Set score = 1 for all 3 position of the 1st codon
                # gap opening and extension penalty is 0
                a[0][0] = 1
                a[1][0] = 1
                a[2][0] = 1

                ########
                # MATRIX FILL
                ########
                for j in range(1,len(a)):        # iterate row, nt
                    nt = g[j:j+3]
                    aa = ntcode[nt]
                    #print j,nt,aa
                    for k in range(1,len(a[j])): # iterate column, aa
                        S = 0
                        if p[k] == aa:
                            S = 1
                        # get max[M(i-1)(j-1), M(i)(j-1), M(i-1)(j)]
                        # but also M(i-4)(j-1) which is the potential codon
                        # position.
                        
                        m1 = a[j-1][k-1]+S
                        m2 = a[j][k-1]
                        m3 = a[j-1][k]
                        if j-4 < 0:
                            m4 = 0
                        else:
                            m4 = a[j-4][k-1]+S
                        
                        a[j][k] = max(m1,m2,m3,m4)
                        #print "",k,p[k],a[j][k]
                
                c = 0
                print("   0  1  2  3")
                for j in a:
                    print(c,j)
                    c += 1
                
                ########
                # TRACEBACK
                ########
                # Assume last position is definitely a match
                nt = len(a)-1    # starting nt idx for the last codon 
                aa = len(a[0])-1 # last aa idx
                
                path = self.walk(a,nt,aa,[[[nt,aa]]])
                
                matched = []
                align = []
                ########
                # Sequence output
                ########
                aseq = gseq = ""
                x = y = 0
                for j in path:
                    j.reverse()
                    for k in j:
                        aa = ntcode[g[k[0]:k[0]+3]]
                        
                        
                        if aa == p[k[1]] and k[1] not in matched:
                            print("%s:%s" % (g[k[0]],aa))
                            print("  ",aseq,gseq)
                            if aseq != "":
                                align.append([aseq,gseq])
                            matched.append(k[1])
                            aseq = aa
                            gseq = g[k[0]]
                        else:
                            print("%s -"  % g[k[0]])
                            gseq += g[k[0]]
                        
                
                print("%s:%s" % (g[k[0]],aa))
                print("  ",aseq,gseq)
                align.append([aseq,g[-3:]])
                
                # correct codons with less than one nt
                for j in range(len(align)-1):
                    print(j,align[j])
                    clen = len(align[j][1])
                    if clen != 3:
                        print(align[j][1])
                        print(align[j+1][1])
                        align[j][1] = align[j][1] + align[j+1][1][:3-clen]
                        align[j+1][1] = align[j+1][1][3-clen:]
                        aa = ntcode[align[j][1]]
                        if aa == align[j][0]:
                            print("  ok")
                        else:
                            print("  problem")
                for j in align:
                    print(j)
                    
        print("Done!")
    
    #
    # Traceback step of dynamic programming
    #           
    def walk(self,a,i,j,path):
        #print i,j
    
        if i-1 >= 0 or j-1 >= 0:
            if i-1 <0 or j-1 <0: # diagonal score
                m1 = 0
            else:
                m1 = a[i-1][j-1]   
            
            if j-1 < 0:           # row score
                m2 = 0
            else:
                m2 = a[i][j-1]
                
            if i-1 < 0:           # column score
                m3 = 0
            else:
                m3 = a[i-1][j]
            #print "m1,m2,m3:",m1,m2,m3
            
            # diag smaller than row and col scores that are identical
            if m1 < m2 and m1 < m3 and m2 == m3:
                p = []
                # go down two path, duplicate
                for x in range(len(path)):
                    # clone
                    p.append(path[x][:]) 
                    # only if the last pair is the same as i,j
                    if p[-1][-1] == [i,j]:
                        p[-1].append([i,j-1])
                    
                    # clone
                    p.append(path[x][:])    
                    # only if the last pair is the same as i,j
                    if p[-1][-1] == [i,j]:
                        p[-1].append([i-1,j])
                path = p
                #print "1:",path
    
                path = self.walk(a,i,j-1,path)
                path = self.walk(a,i-1,j,path)
            else:
                # diag larger or equal to both row and col
                if m1 >= m2 and m1 >= m3:                           
                    for x in range(len(path)):
                        # check which list to append
                        if path[x][-1] == [i,j]:
                            path[x].append([i-1,j-1])
                    #print "2:",path
                    path = self.walk(a,i-1,j-1,path)
                # row > col
                elif m2 > m3:
                    for x in range(len(path)):
                        # check which list to append
                        if path[x][-1] == [i,j]:
                            path[x].append([i,j-1])
                    #print "3:",path
                    path = self.walk(a,i,j-1,path)
                # col > row
                elif m2 < m3:
                    for x in range(len(path)):
                        # check which list to append
                        if path[x][-1] == [i,j]:
                            path[x].append([i-1,j])
                    #print "4:",path
                    path = self.walk(a,i-1,j,path)
        
        return path
    
    def rmlb(self,astr):
        if astr[-2:] == "\r\n":
            astr = astr[:-2]
        elif astr[-1] == "\n":
            astr = astr[:-1]
        return astr

    def help(self):

        print(" -f   function:")
        print("       bt - based on pep sequence to get the cds")
        print("          REQUIRES: pep,cds,exclstop")
        print("       tl - Translate a nt seq. REQUIRES: -cds")
        print("          OPTIONAL: -id, -frame, -discard, sp")
        print("       tl_mindless - take an alignment pair file from ParseBlast")
        print("          parse_align4 and translate the sequence pairs STRAIGHT")
        print("          NEED:align_seq")
        print("       validate  - compare cds and pep file to see if they are")
        print("          synchronized. Requires: -cds, -pep. Optional:")
        print("          frame")
        print("       bl2_pairs - blast pairs of sequences. REQUIRES: fasta,")
        print("          pairs,blast_dir,prog")
        print("       sixpack - 6 frame tl for ONE sequence, REQUIRES: fasta,")
        print("          OPTIONAL:T,T2,m,c,verbose")
        print("       sixpack_simple - 6 frame translation for multiple seq.")
        print("          do not care about stop, NEED: fasta")
        print("       batch_6pack - pass more than one sequences at a time.")
        print("          Same requirements as sixpack")
        print("       suborf -get alternative orfs within orfs that start with")
        print("          Met. NEED: cds, OPT: T, check")
        print("       exclude - generate a list of ORFs after exclusion. NEED:")
        print("          oc,ec")
        print("       rc - reverse complement of a seq. REQUIRES: fasta")
        print("       rc2 - reverse complement based on names, need: fasta ")
        print("       batch_rc - batch rc, NEED: fasta, OPT: name")
        print("       gene_2_cds - need: gene,pep")
        print(" -align_seq output of ParseBlast.parse_align4")
        print(" -pep peptide seq fasta file, for back_translate, the seq header")
        print("      have to specify the starting pos of the corresponding cds")
        print(" -cds fasta file for cds")
        print(" -gene gene sequence")
        print(" -id  seq id, for translating one sequence in a fasta file")
        print(" -frame 0 [default], 1, or 2")
        print(" -fasta seqeunce file")
        print(" -pairs a list of pairs")
        print(" -discard throw away erroneous translations")
        print(" -T   length cutoff, lower bound, default 25")
        print(" -T2  length threshold, upper bound, default 10000")
        print(" -m   1st met as the beginning of protein seq [1, default] or ")
        print("      get all 'sub-ORFs' starting with met in the full length")
        print("      ORFs [2] (NOT WORKING!!!!) or get the longest [0]")
        print(" -c   increment of window size for sixpack, default 30 (kb)")
        print(" -oc  orf coordinates")
        print(" -ec  exclusion coordinates")
        print(" -exclstop - default no [0], or rid of stop [1]")
        print(" -sp  empty for universal code, f for fungi")
        print(" -check check for length agreement.")
        print(" -verbose give additional info [1] or not [0,default]")
        print(" -n   file with names of sequences to modify")
        print("")
        sys.exit(0)


#-------------------------------------------------------------------------------

trans = translate()

if __name__ == '__main__':

    function = pep = cds = id = fasta = pairs = blast_dir = prog = oc = ec = \
               gene = name = align_seq = ""
    frame   = discard = check = verbose = 0
    m       = 1
    T       = 25
    T2      = 10000
    inc     = 30
    manager = FastaManager.fasta_manager()
    parser  = ParseBlast.parser()
    blastu  = BlastUtility.blast_util()
    exclstop = 0
    
    for i in range(1,len(sys.argv),2):
        if sys.argv[i] == "-f":
            function   = sys.argv[i+1]
        elif sys.argv[i] == "-pep":
            pep = sys.argv[i+1]
        elif sys.argv[i] == "-cds":
            cds = sys.argv[i+1]
        elif sys.argv[i] == "-gene":
            gene         = sys.argv[i+1]
        elif sys.argv[i] == "-id":
            id   = sys.argv[i+1]
        elif sys.argv[i] == "-frame":
            frame      = int(sys.argv[i+1])
        elif sys.argv[i] == "-discard":
            discard    = int(sys.argv[i+1])
        elif sys.argv[i] == "-fasta":
            fasta      = sys.argv[i+1]
        elif sys.argv[i] == "-pairs":
            pairs      = sys.argv[i+1]
        elif sys.argv[i] == "-prog":
            prog       = sys.argv[i+1]
        elif sys.argv[i] == "-blast_dir":
            blast_dir  = sys.argv[i+1]
        elif sys.argv[i] == "-T":
            T          = int(sys.argv[i+1])
        elif sys.argv[i] == "-T2":
            T2         = int(sys.argv[i+1])
        elif sys.argv[i] == "-m":
            m          = int(sys.argv[i+1])
        elif sys.argv[i] == "-c":
            inc        = int(sys.argv[i+1])
        elif sys.argv[i] == "-oc":
            oc         = sys.argv[i+1]
        elif sys.argv[i] == "-ec":
            ec         = sys.argv[i+1]
        elif sys.argv[i] == "-exclstop":
            exclstop = int(sys.argv[i+1])
        elif sys.argv[i] == "-check":
            check = int(sys.argv[i+1])
        elif sys.argv[i] == "-verbose":
            verbose = int(sys.argv[i+1])
        elif sys.argv[i] == "-n":
            name = sys.argv[i+1]
        elif sys.argv[i] == "-align_seq":
            align_seq = sys.argv[i+1]
        else:
            print("\nUnknown flag:",sys.argv[i])
            print(" -h for help")
                    
    if function == "bt":
        if pep == "" or cds == "":
            print("\nRequires fasta files for polypeptide and coding sequences\n")
            trans.help()
        trans.bt(pep,cds,exclstop)
    elif function == "tl":
        if cds == "":
            print("\nRequires fasta file\n")
            trans.help()
        trans.translate(cds,id,frame,discard)
    elif function == "tl_mindless":
        if align_seq == "":
            print("\nRequires alignment seq file\n")
            trans.help()
        trans.tl_mindless(align_seq,T)
    elif function == "validate":
        if cds == "" or pep == "":
            print("\nRequires cds and pep files\n")
            trans.help()
        trans.validate(cds,pep,frame)
    elif function == "bl2_pairs":
        if "" in [fasta, pairs, blast_dir, prog]:
            print("\nRequires fasta, pair list, blast_dir, and prog\n")
            trans.help()
        trans.bl2_pairs(fasta,pairs,blast_dir,prog)
    elif function == "pseudo_cds":
        if "" in [pep, cds]:
            print("\nRequires pep and cds\n")
            trans.help()
        trans.pseudo_cds(pep,cds)   
    elif function == "sixpack":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()
        trans.sixpack(fasta,T,T2,m,inc,verbose)
    elif function == "sixpack_simple":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()
        trans.sixpack_simple(fasta)
    elif function == "suborf":
        if "" in [cds]:
            print("\nRequires cds seq\n")
            trans.help()
        trans.suborf(cds,T,check)
    elif function == "exclude":
        if "" in [ec,oc]:
            print("\nRequires excluded and orf coord files\n")
            trans.help()
        trans.exclude(oc,ec)
    elif function == "rc":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()
        trans.rc(fasta,call=1)
    elif function == "rc2":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()
        trans.rc2(fasta)
    elif function == "batch_rc":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()
        trans.batch_rc(fasta,name)
    elif function == "batch_6pack":
        if "" in [fasta]:
            print("\nRequires fasta\n")
            trans.help()        
        trans.batch_6pack(fasta,T,T2,m,inc,verbose)
    elif function == "help":
        trans.help()
    elif function == "gene_2_cds":
        if "" in [gene,pep]:
            print("\nRequires gene and pep seq files\n")
            trans.help()        
        trans.gene_2_cds(gene,pep)
    elif function == "test":
        trans.test()
    else:
        if function == "":
            print("\nNeed to specify function\n")
        else:
            print("\nUnknown function:",function)
        trans.help()

"""
Deprecated 01/04,2005

    def complement2(self,nt):

        c = {"A":"T","T":"A",
              "C":"G","G":"C",
              "N":"N"}
        c_nt = ""
        for i in nt:
            c_nt = c_nt + c[i]

        return c_nt
    # throw an error in mer:
    # TypeError: sequence index must be integer
    def reverse4(self,a):
        t1 = time.time()
        a = a[::-1] 
        t2 = time.time()
        print t2-t1
        return a
    
    # fast but very memory hungry
    def reverse3(self,nt):
        f = lambda s:s and f(s[1:])+s[0]
        return f(nt)

    def reverse(self,nt):

        r_nt = ""
        for i in nt:
            r_nt = i + r_nt

        return r_nt
"""
