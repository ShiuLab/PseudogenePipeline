#This script is designed to compare two files to see how they differ.  There are
#currently two options -o (y or n) [ordered] and -p (y or n)[print, prints the
#different lines]

import sys

class Compare2files:
    def CompareInternal(self, big, small, num, cnt):
        for line in small:
            if line != big[cnt]:
                num += 1
                if p == "y":
                    print "S: ", line[:-1]
                    print "B: ", big[cnt]
                elif p == "n":
                    pass
                else:
                    print "improper -p.  Should be 'y' or 'n'"
            cnt += 1
        for c in range(cnt, len(big)):
            num += 1
            if p == "y":
                print "B: ", big[c]
            elif p == "n":
                pass
            
        return num

    def Ordered(self, num):
        file1Lst = file1.readlines()
        file2Lst = file2.readlines()

        #compare files line-by-line
        cnt = 0
        print "note: sometimes the script screwes up with the first line"
        print "\nnote: S and B represent the smaller and bigger files\n"
        if p == "y":
            print "differing lines: \n"
        if len(file1Lst) < len(file2Lst):
            print "zero"
            small = file1Lst
            big = file2Lst
            print "small: file2"
            print "big:   file1\n"
            num = Comp.CompareInternal(big, small, num, cnt)
            
        elif len(file1Lst) > len(file2Lst):
            small = file2Lst
            big = file1Lst
            print "small: file1"
            print "big:   file2\n"
            num = Comp.CompareInternal(big, small, num, cnt)
                
        else:
            for cnt in range(len(file1Lst)):
                if file1Lst[cnt] != file2Lst[cnt]:
                    num += 1
                    if p == "y":
                        print file1Lst[cnt][:-1]
                    elif p == "n":
                        pass
                    else:
                        print "improper -p.  Should be 'y' or 'n'"
                cnt += 1
        print "%d differing lines." % (num)
            
        
    def Unordered(self, num):
        file1Lst = file1.readlines()
        file2Lst = file2.readlines()
        print len(file1Lst)
        print len(file2Lst)

        print "note: sometimes the script screwes up with the first line"
        print "note: all printed lines are from file1"

        for line in file1Lst:
            print line
            if line not in file2Lst:
                num += 1
                if p == "y":
                    print line
                elif p == "n":
                    pass
                else:
                    print "improper -p.  Should be 'y' or 'n'"
                
        print "%d differing lines." % (num)
        
    
    def Help(self):
        print "\nParameters:"
        print "    file1 - the first file for comparison."
        print "    file2 - the second file for comparison."
        print "    o - ordered.  Whether the two files should be compared"
        print "        line-by-line in order, y or n"
        print "    p - print.  Whether the differing lines should be printed, y or n"




if __name__ == '__main__':
    Comp = Compare2files()
    num = 0

    for i in range(1,len(sys.argv),2):
        if sys.argv[i] == "-file1":
            file1 = open(sys.argv[i+1])
        elif sys.argv[i] == "-file2":
            file2 = open(sys.argv[i+1])
        elif sys.argv[i] == "-o":
            o = sys.argv[i+1]
        elif sys.argv[i] == "-p":
            p = sys.argv[i+1]
        elif sys.argv[i] == "-h":
            Comp.Help()
        else:
            print "UNKNOWN FLAG:",sys.argv[i]
            print "add -h to get help."
            sys.exit(0)

    if o == "y":
        Comp.Ordered(num)
    elif o == "n":
        Comp.Unordered(num)
    else:
        print "improper -o.  Should be 'y' or 'n'"
        print "add -h to get help."
        
