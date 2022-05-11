#This script is a wrapper for taking a Newick format tree file, extracting node
#information and determining relationships from them.
#Created by David E. Hufnagel on Oct 11, 2012
#modules to load: ete2

import sys, os

base = sys.argv[1]
scriptDir = sys.argv[2]
if "/" in scriptDir:
    scriptDir.strip("/")

os.system("python %s/4_1_makeNodesTabFromNewick.py %s %s.nodes" % (scriptDir, base, base))
#os.system("python %s/4_1_parseNodes.py %s.nodes %s.nodes.parsed" % (scriptDir, base, base))
#os.system("python %s/4_1_determineRelatFromNodes.py %s.nodes.parsed %s.nodes.parsed.info" % (scriptDir, base, base))
