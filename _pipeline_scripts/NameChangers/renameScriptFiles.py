#This script is designed to rename all my scripts starting with R_ to 1_ to
#reflect the name change in my organizational system.

import sys
import os

folder = sys.argv[1] #the name of the folder to work in

for fileName in os.listdir(folder):
    if fileName.startswith("R_") and fileName.endswith(".py"):
        newName = "1_" + fileName[2:]
        os.system("mv %s %s" % (fileName, newName))
