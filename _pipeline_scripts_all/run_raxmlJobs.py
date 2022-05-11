import sys,os
ccount=int(sys.argv[1])
os.system('/home/moghegau/packages/RAxML-7.2.6/raxmlHPC-PTHREADS '\
                          '-f d -m PROTGAMMAJTT -T 6 -s tmpRAX%i.phy '\
                          '-n tmpRAX%i.phy.raxout' %(ccount,ccount))
os.system('mv tmpRAX%i.phy _TMP_PHYFILES\n'%ccount)
os.system('mv *_result.tmpRAX%i.phy.raxout _FINAL_RESULTS\n'%ccount)
os.system('mv *_info.tmpRAX%i.phy.raxout _TMP_INFOFILES\n'%ccount)
os.system('mv *_log.tmpRAX%i.phy.raxout _TMP_LOGFILES\n'%ccount)
os.system('mv *_parsimonyTree.tmpRAX%i.phy.raxout _TMP_PARSIMONY\n'%ccount)
os.system('mv *_bestTree.tmpRAX%i.phy.raxout _TMP_BESTTREE\n'%ccount)



