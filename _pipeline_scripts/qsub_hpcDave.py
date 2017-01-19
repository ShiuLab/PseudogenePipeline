
import sys, os, time

class qsub_hpc:
	
	## private method
	def submit(self,jobs,sidx,walltime,memory,jobname,ncpu,nnodes,cusmod,email,logdir,wdir=""):
		
		runtype = 0
		if type(jobs) != list:
			runtype = 1
			inp = open(jobs)
			inl = inp.readlines()
			tmp = []
			for i in inl:
				i = i.strip()
				tmp.append(i)
			jobs = tmp

		if logdir == 0:
			logdir = "-o tmp.o -e tmp.e "
		else:
			logdir = ""

		h = walltime/60
		m = walltime%60
	
		for i in jobs:
			print "  job %i" % sidx
			
			oup = open("%s%i.sh" % (jobname,sidx),"w")
			oup.write("#!/bin/sh -login\n\n#PBS -q main\n")
			oup.write("#PBS -l nodes=%i:ppn=%i,walltime=%i:%i:00,mem=%igb\n" % \
                                               (nnodes,ncpu,h,m,memory))
			if email != "":
				oup.write("#PBS -M %s\n" % email)
				oup.write("#PBS -m a\n")
			if wdir != "":
				oup.write("#PBS -d %s\n" % wd)
			if cusmod != "":
                                csp=cusmod.split(',')
                                for cus1 in csp:
                                        oup.write('module load %s\n'%(cus1))


			oup.write("%s\n" % self.rmlb(i))
			oup.close()
			os.system("chmod 755 %s%i.sh" % (jobname,sidx))			
			os.system("qsub %s%s%i.sh" % (logdir,jobname,sidx))	
			sidx += 1
	
	def quene(self,jcommand,stime,nsub,juser,jrange,walltime,memory,jobname,ncpu,nnodes,cusmod,email,logdir,wdir=""):
		
		jdict = {}
		inp = open(jcommand)
		jobs = inp.readlines()
		if jrange != "*":
			jrange = jrange.split("-")
			jobs  = jobs[int(jrange[0])-1:int(jrange[1])]
		
		oup = open("%s.log" % jcommand,"w")					
		j = 0
		while j < len(jobs): 
			oup.write("%s\t" % time.ctime())
			print "%s\t" % time.ctime()
			
			#####
			# Check number of jobs in quene so far
			#####
			os.system("qstat -u %s > TMPjoblist" % juser)
			inp = open("TMPjoblist")
			inl = inp.readlines()
			if len(inl) != 0:
				currj = len(inl)-5
			else:
				currj = 0
			
			#####
			# Submit more jobs so there is always nsub number of jobs in quene
			#####
			if currj < nsub:
				jseg = jobs[j:j+(nsub-currj)]
				oup.write("submit %i\t" % len(jseg))
				print "current %i, submit %i\n" % (currj,len(jseg))
				self.submit(jseg,j+1,walltime,memory,jobname,ncpu,nnodes,cusmod,email,logdir,wdir)
				j += len(jseg)
			else:
				oup.write("waiting\t")
				print "waiting\t"
			oup.write("submited so far: %i\n" % j)
			print "submited so far: %i\n" % j
			time.sleep(stime)
		
		os.system("rm TMPjoblist")
		print "Done!"

	def qdel(self,duser,drange):
		print [duser,drange]
		#if "" not in [duser,drange]:
		#	print "Specify either user OR range! Quit!...\n"
		#	sys.exit(0)
		
		print "User :",duser
		print "Range:",drange
		
		if duser != "":
			os.system("qstat -u %s > TMP_qdel" % duser)
			inp = open("TMP_qdel")
			inl = inp.readlines()
			dlist = []
			for i in inl:
				dlist.append(i.split(".")[0])
			for i in dlist:
				try:
					print "qdel %s" % i
					os.system("qdel %s" % i)
				except ValueError:
					pass
			os.system("rm TMP_qdel")
		else:
			drange = drange.split("-")
			b = int(drange[0])
			e = int(drange[1])+1 # inclusive
			for i in range(b,e):
				print "qdel %s" % i
				os.system("qdel %i" % i)
			
		print "Done!"

	def check_running(self,user):
		
		os.system("qstat > TMP_check_running_log")
		inp = open("TMP_check_running_log")
		inl = inp.readlines()
		countQ = 0
		countR = 0
		countE = 0
		countH = 0
		countUserH = 0
		countUserQ = 0
		countUserR = 0
		countUserE = 0
		flagStart = 0
		flagUser = 0
		if user != "":
			flagUser = 1
		for i in inl:
			if i[:3] == "---":
				flagStart = 1
			elif flagStart:
				i = i.split(" ")
				tmp = []
				for j in i:
					if j != "":
						tmp.append(j)
				if tmp[4] == "R":
					countR += 1
					if flagUser and tmp[2] == user:
						countUserR += 1
				elif tmp[4] == "Q":
					countQ += 1
					if flagUser and tmp[2] == user:
						countUserQ += 1
				elif tmp[4] == "E":
					countE += 1
					if flagUser and tmp[2] == user:
						countUserE += 1
				elif tmp[4] == "H":
					countH += 1
					if flagUser and tmp[2] == user:
						countUserH += 1
				else:
					print "UNK:",tmp[4]

		print "Running :",countR
		print "Queueing:",countQ
		print "Held    :",countH
		print "Err     :",countE
		if flagUser:
			print "User: %s,  %i running, %i in queue, %i held, %i with error" % \
									(user,countUserR,countUserQ,countUserH,countUserE)

		os.system("rm TMP_check_running_log")
	
	def check_err(self,jobname):
		files = os.listdir("./")
		
		# err log files
		efiles = []
		for i in files:
			if i.find(jobname) != -1 and i.find(".sh.e") != -1:
				efiles.append(i)
		efiles.sort()

		# go through err log files
		print "With error:"
		eidx = []
		for i in efiles:
			s = os.path.getsize("./%s" % i)
			try:
				idx = int(i[i.find(jobname)+len(jobname):i.find(".sh.e")])
				if s != 0:
					inp = open(i)
					inl = inp.readlines()
					print " %s:" % i,inl
					eidx.append(idx)
			except ValueError:
				print " %s:" % i,"job index err"
				continue
		
		print "###############"
		print "Compile new command line. Note that this rely on the presence"
		print "of the shell script file in the same folder as the err log."
		print "In addition, it just look for the last non-empty lines in .sh."
		print "###############"
		print " %i jobs failed" % len(eidx)
		oup = open("cmd_%s_witherr" % (jobname),"w")
		for i in eidx:
			# read the shell script
			inp = open("%s%i.sh" % (jobname,i))
			inl = inp.readlines()
			
			# Look for the last non-empty line
			inl.reverse()
			for j in inl:
				j = j.strip()
				if j != "":
					oup.write("%s\n" % j)
					break
			
		print "Command lines of jobs with err: cmd_%s_witherr" % jobname

		print "Done!"	

	# private method
	def rmlb(self,astr):
		astr = astr.strip()
		return astr

	def help(self):
		print "\nFunctions (-f):"
		print "    submit - create shell script and submit job bsaed on a file"
		print "       where each line is one job. NEED: c, OPT: w,m,J,p,ei,wd"
		print "    quene - submit jobs sequentially. NEED:c,u, OPT:s,n,r,w,m,"
		print "       J,p,e,k"
		print "    qdel  - delete jobs, NEED: r or u (if want to kill all)"
		print "    check_running - check how many job have R stat, OPT: u"
		print "    check_err - check any job error and compile a new cmd line"
		print "       file,  NEED: J"
		print ""
		print "Parameters:"
		print "    c - the command line file, one line per job"
		print "    e - email, default '', if an email is passed, PBS -m will be"
		print "        set to a."
		print "    s - sec between qstat check, default 10"
		print "    n - number of jobs to submit at a time , default 10"
		print "    u - which user to monitor"
		print "    cusmod - custom module name (module load cusmod)"
		print "             Separate multiple module names by comma (mod1,mod2,mod3)"
		print "    r - jobnum1-jobnum2, jobnum-, or all (default)"
		print "    w - walltime, in minutes. Defaulit 10 min."
		print "    wd- working dir"
		print "    m - memory in GB, default 1"
		print "    J - job name"
		print "    nnodes - number of nodes, default 1"
		print "    p - number of cpu (ppn), default 1"		
		print "    o - keep all log files (1) or rename to tmp.o/tmp.e (0)"
		print ""
		sys.exit(0)

if __name__ == '__main__':

	qsub = qsub_hpc()
	f = c = u = e = wd = cus1 = ""
	s = n = w = 10
	m = p = o = nn1 = 1
	r = "*"
	J = "job"
		
	for i in range(1,len(sys.argv),2):
		if sys.argv[i] == "-f":
			f  = sys.argv[i+1]
		elif sys.argv[i] == "-c":
			c  = sys.argv[i+1]
		elif sys.argv[i] == "-e":
			e  = sys.argv[i+1]
		elif sys.argv[i] == "-cusmod":
			cus1  = sys.argv[i+1]
		elif sys.argv[i] == "-nnodes":
			nn1  = int(sys.argv[i+1])
		elif sys.argv[i] == "-s":
			s  = int(sys.argv[i+1])
		elif sys.argv[i] == "-n":
			n  = int(sys.argv[i+1])
		elif sys.argv[i] == "-u":
			u  = sys.argv[i+1]
		elif sys.argv[i] == "-r":
			r  = sys.argv[i+1]
		elif sys.argv[i] == "-w":
			w  = int(sys.argv[i+1])
		elif sys.argv[i] == "-m":
			m  = int(sys.argv[i+1])
		elif sys.argv[i] == "-J":
			J  = sys.argv[i+1]
		elif sys.argv[i] == "-p":
			p  = int(sys.argv[i+1])
		elif sys.argv[i] == "-o":
			o  = int(sys.argv[i+1])
		elif sys.argv[i] == "-wd":
			wd = sys.argv[i+1]
		else:
			print "UNKNOWN FLAG:",sys.argv[i]
			print "Do -h to get help."
			sys.exit(0)

	if f == "submit":
		if "" in [c]:
			print "\nNeed cmd line file, user name\n"
			qsub.help()
		qsub.submit(c,0,w,m,J,p,nn1,cus1,e,o,wd)
	elif f == "quene":
		if "" in [c,u]:
			print "\nNeed cmd line file, user name\n"
			qsub.help()
		qsub.quene(c,s,n,u,r,w,m,J,p,nn1,cus1,e,o,wd)
	elif f == "qdel":
		if r == "" and u == "":
			print "\nNeed range or user name\n"
			qsub.help()
		qsub.qdel(u,r)
	elif f == "check_running":
		if u == "":
			print "\nNeed (u)ser id!\n"
			qsub.help()
		qsub.check_running(u)
	elif f == "check_err":
		if J == "job":
			print "\nYou are using the default job name, make sure this is really what"
			print "you want to check... Go ahead anyway"
		qsub.check_err(J)
	else:
		print "\nUnknown function...\n"
		qsub.help()
