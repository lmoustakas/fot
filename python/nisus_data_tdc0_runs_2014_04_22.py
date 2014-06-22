#!/usr/bin/env python
'''
2014 April 22
ARW & LAM, JPL/Caltech
Run the analysis on all the cosmograil data (not by season).
Each data set has systematics of 0.10, 0.05, 0.02, 0.01 and 0.00 added.
'''
from pylab import *
import os
import sys
import time
#SCRIPT RUNS 25 JOBS AT A TIME.
#IT CHECKS IF JOBS ARE DONE EVERY MINUTE

sys.stdout.flush()
prev_t = time.clock()
njobs=25


def run_queue(numjobs, com):
	os.system(com)
	print com
	bail=0
	while(bail==0):
		time.sleep(30)
		os.system('ps -eLf | grep fot_delay | grep afromero > check')
		lc=0
		for line in file("check"):
		        lc+=1
		print "running:", lc-2, "jobs"
		if(lc<numjobs-2): bail=1
		print '\n'
	os.system('ps -eLf | grep fot_hubble_sim | grep afromero')


rung_list = [0,1,2,3,4,5,6]
pair_list = [1,2,3,4,5,6,7,8]
#ph_un_add_list = ['0p10','0p05','0p02', '0p01', '0p0']

for rung in rung_list:
  for pair in pair_list:
    fnm = 'tdc0_rung%d_pair%d'%(rung,pair)
    t0 = datetime.datetime.now()
    date_string = t0.strftime("_%Y_%m_%d_%H:%M")
    log_fnm = ''.join([fnm,date_string ])

    print '\t', rung, pair, fnm
    #for ph_un_add in ph_un_add_list:
    com = './fot_delay_tdc0.py -i %s.txt -l \'A\' -m \'B\' -su 0. -dtp 0. -dtpmin -1600. -dtpmax 1600. -dmp 0. -dmpmin -10. -dmpmax 10. -sp 1. -spmin 0.00007 -spmax 70. -tp 100. -tpmin 10. -tpmax 10000. -mp -13. -mpmin -30 -mpmax 100. -z 0. -o %s > %s.log &'%(fnm ,fnm,  log_fnm)
    print ''
    print datetime.datetime.now()
    print com
    #run_queue(njobs,com)	



