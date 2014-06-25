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
import glob

#SCRIPT RUNS 25 JOBS AT A TIME.
#IT CHECKS IF JOBS ARE DONE EVERY MINUTE

sys.stdout.flush()
prev_t = time.clock()
#njobs=30
njobs=50


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


rung_list = [0,1,2,3,4]
tdc1_data_dir = os.environ['FOTDIR']+'/data/tdc1/'


#rung_list = rung_list[::-1] # run backwards
for rung in rung_list:
  print tdc1_data_dir + 'rung%d/'%rung + '*.txt'
  fnames = glob.glob(tdc1_data_dir + 'rung%d/'%rung + '*.txt')
  fnames = sorted(fnames)
  counter=0
  for f in fnames:
    counter+=1
    fnm = f.split('/')[len(f.split('/'))-1]
    fnm = fnm.split('.')[0]

    t0 = datetime.datetime.now()
    date_string = t0.strftime("_%Y_%m_%d_%H:%M")
    log_fnm = ''.join([fnm,date_string ])

    print fnm
    #for ph_un_add in ph_un_add_list:
    #com = './fot_delay_tdc1.py -i %s.txt -l \'A\' -m \'B\' -su 0. -dtp 0. -dtpmin -1600. -dtpmax 1600. -dmp 0. -dmpmin -10. -dmpmax 10. -sp 1. -spmin 0.00007 -spmax 70. -tp 3500. -tpmin 10. -tpmax 100000. -mp -13. -mpmin -30 -mpmax 100. -z 0. -o %s > %s.log'%(fnm ,fnm,  log_fnm)
    com = '%s/fot_delay_tdc1.py -i %s.txt -l \'A\' -m \'B\' -su 0. -dtp 0. -dtpmin -1600. -dtpmax 1600. -dmp 0. -dmpmin -10. -dmpmax 10. -sp 1. -spmin 0.00007 -spmax 70. -tp 3500. -tpmin 10. -tpmax 100000. -mp -13. -mpmin -30 -mpmax 100. -z 0. -o %s'%(os.environ['FOTDIR']+'/python', fnm ,fnm)
    print ''
    print datetime.datetime.now()
    f=open('./inputs/rung%d/input.%d'%(rung,counter),'w')
    f.write(com)
    f.close()
    #run_queue(njobs,com)	



