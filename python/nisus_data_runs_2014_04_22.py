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


data_file_list = ['HE0435_Courbin2011', 'HS2209_Eulaers2013', 'J1206_Eulaers2013', 'J1001_Rathnakumar2013', 'J1206_Eulaers2013', 'RXJ1131_Tewes2013']
ph_un_add_list = ['0p10','0p05','0p02', '0p01', '0p0']

for fnm in data_file_list:
  image_list=[]
  if 'mag_A' in open('../data/cosmograil/'+fnm+'.rdb').read(): image_list.append('A')
  if 'mag_B' in open('../data/cosmograil/'+fnm+'.rdb').read(): image_list.append('B')
  if 'mag_C' in open('../data/cosmograil/'+fnm+'.rdb').read(): image_list.append('C')
  if 'mag_D' in open('../data/cosmograil/'+fnm+'.rdb').read(): image_list.append('D')
  print fnm, image_list
  ls = os.listdir("../data/cosmograil/")

  for i in range(0,len(image_list)):
    for j in range(i+1,len(image_list)):
      print '\t',i,j,image_list[i],image_list[j]
      for ph_un_add in ph_un_add_list:
          com = './fot_delay.py -i %s.rdb -l %s -m %s -su %1.3f -dtp 0. -dtpmin -200. -dtpmax 200. -dmp 0.5 -dmpmin -5. -dmpmax 5. -sp 0.07 -spmin 0.00007 -spmax 7. -tp 121. -tpmin 10. -tpmax 1000. -mp -13. -mpmin -30 -mpmax 100. -z 0. -o %s_%s_%s_%s > %s_%s_%s_%s.log &'%(fnm, image_list[i], image_list[j], float(ph_un_add.replace('p','.')),fnm,image_list[i], image_list[j], ph_un_add, fnm, image_list[i], image_list[j], ph_un_add)
          print ''
          print datetime.datetime.now()
          print com
          run_queue(njobs,com)	



