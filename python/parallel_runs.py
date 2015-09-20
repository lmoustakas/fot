#!/usr/bin/env python
if __name__ == "__main__":
    
	import multiprocessing
	import os
	import argparse

	parser=argparse.ArgumentParser(description="Run fot_hubble simulations and analysis in parallel")
	parser.add_argument("-rl", "--run_list", help="file with the commands to be executed", type=str, default = '')
	parser.add_argument("-np", "--numprocesses", help="number of processors active", type=int, default = 1)
	parser.add_argument("-ni", "--niceness", help="set the niceness for these runs", type=int, default = 20)

	args=parser.parse_args()

	os.nice(args.niceness) # don't let these runs take over your cpus.

	runs = [] # initialize array
	for line in file(args.run_list):
		if('&' in line): # ensure processes are not set to run in background. This will make all runs activate simultaneously!
			runs.append(line.split('&')[0])
		else:  # remove end of line character from text file list entries
			runs.append(line.split('\n')[0])

	print 'Your CPU has %d processors'%(multiprocessing.cpu_count())
	print 'You are about to run %d runs over %d active processors at a time'%(len(runs), args.numprocesses)
	variable = raw_input('Are you sure you want to continue? [y/n]: ')
	if(variable!='y'):
		print 'Aborting run'
		exit()

	p = multiprocessing.Pool(processes=args.numprocesses)
	print 'Running processes'
	p.map(os.system, runs)



