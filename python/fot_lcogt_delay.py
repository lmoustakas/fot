#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
Full Of Time
'''

import matplotlib
matplotlib.use('Agg') # Note, this MUST be before importing pylab or matplotlib.pyplot


if __name__ == "__main__":
    
	from fot_library import * 
	import argparse
	parser=argparse.ArgumentParser(description='fot_delay routine to calculate delay inference')
	parser.add_argument("-i","--datafile",	help="time delay challenge data file",type=str)
	parser.add_argument("-im1","--image1",	default='A', help="Image 1 name (e.g. 'A')",type=str)
	parser.add_argument("-im2","--image2",	default='B', help="Image 2 name (e.g. 'B')",type=str)
	parser.add_argument("-su","--systematic",	default=0., help="Additional systematic uncertainty (e.g. 0.1) [mag]",type=float)


	parser.add_argument("-dtp",   "--delay_prior",          default = 0., help="Delay prior for emcee",type=float)
	parser.add_argument("-dtpmin","--delay_prior_min",      default = -1.5, help="Delay prior minimum acceptable value for emcee",type=float)
	parser.add_argument("-dtpmax","--delay_prior_max",      default = +1.5, help="Delay prior maximum acceptable value for emcee",type=float)
	parser.add_argument("-dmp",   "--delta_mag_prior",      default = 0., help="Difference in magniitude between images prior for emcee",type=float)
	parser.add_argument("-dmpmin","--delta_mag_prior_min",  default = -10., help="Difference in magniitude between images prior minimum acceptable value for emcee",type=float)
	parser.add_argument("-dmpmax","--delta_mag_prior_max",  default = +10., help="Difference in magniitude between images prior maximum acceptable value for emcee",type=float)
	parser.add_argument("-sp",     "--sigma_prior",         default = 0.1, help="Quasar light curve sigma prior for emcee",type=float)
	parser.add_argument("-spmin",  "--sigma_prior_min",     default = 7.e-5, help="Quasar light curve sigma prior minimum acceptable value for emcee",type=float)
	parser.add_argument("-spmax",  "--sigma_prior_max",     default = 7.e1, help="Quasar light curve sigma prior maximum acceptable value for emcee",type=float)
	parser.add_argument("-tp",     "--tau_prior",           default = 100., help="Quasar light curve prior for emcee",type=float)
	parser.add_argument("-tpmin",  "--tau_prior_min",       default = 1., help="Quasar light curve prior minimum acceptable value for emcee",type=float)
	parser.add_argument("-tpmax",  "--tau_prior_max",       default = 1.e5, help="Quasar light curve prior maximum acceptable value for emcee",type=float)
	parser.add_argument("-mp",    "--avg_mag_prior",        default = 18., help="Quasar light curve magnitude prior for emcee",type=float)
	parser.add_argument("-mpmin", "--avg_mag_prior_min",    default = 0., help="Quasar light curve magnitude prior minimum acceptable value for emcee",type=float)
	parser.add_argument("-mpmax", "--avg_mag_prior_max",    default = 100., help="Quasar light curve magnitude prior maximum acceptable value for emcee",type=float)
	parser.add_argument("-z",      "--redshift",            default = 0., help="Quasar redshift",type=float)
	parser.add_argument("-o",      "--output_tag",          default = 'tmp', help="Output tag, in quotes",type=str)
	parser.add_argument("-od",     "--outputdir",           default = './', help="Output directory, in quotes",type=str)

	#KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
	for i, arg in enumerate(sys.argv):
 		if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

	args=parser.parse_args()
	# define the right column header names
	print args
	# read in the data
	#print rung, args.datafile
	import os
	fname = args.datafile
	print 'checking path', fname
	if(os.path.isfile(fname)==False):
		print fname
		print 'File not found.'
		exit()

	data=ascii.read(fname)
	#print data
	print len(data), 'rows of data'
	print ''
	#exit()
	t=data['mjd']
	m1=data['mag_%s'%args.image1]
	em1=data['magerr_%s'%args.image1]
	m2=data['mag_%s'%args.image2]
	em2=data['magerr_%s'%args.image2]
    
	emcee_delay_estimator(t, m1, em1, m2, em2 ,args.output_tag+'_dt%s%s'%(args.image1,args.image2), 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  args.sigma_prior,	args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,	args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max,
			  outputdir = args.outputdir)

