#!/usr/bin/env python
#from analysis_library import *
from pylab import *
import numpy
import os
import glob

rcParams['font.size'] = 10
rcParams['figure.facecolor'] = 'white'
rcParams['axes.formatter.limits'] = -3, 4

def get_chain(filename):
	npzfile = numpy.load(filename)
	#print 'npzfile.files', npzfile.files
	samples = npzfile['arr_0']
	return samples

def plot_chain_iterations(chain, ndim, nwalkers, niterations, niteration_cut, labels=''):
	figure(figsize=(15,12))
	for n in range(0,ndim):
	  subplot(ndim,3,3*n+1)
	  for k in range(0,nwalkers):
	    plot((chain[k*niterations:(k+1)*niterations,n]))
	    xticks(rotation=30)
	    if(len(labels)==ndim):
		ylabel(labels[n])
	  xlabel('iteration')
	  subplot(ndim,3,3*n+2)
	  for k in range(0,nwalkers):
	    plot((chain[k*niterations+niteration_cut:(k+1)*niterations,n]))
	    xticks(rotation=30)
	    if(len(labels)==ndim):
		ylabel(labels[n])
	  xlabel('iteration')
	  subplot(ndim,3,3*n+3)
	  hist(ravel((chain[:,n].reshape((nwalkers, niterations)))[:,niteration_cut:niterations]), bins=100)
	  xticks(rotation=30)
	  if(len(labels)==ndim):
	     xlabel(labels[n])
	subplots_adjust(bottom=0.1, top=0.97, wspace=0.5, hspace=0.5)
	rcParams['font.size'] = 12
	rcParams['figure.facecolor'] = 'white'

def get_conf_interval(chain, dimension, nwalkers, niterations, niteration_cut, conf_interval=0.68, display=False, label=''):
	if (conf_interval < 0. or conf_interval>1.):
		print 'ERROR, conf_interval is a value between 0. and 1. Entry is', conf_interval 
	x = sort(ravel((chain[:,dimension].reshape((nwalkers, niterations)))[:,niteration_cut:niterations]))
	c = arange(len(x), dtype=float64)/len(x)

	newlist = array([ max((x[n:]-x[n])[c[n:]-c[n]<conf_interval]) for n in range(0,int((1.-conf_interval)*nwalkers*niteration_cut))[::200]])
	i1 = x[numpy.where((newlist == min(newlist)))[0]*200][0]
	i2=i1+min(newlist)
	if(display==True):
		figure()
		subplot(311)
		plot(x,c,',')
		plot([i1,i1],[0.,1.],'r--', lw=2)
		plot([(i1+i2)/2.,(i1+i2)/2.],[0.,1.],'r-', lw=2)
		plot([i2,i2],[0.,1.],'r--', lw=2)
		subplot(312)
		a,b,c = hist(x, bins=50)
		plot([i1,i1],[0.,max(a)],'r--', lw=2)
		plot([(i1+i2)/2.,(i1+i2)/2.],[0.,max(a)],'r-', lw=2)
		plot([i2,i2],[0.,max(a)],'r--', lw=2)
		ylim(0.,max(a))
		show()
	#print 'done calculating minimum'
	#print len(newlist)
	#print numpy.where((newlist == min(newlist)))[0]
	return i1,i2,x


if __name__ == "__main__":
    
    import argparse
    parser=argparse.ArgumentParser(description='fot_delay routine to calculate delay inference')
    parser.add_argument("-fbn","--file_base_name",	help="for example /aurora_nobackup/lenssim/romerowo/fot/outputs/tdc1_2014_08_08/chains/rung0/chain_samples_tdc1_rung0_double_pair315",type=str)
    parser.add_argument("-fout","--output_filename",	default='ana_tdc1.out', help="output file name. outputs is appended to this file",type=str)
    parser.add_argument("-display","--display",	default=0, help="display (0 = no display, 1 = display delay and sigma distributions, 2 = display chain iterations",type=int)

    
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    # define the right column header names
    #print args
    

#fnm_lc1 = max(glob.iglob('%s_lc1*.npz'%(args.file_base_name)), key=os.path.getctime)
#fnm_lc2 = max(glob.iglob('%s_lc2*.npz'%(args.file_base_name)), key=os.path.getctime)
print args.file_base_name
#print glob.iglob('%s_dt*.npz'%(args.file_base_name))
#fnm_dt = max(glob.iglob('%s_dt*.npz'%(args.file_base_name)), key=os.path.getctime)
fnm_dt = '../outputs/chain_samples_tdc1_rung0_double_pair182_dt_2014_08_17_21:58.npz'
#print fnm_lc1
#print fnm_lc2
#print fnm_dt

#chain_lc1 = get_chain(fnm_lc1)
#chain_lc2 = get_chain(fnm_lc2)
chain_dt = get_chain(fnm_dt)


ndim, nwalkers, niterations, niteration_cut = [8,100,4000,2000] # these should be set in the command line
if(args.display>1): plot_chain_iterations(chain_dt, ndim, nwalkers, niterations, niteration_cut)
show()
delay_lower, delay_upper, delay = get_conf_interval(chain_dt, 0, nwalkers, niterations, niteration_cut)
delay_lower_95, delay_upper_95, delay_95 = get_conf_interval(chain_dt, 0, nwalkers, niterations, niteration_cut, conf_interval=0.95)
delay_lower_98, delay_upper_98, delay_98 = get_conf_interval(chain_dt, 0, nwalkers, niterations, niteration_cut, conf_interval=0.98)

sig_dt_lower, sig_dt_upper, sig_dt = get_conf_interval(chain_dt, 5, nwalkers, niterations, niteration_cut)

ndim, nwalkers, niterations, niteration_cut = [3,100,2000,1000] # these should be set in the command line
if(args.display>1): plot_chain_iterations(chain_lc1, ndim, nwalkers, niterations, niteration_cut)
if(args.display>1): plot_chain_iterations(chain_lc2, ndim, nwalkers, niterations, niteration_cut)

sig_lc1_lower, sig_lc1_upper, sig_lc1 = get_conf_interval(chain_lc1, 0, nwalkers, niterations, niteration_cut)
sig_lc2_lower, sig_lc2_upper, sig_lc2 = get_conf_interval(chain_lc2, 0, nwalkers, niterations, niteration_cut)

s=((args.file_base_name).split('/')[-1]).split('chain_samples_')[-1]+'\t'
#s+='%1.3f\t'%(0.5*(delay_lower + delay_upper))
s+='%1.3f\t'%(delay_lower)
s+='%1.3f\t'%(delay_upper)
s+='%1.3f\t'%(delay_lower_95)
s+='%1.3f\t'%(delay_upper_95)
s+='%1.3f\t'%(delay_lower_98)
s+='%1.3f\t'%(delay_upper_98)
#s+='%1.3e\t'%(0.5*(sig_dt_lower + sig_dt_upper))
s+='%1.3e\t'%(sig_dt_lower)
s+='%1.3e\t'%(sig_dt_upper)
#s+='%1.3e\t'%(0.5*(sig_lc1_lower + sig_lc1_upper))
s+='%1.3e\t'%(sig_lc1_lower)
s+='%1.3e\t'%(sig_lc1_upper)
#s+='%1.3e\t'%(0.5*(sig_lc2_lower + sig_lc2_upper))
s+='%1.3e\t'%(sig_lc2_lower)
s+='%1.3e\n'%(sig_lc2_upper)

with open(args.output_filename, "a") as myfile:
    myfile.write(s)


print s.split('\n')[0]

if(args.display>0):
	figure()
	subplot(211)
	delay_range = arange(min(delay),max(delay),0.05)
	h,b=numpy.histogram(delay, delay_range)
	#h/=numpy.sum(h)
	plot(b[:-1],h, 'k', drawstyle='steps-mid')

	plot([(delay_lower+delay_upper)/2.,(delay_lower+delay_upper)/2.],[0.,max(h)],'k-', lw=2)
	plot([delay_lower,delay_lower],[0.,max(h)],'k--', lw=2)
	plot([delay_upper,delay_upper],[0.,max(h)],'k--', lw=2)
	plot([delay_lower_95,delay_lower_95],[0.,max(h)],'k-.', lw=2)
	plot([delay_upper_95,delay_upper_95],[0.,max(h)],'k-.', lw=2)
	plot([delay_lower_98,delay_lower_98],[0.,max(h)],'k:', lw=2)
	plot([delay_upper_98,delay_upper_98],[0.,max(h)],'k:', lw=2)
	ylim(0.,max(h))
	xlabel('delay, days')
	ylabel('counts')
	subplot(212)
	low = min([min(sig_dt),min(sig_lc1), min(sig_lc2)])
	high = max([max(sig_dt),max(sig_lc1), max(sig_lc2)])
	sig_range = arange(low,high,(high-low)/100.)

	h_sdt,b_sdt=numpy.histogram(sig_dt, bins=sig_range, density = 1)
	norm_sdt = numpy.sum(h_sdt)
	#print norm_sdt
	plot(b_sdt[:-1],array(h_sdt)/norm_sdt, 'k', drawstyle='steps-mid')

	h_slc1,b_slc1=numpy.histogram(sig_lc1, bins=sig_range, density = 1)
	norm_slc1 = numpy.sum(h_slc1)
	plot(b_slc1[:-1],h_slc1/norm_slc1, 'b', drawstyle='steps-mid')

	h_slc2,b_slc2=numpy.histogram(sig_lc2, bins=sig_range, density = 1)
	norm_slc2 = numpy.sum(h_slc2)
	plot(b_slc2[:-1],h_slc2/norm_slc2, 'r', drawstyle='steps-mid')

	xlabel(r'$\sigma$, $mag/\sqrt{day}$')
	ylabel('counts')
	#show()
	'''
	plot([sig_dt_lower,sig_dt_lower],[0.,max(h)],'k--', lw=2)
	plot([(sig_dt_lower+sig_dt_upper)/2.,(sig_dt_lower+sig_dt_upper)/2.],[0.,max(h)],'k-', lw=2)
	plot([sig_dt_upper,sig_dt_upper],[0.,max(h)],'k--', lw=2)

	plot([sig_lc1_lower,sig_lc1_lower],[0.,max(h)],'b--', lw=2)
	plot([(sig_lc1_lower+sig_lc1_upper)/2.,(sig_lc1_lower+sig_lc1_upper)/2.],[0.,max(h)],'b-', lw=2)
	plot([sig_lc1_upper,sig_lc1_upper],[0.,max(h)],'b--', lw=2)


	h,b=numpy.histogram(sig_lc2, bins=50)
	h/=len(h)
	plot(b[:-1],h, 'k', drawstyle='steps-mid')
	plot([sig_lc2_lower,sig_lc2_lower],[0.,max(h)],'r--', lw=2)
	plot([(sig_lc2_lower+sig_lc2_upper)/2.,(sig_lc2_lower+sig_lc2_upper)/2.],[0.,max(h)],'r-', lw=2)
	plot([sig_lc2_upper,sig_lc2_upper],[0.,max(h)],'r--', lw=2)
	ylim(0.,max(h))
	'''
if(args.display!=0): show()


