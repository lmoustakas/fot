#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
"Full of time" project, developing a the analysis tools for
a) generating realistic quasar light curves (with uncertainties)
b) calculating quasar light curve structural parameters through Bayesian inference
c) generating time delayed and flux-offset quasar light curves (with uncertainties)
d) calculating the time delays and light curve structure parameters through Bayesian inference

This is the central library of functions. 

Let us set the environment variable FOTDIR to the local location of the repository. 
'''
import matplotlib
matplotlib.use('Agg') # Note, this MUST be before importing pylab or matplotlib.pyplot
from pylab import *
import os
from astropy.io import ascii
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
import datetime 

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['legend.borderpad']=0.1
rcParams['legend.borderaxespad']=0.
rcParams['axes.formatter.limits']=-4,7
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.5
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.15
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.88
rcParams['figure.subplot.bottom']=0.15
rcParams['figure.figsize'] = 10, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1
rcParams['xtick.major.pad'] = 1
rcParams['ytick.major.pad'] = 1

def read_cosmograil_data(fname,mags,magerrs):
    '''
    Read in data from the cosmograil release.
    Basic example: 
    time,m,me=read_cosmograil_data('RXJ1131_Tewes2013.rdb',['mag_A','mag_B'],['magerr_A','magerr_B'])
    Then, under the returned "m", for example, the magnitudes can be accessed via m['mag_A'], and so on. 
    '''
    data=ascii.read(os.environ['FOTDIR']+'/data/cosmograil/'+fname,data_start=2)
    return data['mhjd'],data[mags],data[magerrs]

def crude_lc_param_estimate(_t, _x):
  avg_lc = cumsum(_x)[len(_x)-1]/len(_x)
  tau=mean(diff(_t)*(avg_lc-_x[0:len(_x)-1]))/mean(diff(_x))
  #tau=-diff(_t)/log(1.-(diff(_x)/(avg_lc-_x[0:len(_x)-1])))
  #print 'tau',tau
  tau=abs(tau)
  if(tau==0.): tau=100.
  #print tau
  avg_dt=cumsum(diff(_t))[len(_t)-2]/(len(_t)-1.)
  var=cumsum((_x-avg_lc)**2)[len(_x)-1]/len(_x)
  sig_sq=mean((_x[1:len(_x)]-avg_lc)**2)/mean(1.-exp(-2.*(diff(_t))/tau))*2./tau
  sigma=sqrt(sig_sq)
  return sigma,tau, avg_lc

def optimize_lc_param_estimate(_t, _x, _e):
  nll = lambda *args: -kelly_ll(*args)
  result = op.fmin(nll, crude_lc_param_estimate(_t, _x), args=(_t, _x, _e))
  return result

def diff_check(_sigma, _tau, _avg_mag,_t,_x,_e):
  x_hat,dif_err = kelly_diff([_sigma, _tau, _avg_mag],_t,_x,_e)
  chi=((_x-_avg_mag)-x_hat)/dif_err
  return mean(chi), cumsum(chi**2)[len(chi)-1]/len(chi)-mean(chi)**2

def kelly_estimates(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'qsr_ll.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_diff ( [sigma, tau, b], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, tau, avg_mag = theta
  
  x=flux_array
  t=time_array
  #ph_err_array*=0.
  x_star=x-avg_mag
  #print 'in kelly_diff: avg_mag=',avg_mag 
  x_hat=[0.]
  Omega=[tau*sig**2/2.]
  a=[0.]
  ll=0.
  for i in range(1,len(t)+1):
    if(i==1): a.append(0.)
    if(i>1): a.append(np.exp(-abs(t[i-1]-t[i-2])/tau)) #t array start at 0, the rest start at 1
    Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)))
    x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)*(x_star[i-2]-x_hat[i-1]))
  #print 'mean a %1.2f'%(np.mean(a[1:len(t)+1]))
  #print 'mean Omega %1.2e'%(np.mean(Omega[1:len(t)+1]))
  return x_hat[1:len(t)+1], np.sqrt(Omega[1:len(t)+1]+ph_err_array**2)

################################################################################################
#from pylab import * #import for testing only.
def kelly_ll(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'qsr_ll.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_ll ( [sigma, tau, b], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, tau, avg_mag = theta
  if(tau<=0.): return -np.inf
  x=flux_array
  t=time_array
  x_hat, err = kelly_estimates(theta, time_array, flux_array, ph_err_array)
  x_star=x-avg_mag
  ll = np.cumsum( -((x_hat-x_star)**2 / (err**2)) - np.log(2*np.pi*(err**2)))  #ph_err_array*=0.
  return ll[len(x)-1]

def hojjati_ll(time_array, flux_array, ph_err_array, sig, tau, b):
  Nsamp=len(flux_array)
  # QUASAR DAMPED RANDOM WALK MATRIX
  M=(np.identity(Nsamp))
  print sig, tau
  for a in range(0,Nsamp):
    for b in range(0,a):
      #M[a][b]=0.
      M[a][b]=exp(-abs(time_array[a]-time_array[b])/tau)
      M[b][a]=M[a][b]	
  M*=sig**2    
  # ADD PHOTOMETRIC ERRORS
  for a in range(0,Nsamp):
    M[a][a]+=ph_err_array[a]**2
  #use slogdet
  M_det = np.linalg.det(M)
  log_M_det=cumsum(log(np.linalg.eigvals(M)))[Nsamp-1]
  if(M_det!=0.): log_M_det = log(M_det)
  M_inv = np.linalg.inv(M)
  ll=-array(matrix(flux_array[0:Nsamp])*matrix(M_inv)*matrix(flux_array[0:Nsamp]).T)[0][0]-log_M_det-Nsamp*log(2.*np.pi)
  return ll

def merge(_t, _x1, _e1, _x2, _e2, _dt, _dmag):
  _t2  = _t.copy()  - _dt 
  _x2  = _x2.copy() - _dmag 
  
  t_cat=np.concatenate([_t,_t2])
  x_cat=np.concatenate([_x1,_x2])
  e_cat=np.concatenate([_e1,_e2])
    
  #time order the array
  t_cat_sort  =  np.array([x for (x,y,z) in sorted(zip(t_cat,x_cat,e_cat))])
  x_cat_sort  =  np.array([y for (x,y,z) in sorted(zip(t_cat,x_cat, e_cat))])
  e_cat_sort  =  np.array([z for (x,y,z) in sorted(zip(t_cat,x_cat, e_cat))])
  
  return  t_cat_sort, x_cat_sort, e_cat_sort

def kelly_delay_ll(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  delay, delta_mag, sig, tau, avg_mag = theta
  t, x, e = merge(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, delta_mag)
  return kelly_ll([sig,tau,avg_mag], t, x, e)

def emcee_delay_estimator(t, lc1, e1, lc2, e2, output_tag):
  outputdir=os.environ['FOTDIR']+'/outputs/'
  t0 = datetime.datetime.now()
  print 'Process Started on', t0
  date_string = t0.strftime("_%Y_%m_%d_%H:%M")
  #print date_string
  output_tag = ''.join([output_tag,date_string ])
  print 'The output_tag for this run is:', output_tag

  crude_sig1, crude_tau1, crude_avg_lc1 = crude_lc_param_estimate(t,lc1)
  crude_sig2, crude_tau2, crude_avg_lc2 = crude_lc_param_estimate(t,lc2)

  print ''
  print 'CRUDE LIGHT CURVE PARAMETERS'
  print '----------------------------'
  print 'parameter lc1   \tlc2'
  print '----------------------------'
  print 'sigma     %+1.2f \t%+1.2f'%(crude_sig1,    crude_sig2)
  print 'tau       %+1.2f \t%+1.2f'%(crude_tau1,    crude_tau2)
  print 'avg_mag   %+1.2f \t%+1.2f'%(crude_avg_lc1, crude_avg_lc2)
  print ''

  print ''
  print 'BEGINNING LIGHT CURVE PARAMETER OPTIMIZATION'
  print ''
  op_sig1, op_tau1, op_avg_lc1 = optimize_lc_param_estimate(t,lc1,e1)
  op_sig2, op_tau2, op_avg_lc2 = optimize_lc_param_estimate(t,lc2,e2)

  op_sig1=abs(op_sig1)
  op_sig2=abs(op_sig2)
  print ''
  print 'OPTIMIZED LIGHT CURVE PARAMETERS'
  print '--------------------------------'
  print 'parameter lc1   \tlc2'
  print '--------------------------------'
  print 'sigma     %+1.2f \t%+1.2f'%(op_sig1,    op_sig2)
  print 'tau       %+1.2f \t%+1.2f'%(op_tau1,    op_tau2)
  print 'avg_mag   %+1.2f \t%+1.2f'%(op_avg_lc1, op_avg_lc2)
  print ''

  print '\nProcess Started on', t0
  print 'It is now         ', datetime.datetime.now()
  #ESTIMATE DELAY
  print ''
  print 'time series is %1.2f days long'%(t[len(t)-1]-t[0])
  print 'with  %d data points'%(len(t))
  nbins=int(len(t)/50.)
  fig=figure()
  counts, bin_vals, somemore =hist(log10(diff(t)), nbins)
  xlabel('Time Between Measurments, log10(days)')
  fig.savefig(outputdir+"delta_t_histogram_%s.png"%(output_tag))
  
  # ESTIMATE THE MODE OF THE DISTRIBUTION
  n_max=argmax(counts)
  mode=bin_vals[n_max]
  print 'The most common sample time distance is %1.2f days within the range of [%1.2f,%1.2f] days'%(10**mode, 10**(0.5*(bin_vals[n_max-1]+bin_vals[n_max])), 10**(0.5*(bin_vals[n_max]+bin_vals[n_max+1])))
  print ''
  
  #delta_delay=(10**mode/50.)
  #print 'Scanning for initial delay estimate with delay step %1.2f'%(delta_delay)
  #del_array=arange(-1000.,1000.1,delta_delay)
  #del_array=arange(-1000.,1000.1,10.)
  data_record_length = t[len(t)-1]-t[0]
  del_array=arange(-0.5*data_record_length,+0.5*data_record_length,0.1)
  count=-1
  ll=[]
  max_val  = -1.e10
  c_sig    = 0.5*(op_sig1+op_sig2)
  c_tau    = 0.5*(op_tau1+op_tau2)
  c_avg_lc = (op_avg_lc1)
  delta_mag=(op_avg_lc2-op_avg_lc1)
  for delta in del_array:
    count+=1
    ct, cx, ce = merge(t, lc1, e1, lc2, e2, delta, delta_mag)
    val = kelly_ll([c_sig, c_tau, c_avg_lc] , ct, cx, ce)
    if(val>max_val):
      max_val=val
    if(count==int(0.01*len(del_array))): print '1% done'
    if(count==int(0.10*len(del_array))): print '10% done'
    if(count==int(0.25*len(del_array))): print '25% done'
    if(count==int(0.50*len(del_array))): print '50% done'
    if(count==int(0.75*len(del_array))): print '75% done'
    ll.append(val)
  print 'DONE'
  print '\nProcess Started on', t0
  print 'It is now         ', datetime.datetime.now()
  print ''
  #PLOT DELAY SEARCH
  fig=figure()
  subplot(311)
  plot(del_array,ll,'b-')
  xlabel('delay, days')
  ylabel('likelihood')
  title('Delay Likelihood Scan\nsigma=%1.2f, tau=%1.2f, avg_mag=%1.2f'%(c_sig, c_tau, c_avg_lc))
  k_max=argmax(ll)
  plot([del_array[k_max],del_array[k_max]],[min(ll),max(ll)],'k--')
  print 'BEST DELAY=%1.5e'%(del_array[k_max])
  subplot(312)
  delta_delay=del_array[1]-del_array[0]
  k_lo = k_max - int(50/(delta_delay)+1)
  k_up = k_max + int(50/(delta_delay)+1)
  plot(del_array[k_lo:k_up],ll[k_lo:k_up],'b-')
  xlabel('delay, days')
  ylabel('likelihood')
  #title('Delay Likelihood Scan\nsigma=%1.2f, tau=%1.2f, avg_mag=%1.2f'%(c_sig, c_tau, c_avg_lc))
  plot([del_array[k_max],del_array[k_max]],[min(ll),max(ll)],'k--')
  ylim(min(ll[k_lo:k_up]),ll[k_max]+100.)
  subplot(313)
  delta_delay=del_array[1]-del_array[0]
  k_lo = k_max - int(5/(delta_delay)+1)
  k_up = k_max + int(5/(delta_delay)+1)
  plot(del_array[k_lo:k_up],ll[k_lo:k_up],'b-')
  xlabel('delay, days')
  ylabel('likelihood')
  #title('Delay Likelihood Scan\nsigma=%1.2f, tau=%1.2f, avg_mag=%1.2f'%(c_sig, c_tau, c_avg_lc))
  plot([del_array[k_max],del_array[k_max]],[min(ll),max(ll)],'k--')
  ylim(max(ll)-1000.,max(ll))
  fig.savefig(outputdir+"delay_search_%s.png"%(output_tag))
  delay = del_array[k_max]


  #ESTIMATES OF QSR LIGHT CURVE
  delta_mag = (op_avg_lc2-op_avg_lc1)
  sigma = abs(c_sig)
  tau = c_tau
  avg_mag=op_avg_lc1
  print ''
  print 'INITIAL PARAMETERS (average of lc1 and lc2)'
  print '-------------------------------------------'
  print 'delay      = %1.2f'%(delay)
  print 'delta_mag  = %1.2f'%(delta_mag)
  print 'avg_mag    = %1.2f'%(avg_mag)
  print 'tau        = %1.2f'%(tau)
  print 'sigma      = %1.2f'%(sigma)
  print ''

  ct, cx, ce = merge(t, lc1, e1, lc2, e2, delay, delta_mag)
  c_sig, c_tau, c_avg_lc = optimize_lc_param_estimate(ct,cx,ce)
  sigma   = abs(c_sig)
  tau     = c_tau
  avg_mag = c_avg_lc
  print ''
  print 'INITIAL PARAMETERS (optimization of merged light curve)'
  print '-------------------------------------------------------'
  print 'delay      = %1.2f'%(delay)
  print 'delta_mag  = %1.2f'%(delta_mag)
  print 'avg_mag    = %1.2f'%(avg_mag)
  print 'tau        = %1.2f'%(tau)
  print 'sigma      = %1.2f'%(sigma)
  print ''
  #exit()

  #PLOT LIGHT CURVES
  fig = figure()
  ax=subplot(311)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t,lc2,e2, fmt='r.', ms=3, label='light curve 2')
  legend(loc=1)
  xlabel('time, day')
  ylabel('magnitude')
  #plot(ct,cx,'k.', ms=3)
  ax=subplot(312)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t-delay, lc2-delta_mag, e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
 #plot(ct,cx,'k.', ms=3)
  ax=subplot(313)
  ct, cx, ce = merge(t, lc1, e1, lc2, e2, delay, delta_mag)
  errorbar(ct,cx, ce,fmt='b.', ms=3, label = 'merged light curve')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
  fig.savefig(outputdir+"op_light_curves_%s.png"%(output_tag))
  #show()

  def lnprior(_theta):
    _delay, _delta_mag, _sigma, _tau, _avg_mag = _theta
    if delay-0.5*data_record_length<_delay<delay+0.5*data_record_length and delta_mag-1.<_delta_mag<delta_mag+1. and sigma*0.1<_sigma<sigma*10. and tau*0.1<_tau<tau*10. and avg_mag-1.<_avg_mag<avg_mag+1.:
	  return 0.0
    return -np.inf
      
  def lnprob(theta, t, lc1, err1, lc2, err2):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_delay_ll(theta, t, lc1, err1, lc2, err2)

  
  print ''
  print '-------------'
  print 'RUNNING EMCEE'
  print '-------------'
  print ''
  #RUN EMCEE
  ndim, nwalkers = 5, 100
  n_burn_in_iterations = 100
  n_iterations = 1000
  print 'pos'
  true_vals=[delay, delta_mag, sigma, tau, avg_mag]
  #print 'np.random.randn(ndim)',np.random.randn(ndim)
  pos = [true_vals + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
  #pos = [[uniform(delay-100,delay+100.),uniform(delta_mag-1.,delta_mag+1.), uniform(0.,sigma*10.), uniform(tau*0.1, tau*10.), uniform(avg_mag-0.2+avg_mag-0.2)] for i in range(nwalkers)]
  #pos = [true_vals + [uniform(-1.e3, 1.e3), uniform(-10.,10.), uniform(0.1,10.), uniform(0.1,800.), uniform(-50.,50.)] for i in range(nwalkers)]
  #uniform([low, high, size])
  r=np.random.randn(ndim)
  #pos = [true_vals + [10**mode*r[0], 0.1*delta_mag*r[1], abs((sigma*r[2]+sigma)), ] for i in range(nwalkers)]
  print 'sampler'
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, lc1, e1, lc2, e2))
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'burn-in sampler.run_mcmc'
  #sampler.run_mcmc(pos, 500)
  sampler.run_mcmc(pos, n_burn_in_iterations)

  samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_delay, mc_delta_mag, mc_sigma, mc_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  print '\tErrors based on 16th, 50th, and 84th percentile'
  print '\tParam     \tMC_rec\terr_up\terr_low'
  print '\tdelay     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_delay[0], mc_delay[1], mc_delay[2])
  print '\tdelta_mag \t%1.2e\t+%1.2e\t-%1.2e'%(mc_delta_mag[0], mc_delta_mag[1], mc_delta_mag[2])
  print '\tsigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print '\ttau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_tau[0], mc_tau[1], mc_tau[2])
  print '\tavg_mag   \t%1.2e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])

  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  fig.savefig(outputdir+"burn_in_triangle_%s.png"%(output_tag))

  sampler.reset()
  
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.run_mcmc'
  #sampler.run_mcmc(pos, 100)
  sampler.run_mcmc(pos, n_iterations)
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.chain'
  samples = sampler.chain[:, int(0.1*n_iterations):, :].reshape((-1, ndim))
  print len(samples)
  print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

  mc_delay, mc_delta_mag, mc_sigma, mc_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  print 'Errors based on 16th, 50th, and 84th percentile'
  print 'Param     \tMC_rec\terr_up\terr_low'
  print 'delay     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_delay[0], mc_delay[1], mc_delay[2])
  print 'delta_mag \t%1.2e\t+%1.2e\t-%1.2e'%(mc_delta_mag[0], mc_delta_mag[1], mc_delta_mag[2])
  print 'sigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print 'tau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_tau[0], mc_tau[1], mc_tau[2])
  print 'avg_mag   \t%1.2e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])
  #figure(2)
  #for i in range(ndim):
      #figure()
      #hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
      #title("Dimension {0:d}".format(i))
  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  fig.savefig(outputdir+"triangle_%s.png"%(output_tag))
  
  fig = figure()
  ax=subplot(311)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t,lc2,e2, fmt='r.', ms=3, label='light curve 2')
  legend(loc=1)
  xlabel('time, day')
  ylabel('magnitude')
  #plot(ct,cx,'k.', ms=3)
  ax=subplot(312)
  errorbar(t,lc1,e1, fmt='b.', ms=3, label='light curve 1')
  errorbar(t-mc_delay[0], lc2-mc_delta_mag[0], e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
 #plot(ct,cx,'k.', ms=3)
  ax=subplot(313)
  ct, cx, ce = merge(t, lc1, e1, lc2, e2, mc_delay[0], mc_delta_mag[0])
  errorbar(ct,cx, ce,fmt='b.', ms=3, label = 'merged light curve')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
  fig.savefig(outputdir+"mc_light_curves_%s.png"%(output_tag))

  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  #show()
#################################
#END def emcee_delay_estimator(...)

#########################################
#SIMULATION FUNCTIONS
#########################################

def hubble_obs_sched(num_orbits):
  t_sample=[0.] #initialize time sampling array
  t=[0.] #initialize observation time array
  T_orbit=90. #in minutes
  frac_obs=0.6   #fraction of orbit where source is visible
  t_int=20. # integration time in minutes
  time_uncertainty = 5. #uncertainty on when the satellite is ready to observe in minutes
  min_per_day=1./(60.*24) #convert from minutes to days
  bail=0 #flag for bailing out of the time sampling loop
  while (bail==0):
    t_val=t_sample[len(t_sample)-1]+t_int+(random.gauss(0.,time_uncertainty)) #in minutes
    t_sample.append(t_val)
    if(np.mod(t_val,T_orbit)<frac_obs*T_orbit):
      t.append(t_val)
    if(t_val>num_orbits*T_orbit): bail=1
  t=np.array(t)*min_per_day
  return t

def X_int(delta_t,X_0, tau, avg_mag, sigma, Nsteps):
  #print delta_t, tau
  num_steps=10*int(1+delta_t/tau)
  if(num_steps>Nsteps):
    num_steps=Nsteps
    print 'X_int WARNING: accurate statistics requires Nsteps>1000 steps delta_t=%1.2e, tau=%1.2e'%(delta_t, tau)
  ds = delta_t/Nsteps
  s=arange(Nsteps, dtype=float64)*ds
  dB=array([random.gauss(0.,sqrt(ds)) for _ in xrange(Nsteps)])
  #print '%1.2e\t%1.2e\t%1.2e'%(cumsum(dB)[Nsteps-1]/Nsteps, sqrt(cumsum(dB**2)[Nsteps-1]/Nsteps-(cumsum(dB)[Nsteps-1]/Nsteps)**2), ds)
  x1 = exp(-delta_t/tau)*X_0
  x2 = avg_mag*(1.-exp(-delta_t/tau))
  x3 = sigma*cumsum( exp(-(delta_t-s)/tau)*dB )[Nsteps-1]
  return x1+x2+x3
  
def drw_lightcurve(time_array, X_0, tau, sigma, avg_mag, redshift, Nsteps): 
  # convert parameters to rest frame
  tau_rest=tau*1./(1.+redshift)
  sigma_rest=sigma*sqrt(1.+redshift)
  avg_mag_rest=avg_mag
  t_rest=time_array*1./(1.+redshift)
  X=[X_0]
  for k in range(1,len(time_array)):
    X.append(X_int(t_rest[k]-t_rest[k-1], X[k-1], tau_rest, avg_mag_rest, sigma_rest, Nsteps))
  X=array(X)
  return X 
  
def delayed_lightcurve(time_array, delay, delta_mag, redshift, tau, avg_mag, sigma, Nsteps):
    t1=time_array.copy()
    t2=time_array.copy() - delay
    id1=zeros(len(t1))
    id2=ones(len(t1))
    
    t_cat=concatenate([t1,t2])
    id_cat=concatenate([id1,id2])
    
    #time order the array
    t_cat_sort  =  array([x for (x,y) in sorted(zip(t_cat,id_cat))])
    id_cat_sort =  array([y for (x,y) in sorted(zip(t_cat,id_cat))])
    
    X_0 = random.gauss( avg_mag, sigma*sqrt(tau/2.) )
    lc=drw_lightcurve(t_cat_sort, X_0, tau, sigma, avg_mag, redshift, Nsteps)
    
    lc1 = array([y for (x,y) in zip(id_cat_sort,lc) if x==0.])
    lc2 = array([y for (x,y) in zip(id_cat_sort,lc) if x==1.])
    
    lc2+=delta_mag
    return lc1,lc2
