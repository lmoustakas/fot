# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
import datetime 
from delay_ll  import merge
from delay_ll  import kelly_delay_ll
from qsr_ll    import kelly_ll
from qsr_ll    import kelly_estimates
from read_data import read_data

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

# DEFINE THE LIKELIHOOD FUNCTION

def lnprior(theta):
  delay, delta_mag, sigma, tau, avg_mag = theta
  if -1.e3<delay<1.e3 and -10<delta_mag<10. and 0.1<sigma<10. and 0.<tau<800. and -50.<avg_mag<50.:
        return 0.0
  return -np.inf
    
def lnprob(theta, t, lc1, err1, lc2, err2):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_delay_ll(theta, t, lc1, err1, lc2, err2)

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

def emcee_delay_estimator(t, lc1, e1, lc2, e2, output_tag): 
  t0 = datetime.datetime.now()
  print 'Process Started on', t0

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
  fig.savefig("delta_t_histogram_%s.png"%(output_tag))
  
  # ESTIMATE THE MODE OF THE DISTRIBUTION
  n_max=argmax(counts)
  mode=bin_vals[n_max]
  print 'The most common sample time distance is %1.2f days within the range of [%1.2f,%1.2f] days'%(10**mode, 10**(0.5*(bin_vals[n_max-1]+bin_vals[n_max])), 10**(0.5*(bin_vals[n_max]+bin_vals[n_max+1])))
  print ''
  
  delta_delay=(10**mode/50.)
  print 'Scanning for initial delay estimate with delay step %1.2f'%(delta_delay)
  #del_array=arange(-1000.,1000.1,delta_delay)
  del_array=arange(-1000.,1000.1,10.)
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
  fig.savefig("delay_search_%s.png"%(output_tag))
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
  errorbar(t,lc1,e1, fmt='g.', ms=3, label='light curve 1')
  errorbar(t-delay, lc2-delta_mag, e2, fmt='r.', ms=3, label='del. and mag. shift lc2')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
 #plot(ct,cx,'k.', ms=3)
  ax=subplot(313)
  ct, cx, ce = merge(t, lc1, e1, lc2, e2, delay, delta_mag)
  errorbar(ct,cx, ce,fmt='k.', ms=3, label = 'merged light curve')
  xlabel('time, day')
  ylabel('magnitude')
  legend(loc=1)
  fig.savefig("light_curves_%s.png"%(output_tag))
  #show()

  print ''
  print '-------------'
  print 'RUNNING EMCEE'
  print '-------------'
  print ''
  #RUN EMCEE
  ndim, nwalkers = 5, 100
  n_burn_in_interations = 100
  n_interations = 1000
  print 'pos'
  true_vals=[delay, delta_mag, sigma, tau, avg_mag]
  print 'np.random.randn(ndim)',np.random.randn(ndim)
  #pos = [true_vals*(ones(ndim) + 1e-1*np.random.randn(ndim)) for i in range(nwalkers)]
  pos = [[uniform(delay-10,delay+10.),uniform(delta_mag-1.,delta_mag+1.), uniform(0.,sigma*10.), uniform(tau*0.1, tau*10.), uniform(avg_mag-0.2+avg_mag-0.2)] for i in range(nwalkers)]
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
  sampler.run_mcmc(pos, n_burn_in_interations)

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
  fig.savefig("burn_in_triangle_%s.png"%(output_tag))


  sampler.reset()
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.run_mcmc'
  #sampler.run_mcmc(pos, 100)
  sampler.run_mcmc(pos, n_interations)
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  print 'sampler.chain'
  samples = sampler.chain[:, 500:, :].reshape((-1, ndim))
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
  fig.savefig("triangle_%s.png"%(output_tag))
  print 'Process Started on', t0
  print 'It is currently   ', datetime.datetime.now()
  #show()

#################################
#END def emcee_delay_estimator(...)

t, mag1, e1, mag2, e2, mag3, e3, mag4, e4 = read_data('../data/cosmograil/RXJ1131_Tewes2013.rdb') 
e1=array([max(0.1,x) for x in e1])
e2=array([max(0.1,x) for x in e2])
e3=array([max(0.1,x) for x in e3])
e4=array([max(0.1,x) for x in e4])
lc1=mag2
err1=e2
lc2=mag4
err2=e4
emcee_delay_estimator(t, lc1, err1, lc2, err2, 'RXJ1131_curves_BD')
#emcee_delay_estimator(t, mag2, e2, mag3, e3, 'RXJ1131_curves_BC')
