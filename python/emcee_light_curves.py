# -*- coding: utf-8 -*-
#!/usr/bin/env python
from pylab import *
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
from drw_lightcurve import drw_lightcurve
from hubble_observation_schedule import hubble_obs_sched
from qsr_ll import kelly_ll
from qsr_ll    import kelly_estimates
from read_data import read_data

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.5
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.1
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.88
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 10, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1

# DEFINE THE LIKELIHOOD FUNCTION

#def lnprior(theta):
  #sigma, tau, avg_mag = theta
  #if 0.01<sigma<10.0 and 10.<tau<800 and -100<avg_mag<100:
        #return 0.0
  #return -np.inf
'''
  if 1.e-10<sigma<1.e10 and 1.e-10<tau<1.e10 and -1000.<b<1000.:
        return 0.0
  return -np.inf
''' 

#test lnprior
'''
for k in arange(-1.e4,1.e4, 1.e2):
  theta = [0.56, 121, k]
  print k,lnprior(theta)
exit()
'''

#M_BH = 10**(7.78/8.) #in 10^8 M_sun
#R = 1./(1.1*M_BH)*121 # in 100 R_suns

#z = 0.658 # redshift
#avg_mag=19.5
#tau=121.
#sigma=0.56 #for RX J1131 at 1000 A

#X_0=avg_mag #initial condition of light curve
##print 'tau',tau
#Nsteps=1000
#Nsamples=1000

#t=hubble_obs_sched(80)
#X = drw_lightcurve(t, X_0, tau, sigma, avg_mag, z, Nsteps)
#error=0.02
#err=array([random.gauss(0.,error) for _ in xrange(len(X))])
#X_data=X+err
#figure()
#subplot(111)
#errorbar(t,X_data,yerr=error, fmt='k.', label='2% error')
#plot(t,X, 'rs', ms=4, mec='r', label='Light Curve')
#legend(loc=2)
#xlabel('time, days')
#ylabel('Rel. Mag.')

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


'''
#test lnprob
#b_array=concatenate([-10**arange(5.,-5.,-1.),10**arange(-5.,5.,1.)])
b_array=arange(-100.,100.)
figure()
for b_val in b_array:
  theta = [sigma, tau, b_val]
  print b_val,lnprob(theta, t, X_data, error*ones(len(X_data)))
  plot([b_val],[lnprob(theta, t, X_data, error*ones(len(X_data)))] ,'k.')
show()
exit()
'''

'''
# commenting out the sanity check code block 
# NUMERICALLY OPTIMIZE THE SOLUTION
print 'NUMERICAL OPTIMIZATION OF THE SOLUTION'
nll = lambda *args: -kelly_ll(*args)
#result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
result = op.fmin(nll, [sigma, tau, b], args=(t, X_data, error*ones(len(X_data))))
#m_ml, b_ml, lnf_ml = result["x"]
sigma_ml, tau_ml, b_ml = result
#if(result[2]<0.): result[2]*=-1.

print 'sigma=%1.2e\tsigma_ml=%1.2e'%((sigma),sigma_ml)
print 'tau=%1.2e\ttau_ml=%1.2e'%((tau),tau_ml)
print 'b=%1.2e\tb_ml=%1.2e'%(b,b_ml)
#exit()
'''
def emcee_delay_estimator(t, lc, err, output_tag): 
  t0 = datetime.datetime.now()
  print 'Process Started on', t0

  crude_sig, crude_tau, crude_avg_mag = crude_lc_param_estimate(t,lc)

  print ''
  print 'CRUDE LIGHT CURVE PARAMETERS'
  print '----------------------------'
  print 'parameter value'
  print '----------------------------'
  print 'sigma     %+1.2f '%(crude_sig)
  print 'tau       %+1.2f '%(crude_tau)
  print 'avg_mag   %+1.2f '%(crude_avg_mag)
  print ''

  print ''
  print 'BEGINNING LIGHT CURVE PARAMETER OPTIMIZATION'
  print ''
  op_sig, op_tau, op_avg_mag = optimize_lc_param_estimate(t,lc,err)

  op_sig=abs(op_sig)
  print ''
  print 'OPTIMIZED LIGHT CURVE PARAMETERS'
  print '----------------------------'
  print 'parameter value'
  print '----------------------------'
  print 'sigma     %+1.2f '%(op_sig)
  print 'tau       %+1.2f '%(op_tau)
  print 'avg_mag   %+1.2f '%(op_avg_mag)
  print ''
  print 'Process Started on', t0
  print 'It is now         ', datetime.datetime.now()
  print ''
  print '-----------------'
  print 'RUNNING EMCEE'
  print '-----------------'
  sigma   = op_sig
  tau     = op_tau
  avg_mag = op_avg_mag
  def lnprior(_theta):
    _sigma, _tau, _avg_mag = _theta
    if sigma*0.1<_sigma<sigma*10. and tau*0.1<_tau<tau*10. and avg_mag-3.<_avg_mag<avg_mag+3.:
	  return 0.0
    return -np.inf
	
  def lnprob(theta, t, lc, err):
      lp = lnprior(theta)
      if not np.isfinite(lp):
	  return -np.inf
      return lp + kelly_ll(theta, t, lc, err)

  ndim, nwalkers = 3, 100
  n_burn_in_iterations = 100
  n_iterations = 1000

  print 'pos'
  true_vals=[op_sig, op_tau, op_avg_mag]
  pos = [true_vals + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]
  print 'sampler'
  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, lc, err))
  print 'burn-in sampler.run_mcmc'
  sampler.run_mcmc(pos, n_burn_in_iterations)
  
  samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
  print("\tMean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_sigma, mc_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))
  print '\tErrors based on 16th, 50th, and 84th percentile'
  print '\tParam     \tMC_rec\terr_up\terr_low'
  print '\tsigma     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_sigma[0], mc_sigma[1], mc_sigma[2])
  print '\ttau       \t%1.2e\t+%1.2e\t-%1.2e'%(mc_tau[0], mc_tau[1], mc_tau[2])
  print '\tavg_mag   \t%1.2e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])

  #fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
			#truths=[mc_delay[0], mc_delta_mag[0], mc_sigma[0], mc_tau[0], mc_avg_mag[0]])
  fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"])
  fig.savefig("burn_in_triangle_%s.png"%(output_tag))

  
  sampler.reset()
  print ''
  print 'Process Started on', t0
  print 'It is now         ', datetime.datetime.now()
  print ''
  print 'sampler.run_mcmc'
  sampler.run_mcmc(pos, n_iterations)
  print 'sampler.chain'
  samples = sampler.chain[:, int(0.1*n_iterations):, :].reshape((-1, ndim))
  print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
  mc_sigma, mc_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

  print 'Errors based on 16th, 50th, and 84th percentile'
  print 'Param   \tMC_rec\terr_up\terr_low'
  print 'sig     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_sigma[0],   mc_sigma[1],   mc_sigma[2])
  print 'tau     \t%1.2e\t+%1.2e\t-%1.2e'%(mc_tau[0],     mc_tau[1],     mc_tau[2])
  print 'avg_mag \t%1.2e\t+%1.2e\t-%1.2e'%(mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2])

  fig = triangle.corner(samples, labels=["$\sigma$", r"$\tau$", "avg_mag"])
  fig.savefig("triangle_%s.png"%(output_tag))
  #show()
  print 'Process Started on', t0
  print 'It ended on        ', datetime.datetime.now()

  
t, mag1, e1, mag2, e2, mag3, e3, mag4, e4 = read_data('../data/cosmograil/RXJ1131_Tewes2013.rdb') 
e1=array([max(0.1,x) for x in e1])
e2=array([max(0.1,x) for x in e2])
e3=array([max(0.1,x) for x in e3])
e4=array([max(0.1,x) for x in e4])
t   = t[0:70]
lc  = mag1[0:70]
err = e1[0:70]
#lc2=mag4
#err2=e4
emcee_delay_estimator(t, lc, err, 'RXJ1131_curve_A')
#emcee_delay_estimator(t, mag2, e2, mag3, e3, 'RXJ1131_curves_BC')
