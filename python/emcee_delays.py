# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
import time 
from hubble_observation_schedule import hubble_obs_sched
from delayed_lightcurve import delayed_lightcurve
from delay_ll import kelly_delay_ll

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
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
  if -1.e3<delay<1.e3 and -10<delta_mag<10. and 0.1<sigma<10. and 10.<tau<800. and -100.<avg_mag<100.:
        return 0.0
  return -np.inf
    
def lnprob(theta, t, lc1, err1, lc2, err2):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_delay_ll(theta, t, lc1, err1, lc2, err2)

#HUBBLE-LIKE OBSERVATION SCHEDULE 
t0 = time.clock()
delay=1.5
delta_mag=-1.

tau=121.
sigma=0.56
avg_mag=18.5
redshift = 0.658 # redshift
Nsteps=1000

t=hubble_obs_sched(80)
#t=arange(0.,1825.,5.)
lc1,lc2=delayed_lightcurve(t, delay, delta_mag, redshift,  tau, avg_mag, sigma, Nsteps)

error=0.02
err1=array([random.gauss(0.,error) for _ in xrange(len(t))])
err2=array([random.gauss(0.,error) for _ in xrange(len(t))])
#X_data=X+err
#t*=(1+z) # convert back to time at the observer
#print lc1
figure()
ax=subplot(221)
plot(t,lc1, 'bo', ms=6,  label='light curve 1')
plot(t,lc2, 'r.', ms=3, label='light curve 2')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delayed Light Curves', fontsize=18)
#xlim(0., 6.)
subplot(223, sharex=ax, sharey=ax)
plot(t,lc1,       'bo', ms=6, label='light curve 1')
plot(t-delay,lc2-delta_mag, 'r.', ms=3, label='light curve 2')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delay Aligned Light Curves', fontsize=18)
#xlim(0., 6.)

lc1+=err1
lc2+=err2

subplot(222, sharex=ax, sharey=ax)
plot(t,lc1, 'bo', ms=6,  label='2% error')
plot(t,lc2, 'r.', ms=3, label='2% error')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delayed Light Curves\nwith Photometric Errors', fontsize=18)
#xlim(0., 6.)
subplot(224, sharex=ax, sharey=ax)
plot(t,lc1,       'bo', ms=6, label='2% error')
plot(t-delay,lc2-delta_mag, 'r.', ms=3, label='2% error')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delay Aligned Light Curves\nwith Photometric Errors', fontsize=18)
#xlim(0., 6.)
#show()

'''
# NUMERICALLY OPTIMIZE THE SOLUTION
nll = lambda *args: -kelly_delay_ll_3(*args)
#result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
result = op.fmin(nll, [delay, log10(sigma), log10(tau), b], args=(t, lc1, err1, lc2, err2))
#m_ml, b_ml, lnf_ml = result["x"]
delay_ml, log10_sigma_ml, log10_tau_ml, b_ml = result
if(result[2]<0.): result[2]*=-1.

print 'delay=%1.2e\tdelay_ml=%1.2e'%(delay,delay_ml)
print 'log10_sigma=%1.2e\tlog10_sigma_ml=%1.2e'%(log10(sigma),log10_sigma_ml)
print 'log10_tau=%1.2e\tlog10_tau_ml=%1.2e'%(log10(tau),log10_tau_ml)
print 'b=%1.2e\tb_ml=%1.2e'%(b,b_ml)
'''
#RUN EMCEE
ndim, nwalkers = 5, 100
print 'pos'
true_vals=[delay, delta_mag, sigma, tau, avg_mag]
pos = [true_vals + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]
print 'sampler'
e1=error*ones(len(t))
e2=error*ones(len(t))
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, lc1, e1, lc2, e2))
print 'burn-in sampler.run_mcmc'
#sampler.run_mcmc(pos, 500)
sampler.run_mcmc(pos, 1000)
sampler.reset()
print 'sampler.run_mcmc'
#sampler.run_mcmc(pos, 100)
sampler.run_mcmc(pos, 100)
print 'sampler.chain'
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

mc_delay, mc_delta_mag, mc_sigma, mc_tau, mc_avg_mag = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

print 'Errors based on 16th, 50th, and 84th percentile'
print 'Param     \tTrue\t\tMC_rec\terr_up\terr_low\tdiff'
print 'delay     \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(delay, mc_delay[0], mc_delay[1], mc_delay[2],mc_delay[0]-delay)
print 'delta_mag \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(delta_mag, mc_delta_mag[0], mc_delta_mag[1], mc_delta_mag[2],mc_delta_mag[0]-delta_mag)
print 'sigma     \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(sigma, mc_sigma[0], mc_sigma[1], mc_sigma[2],mc_sigma[0]-sigma)
print 'tau       \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(tau, mc_tau[0], mc_tau[1], mc_tau[2],mc_tau[0]-tau)
print 'avg_mag   \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(avg_mag, mc_avg_mag[0], mc_avg_mag[1], mc_avg_mag[2],mc_avg_mag[0]-avg_mag)
#figure(2)
#for i in range(ndim):
    #figure()
    #hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
    #title("Dimension {0:d}".format(i))
fig = triangle.corner(samples, labels=["delay", "delta_mag", "$\sigma$", r"$\tau$", "avg_mag"],
                      truths=[delay, delta_mag, sigma, tau, avg_mag])
fig.savefig("triangle.png")
print time.clock() - t0, "seconds process time"
show()
