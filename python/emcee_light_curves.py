# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
import scipy.optimize as op
import emcee
import triangle
import random
from drw_lightcurve import drw_lightcurve
from hubble_observation_schedule import hubble_obs_sched
from qsr_ll import kelly_ll

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

def lnprior(theta):
  sigma, tau, b = theta
  if 1.e-10<sigma<1.e10 and 1.e-10<tau<1.e10 and -1000.<b<1000.:
        return 0.0
  return -np.inf

#test lnprior
'''
for k in arange(-1.e4,1.e4, 1.e2):
  theta = [0.56, 121, k]
  print k,lnprior(theta)
exit()
'''

def lnprob(theta, t, lc, err):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + kelly_ll(theta, t, lc, err)
#M_BH = 10**(7.78/8.) #in 10^8 M_sun
#R = 1./(1.1*M_BH)*121 # in 100 R_suns

z = 0.658 # redshift
b=0.0
tau=121.
sigma=0.56 #for RX J1131 at 1000 A

X_0=b*tau #initial condition of light curve
print 'tau',tau
Nsteps=1000
Nsamples=1000

t=hubble_obs_sched(80)
X = drw_lightcurve(t, X_0, tau, sigma, b, z, Nsteps)
error=0.02
err=array([random.gauss(0.,error) for _ in xrange(len(X))])
X_data=X+err
figure()
subplot(111)
errorbar(t,X_data,yerr=error, fmt='k.', label='2% error')
plot(t,X, 'rs', ms=4, mec='r', label='Light Curve')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')

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

#RUN EMCEE
ndim, nwalkers = 3, 100
#ndim, nwalkers = 3, 10
print 'pos'
true_vals=[sigma, tau, b]
#print 1e-4*np.random.randn(ndim)
pos = [true_vals + 1e-3*np.random.randn(ndim) for i in range(nwalkers)]
#print pos
print 'sampler'
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(t, X_data, error*ones(len(X))))
print 'burn-in sampler.run_mcmc'
#sampler.run_mcmc(pos, 500)
sampler.run_mcmc(pos, 100)
sampler.reset()
print 'sampler.run_mcmc'
sampler.run_mcmc(pos, 1000)
#sampler.run_mcmc(pos, 100)
#for i in range(ndim):
    #figure()
    #hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
    #title("Dimension {0:d}".format(i))
#show()
print 'sampler.chain'
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
#exit()
#PROBLEM MUST BE HERE!
#samples[:, 2] = np.exp(samples[:, 2])
mc_sigma, mc_tau, mc_b = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84],axis=0)))

print 'Errors based on 16th, 50th, and 84th percentile'
print 'Param \tTrue\t\tMC_rec\terr_up\terr_low\tdiff'
print 'sig\t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(sigma, mc_sigma[0], mc_sigma[1], mc_sigma[2],mc_sigma[0]-sigma)
print 'tau\t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(tau, mc_tau[0], mc_tau[1], mc_tau[2],mc_tau[0]-tau)
print 'b     \t%1.2e\t%1.2e\t+%1.2e\t-%1.2e\t%1.2e'%(b, mc_b[0], mc_b[1], mc_b[2],mc_b[0]-b)
#figure(2)
#for i in range(ndim):
    #figure()
    #hist(sampler.flatchain[:,i], 100, color="k", histtype="step")
    #title("Dimension {0:d}".format(i))
fig = triangle.corner(samples, labels=["$\sigma$", r"$\tau$", "$b$"],
                      truths=[sigma, tau, b])
fig.savefig("triangle.png")
show()
