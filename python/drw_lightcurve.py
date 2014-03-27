# -*- coding: utf-8 -*-
from pylab import *
import random

#This follows the stochastic integral equation A3 of Kelley 2009, arXiv:0903.5315v1
def X_int(delta_t,X_0, tau, b, sigma, Nsteps):
  ds = delta_t/Nsteps
  s=arange(Nsteps, dtype=float64)/(Nsteps-1)*delta_t
  #print t, s[len(s)-1]
  dB=array([random.gauss(0.,ds) for _ in xrange(Nsteps)])
  #print '%1.2e\t%1.2e\t%1.2e'%(cumsum(dB)[Nsteps-1]/Nsteps, sqrt(cumsum(dB**2)[Nsteps-1]/Nsteps-(cumsum(dB)[Nsteps-1]/Nsteps)**2), ds)
  x1=exp(-delta_t/tau)*X_0
  x2=b*tau*(1-exp(-delta_t/tau))
  x3=sigma*cumsum(exp(-(delta_t-s)/tau)*dB)[Nsteps-1]*sqrt(Nsteps)
  return x1+x2+x3
  
#This follows the stochastic differential equation 1 of Kelley 2009, arXiv:0903.5315v1
# this one is slower than the integral implementation by a factor of 7
def X_diff(t, X_0, tau, b, sigma, Nsteps):
  ds = t/Nsteps
  eps=array([random.gauss(0.,1.) for _ in xrange(Nsteps)])
  X=X_0
  for i in range(0,Nsteps):
    dX = -(1./tau)*X*ds + sigma*sqrt(ds)*eps[i] + b*ds
    X+=dX
  return X
  

