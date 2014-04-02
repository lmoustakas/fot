# -*- coding: utf-8 -*-
#from pylab import *
import numpy as np

################################################################################################

def kelly_ll(theta, time_array, flux_array, ph_err_array):
  #see Equations 6-12 of Kelly, 2009
  if(len(theta)!=3):
    print 'qsr_ll.py'
    print 'USAGE:***********************************************************************' 
    print 'kelly_ll ( [sigma, tau, b], time_array, flux_array, measurement_error_array )'
    print 'exiting'
    exit()
  sig, tau, b = theta
  x=flux_array
  t=time_array
  ph_err_array*=0.
  x_star=x-b*tau
  x_hat=[0.]
  Omega=[tau*sig**2/2.]
  a=[0.]
  ll=-(x_hat[0] - x_star[0])**2/(Omega[0]+ph_err_array[0]**2) - np.log(2*np.pi*(Omega[0]+ph_err_array[0]**2))
  #print sig, tau
  for i in range(1,len(t)):
    x_star[i]=x[i]-b*tau
    a.append(np.exp(-abs(t[i]-t[i-1])/tau))
    Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)))
    x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)*(x_star[i-1]-x_hat[i-1]))
    ll += -(x_hat[i] - x_star[i])**2/(Omega[i]+ph_err_array[i]**2) - np.log(2*np.pi*(Omega[i]+ph_err_array[i]**2))
  return ll

################################################################################################

def hojjati_ll(time_array, flux_array, ph_err_array, sig, tau, b):
  #ll=-999.
  Nsamp=len(flux_array)
  #Nsamp=10
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
  #print '\t%1.2e\t%1.2e'%(M_det,log_M_det)
  #print np.linalg.eigvals(M)
  #print cumprod(np.linalg.eigvals(M))
  #print cumsum(log(np.linalg.eigvals(M)))[Nsamp-1]
  #print ''
  M_inv = np.linalg.inv(M)
  #if(M_det!=0. and M_det==M_det):
  ll=-array(matrix(flux_array[0:Nsamp])*matrix(M_inv)*matrix(flux_array[0:Nsamp]).T)[0][0]-log_M_det-Nsamp*log(2.*np.pi)
    #ll=-array(matrix(flux_array[0:Nsamp])*matrix(M_inv)*matrix(flux_array[0:Nsamp]).T)[0][0]
    #ll=-log(M_det)-Nsamp*log(2.*pi)
  return ll
  
