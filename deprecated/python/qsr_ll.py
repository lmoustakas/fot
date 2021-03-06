# -*- coding: utf-8 -*-
import numpy as np

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
  #print len(x_hat), len(err), len(x_star)
  '''
  #THESE LINES ARE FOR TESTING PURPOSES ONLY
  ####################################################
  from pylab import *
  figure()
  ax=subplot(221)
  errorbar(time_array,x_star,ph_err_array,fmt='b.')
  errorbar(time_array,x_hat, err, fmt='r.')
  subplot(222, sharex=ax)
  errorbar(time_array,x_star-x_hat,err,fmt='b.')
  subplot(223, sharex=ax)
  plot(time_array,(x_star-x_hat)/err,'b.')
  subplot(224)
  hist((x_star-x_hat)/err, bins=100, range=[-30.,30.])
  ####################################################
  '''
  ll = np.cumsum( -((x_hat-x_star)**2 / (err**2)) - np.log(2*np.pi*(err**2)))  #ph_err_array*=0.
  return ll[len(x)-1]
  #x_hat=[0.]
  #Omega=[tau*sig**2/2.]
  #a=[0.]
  #ll=0.
  #for i in range(1,len(t)+1):
    #if(i==1): a.append(0.)
    #if(i>1): a.append(np.exp(-abs(t[i-1]-t[i-2])/tau)) #t array start at 0, the rest start at 1
    #Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)))
    #x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)*(x_star[i-2]-x_hat[i-1]))
    #ll += -(x_hat[i] - x_star[i-1])**2/(Omega[i]+ph_err_array[i-1]**2) - np.log(2*np.pi*(Omega[i]+ph_err_array[i-1]**2))
  #return ll
  
'''
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
  #ll=-(x_hat[0] - x_star[0])**2/(Omega[0]+ph_err_array[0]**2) - np.log(2*np.pi*(Omega[0]+ph_err_array[0]**2))
  ll=0.
  #print sig, tau
  for i in range(0,len(t)):
    x_star[i]=x[i]-b*tau
    a.append(np.exp(-abs(t[i]-t[i-1])/tau))
    Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)))
    x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)*(x_star[i-1]-x_hat[i-1]))
    ll += -(x_hat[i] - x_star[i])**2/(Omega[i]+ph_err_array[i]**2) - np.log(2*np.pi*(Omega[i]+ph_err_array[i]**2))
  return ll
'''
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
  
