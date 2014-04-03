# -*- coding: utf-8 -*-
from pylab import *
import numpy as np

def kelly_delay_ll(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  delay, delta_mag, sig, tau, avg_mag = theta
  #see Equations 6-12 of Kelly, 2009
  t1=time_array.copy()
  x1=flux_array1.copy()
  t2=time_array.copy()  - delay
  x2=flux_array2.copy() - delta_mag 
  
  t_cat=concatenate([t1,t2])
  x_cat=concatenate([x1,x2])
  err_cat=concatenate([ph_err_array1,ph_err_array2])
    
  #time order the array
  t_cat_sort  =  array([x for (x,y,z) in sorted(zip(t_cat,x_cat,err_cat))])
  x_cat_sort  =  array([y for (x,y,z) in sorted(zip(t_cat,x_cat, err_cat))])
  err_cat_sort  =  array([z for (x,y,z) in sorted(zip(t_cat,x_cat, err_cat))])
  
  x=array(x_cat_sort)
  t=array(t_cat_sort)
  ph_err_array=array(err_cat_sort)
  
  x_star=x-avg_mag
  x_hat=[0.]
  Omega=[tau*sig**2/2.]
  a=[0.]
  ll=0.
  for i in range(1,len(t)):
    #a.append(exp(-abs(t[i]-t[i-1])/tau))
    #Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)))
    #x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-1]**2)*(x_star[i-1]-x_hat[i-1]))
    #ll += -(x_hat[i] - x_star[i])**2/(Omega[i]+ph_err_array[i]**2) - log(2*pi*(Omega[i]+ph_err_array[i]**2))
    if(i==1): a.append(0.)
    if(i>1): a.append(np.exp(-abs(t[i-1]-t[i-2])/tau)) #t array start at 0, the rest start at 1
    Omega.append(Omega[0]*(1-a[i]**2) + a[i]**2*Omega[i-1]*(1.-Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)))
    x_hat.append(a[i]*x_hat[i-1] + a[i]*Omega[i-1]/(Omega[i-1]+ph_err_array[i-2]**2)*(x_star[i-2]-x_hat[i-1]))
    ll += -(x_hat[i] - x_star[i-1])**2/(Omega[i]+ph_err_array[i-1]**2) - np.log(2*np.pi*(Omega[i]+ph_err_array[i-1]**2))
    
  return ll
  
#def kelly_delay_ll_2(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  #delay, sig, tau, b = theta
  #return kelly_delay_ll(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, sig, tau, b)

#def kelly_delay_ll_3(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  #delay, log_sig, log_tau, b = theta
  #sig=10**log_sig
  #tau=10**log_tau
  #return kelly_delay_ll(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, sig, tau, b)

  