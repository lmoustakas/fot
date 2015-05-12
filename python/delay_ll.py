# -*- coding: utf-8 -*-
import numpy as np
from qsr_ll import kelly_ll

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

###################################################################################################

def kelly_delay_ll(theta, time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2):
  delay, delta_mag, sig, tau, avg_mag = theta
  t, x, e = merge(time_array, flux_array1, ph_err_array1, flux_array2, ph_err_array2, delay, delta_mag)
  return kelly_ll([sig,tau,avg_mag], t, x, e)
  
###################################################################################################

