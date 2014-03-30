# -*- coding: utf-8 -*-
from pylab import *
import numpy
from drw_lightcurve import X_int

def delayed_lightcurve(time_array, delay, delta_mag, tau, b, sigma, Nsteps):
    t1=time_array.copy()
    t2=time_array.copy() + delay
    id1=zeros(len(t1))
    id2=ones(len(t1))
    
    t_cat=concatenate([t1,t2])
    id_cat=concatenate([id1,id2])
    
    #time order the array
    t_cat_sort  =  [x for (x,y) in sorted(zip(t_cat,id_cat))]
    id_cat_sort =  [y for (x,y) in sorted(zip(t_cat,id_cat))]
    

    lc=[0.]
    for k in range(1,len(t_cat_sort)):
      lc.append(X_int(t_cat_sort[k]-t_cat_sort[k-1],lc[k-1], tau, b, sigma, Nsteps))
    
    lc1 = [y for (x,y) in zip(id_cat_sort,lc) if x==0.]
    lc2 = [y for (x,y) in zip(id_cat_sort,lc) if x==1.]
    
    
    #subplot(311)
    #plot(t_cat_sort,lc,'o-', lw=1, ms=6)
    #subplot(312)
    #plot(t1,lc1,'g.-')
    #plot(t2,lc2,'r.-')
    #subplot(313)
    #plot(t1,lc1,'g.-')
    #plot(t1,lc2,'r.-')
    
    return lc1,lc2
    #show()
    