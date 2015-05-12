import numpy as np
import random

def hubble_obs_sched(num_orbits):
  t_sample=[0.] #initialize time sampling array
  t=[0.] #initialize observation time array
  T_orbit=90. #in minutes
  frac_obs=0.6   #fraction of orbit where source is visible
  t_int=20. # integration time in minutes
  time_uncertainty = 5. #uncertainty on when the satellite is ready to observe in minutes
  min_per_day=1./(60.*24) #convert from minutes to days
  bail=0 #flag for bailing out of the time sampling loop
  while (bail==0):
    t_val=t_sample[len(t_sample)-1]+t_int+(random.gauss(0.,time_uncertainty)) #in minutes
    t_sample.append(t_val)
    if(np.mod(t_val,T_orbit)<frac_obs*T_orbit):
      t.append(t_val)
    if(t_val>num_orbits*T_orbit): bail=1
  t=np.array(t)*min_per_day
  return t