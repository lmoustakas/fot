# -*- coding: utf-8 -*-
from pylab import *
from delayed_lightcurve import delayed_lightcurve
from hubble_observation_schedule import hubble_obs_sched
import random

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.2
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.1
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.9
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 8, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 7
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1



delay=1.5
delta_mag=1.

tau=121.
sigma=0.56
avg_mag=18.5
redshift = 0.658 # redshift
Nsteps=1000

time_array=hubble_obs_sched(80)
lc1,lc2=delayed_lightcurve(time_array, delay, delta_mag, redshift,  tau, avg_mag, sigma, Nsteps)

#figure(1)
#plot(time_array,'.')

figure(2)
#subplot(311)
#plot(concatenate(t,lc,'o-', lw=1, ms=6)
subplot(211)
plot(time_array,lc1,'b.-', label='light curve 1')
plot(time_array,lc2,'r.-', label='light curve 2')
legend(loc=1)
ylabel('flux, arb. u.')
xlabel('time, days')
subplot(212)
plot(time_array,lc1,'b.-', label='light curve 1')
plot(time_array-delay,lc2-delta_mag,'r.-', label='delay shifted\nlight curve 2')
legend(loc=1)
ylabel('flux, arb. u.')
xlabel('time, days')
#subplot(313)
#plot(t1,lc1,'g.-')
#plot(t1,lc2,'r.-')


show()