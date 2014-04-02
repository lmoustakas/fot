# -*- coding: utf-8 -*-
from pylab import *
import random

from drw_lightcurve import drw_lightcurve
from hubble_observation_schedule import hubble_obs_sched

rcParams['font.size']=14
rcParams['legend.fontsize']=14
rcParams['legend.borderaxespad']=0.
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.7
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.1
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.9
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 8, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1


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
figure(1)
subplot(211)
plot(t, ones(len(t)), '.')
xlabel('time, days')
title('Hubble-like obs. schedule, 80 orbits')
ylim(0.,2.)
xlim(0.,max(t))
subplot(212)
plot(t, ones(len(t)), '.')
xlabel('time, days')
xlim(0.,1.)
ylim(0.,2.)
title('Hubble-like obs. schedule, 1 day zoom in')
X = drw_lightcurve(t, X_0, tau, sigma, b, z, Nsteps)
figure(2)
ax=subplot(411)
plot(t,X, 'ro', ms=3, mec='r', label='Light Curve')
#errorbar(t,X, yerr=0.015*abs(X), fmt=',', ms=3)
val_max=max(X)
val_min=min(X)
if(val_max>0.): val_max*=3.
if(val_max<0.): val_max*=0.3
if(val_min>0.): val_min*=0.5
if(val_min<0.): val_min*=1.5
ylim(val_min,val_max)
legend(loc=1)
title('Simulated Observations\n$\sigma=%1.2f$, $\\tau=%1.0f$ days'%(sigma,tau))
xlabel('time, days')
ylabel('Rel. Mag.')
ylim(val_min,val_max)
subplot(412, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.015) for _ in xrange(len(X))])
errorbar(t,X+err,yerr=0.015, fmt='k.', label='1.5% error')
plot(t,X,'ro', ms=3, mec='r')
legend(loc=1)
xlabel('time, days')
ylabel('Rel. Mag.')
ylim(val_min,val_max)
subplot(413, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.03) for _ in xrange(len(X))])
errorbar(t,X+err,yerr=0.03, fmt='k.', label='3% error')
plot(t,X,'ro', ms=3, mec='r')
legend(loc=1)
xlabel('time, days')
ylabel('Rel. Mag.')
ylim(val_min,val_max)
subplot(414, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.1) for _ in xrange(len(X))])
errorbar(t,X+err ,yerr=0.1 , fmt='k.', label='10% error')
plot(t,X,'ro', ms=3, mec='r')
xlabel('time, days')
ylabel('Rel. Mag.')
legend(loc=1)
#title('Rel. Mag. vs. Time\n for tau = 1e-2, 1e-1, 1e0, 1e1, 1e2\n from top to bottom\n X(0)=0, b=0, sigma=1, 1000 integration steps ')
ylim(val_min,val_max)
show()