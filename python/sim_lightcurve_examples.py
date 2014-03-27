# -*- coding: utf-8 -*-
from pylab import *
from drw_lightcurve import *
import time

rcParams['font.size']=14
rcParams['legend.fontsize']=14
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


M_BH = 10**(7.78/8.) #in 10^8 M_sun
R = 1./(1.1*M_BH)*121 # in 100 R_suns
z = 2. # redshift

X_0=0. #initial condition of light curve
tau=1.1 * M_BH * R
print 'tau',tau
b=0.0
sigma=0.56 #for RX J1131 at 1000 A
Nsteps=1000
Nsamples=1000

X=[0.]
t_sample=[0.]
t=[0.]
T_orbit=90. #minutes
t_obs=0.6   #fraction of orbit where source is visible
t_int=20. # minutes
dt=90./(60.*24)
counter=0
bail=0
while (bail==0):
  counter+=1
  t_val=t_sample[len(t_sample)-1]+20.+(random.gauss(0.,5.)) #in minutes
  #print t_val, mod(t_val,90.), t_obs*T_orbit
  t_sample.append(t_val)
  if(mod(t_val,90.)<t_obs*T_orbit):
    t.append(t_val)
  if(t_val>80.*90.): bail=1
t=array(t)/(60.*24.)

figure(1)
plot(t, ones(len(t)), '.')
plot(array(t_sample)/(60.*24.), 2.*ones(len(t_sample)), '.')
xlabel('time, days')
ylim(0.,3.)
#show()
#exit() 
for k in range(1,len(t)):
  #t.append(t[k-1]+1.+abs(random.gauss(0.,dt)))
  #t.append(t[k-1]+dt)
  X.append(X_int(t[k]-t[k-1], X[k-1], tau, b, sigma, Nsteps))
X=array(X)
#X+=10.
figure(2)
ax=subplot(411)
plot(t,X, 'ro', ms=3, mec='r', label='Light Curve')
#errorbar(t,X, yerr=0.015*abs(X), fmt=',', ms=3)
legend(loc=2)
title('Simulated Observations\n$\sigma=%1.2f$, $\\tau=%1.0f$ days'%(sigma,tau))
xlabel('time, days')
ylabel('Rel. Mag.')
subplot(412, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.015) for _ in xrange(len(X))])
errorbar(t,X+err,yerr=0.015, fmt='k.', label='1.5% error')
plot(t,X,'ro', ms=3, mec='r')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')
subplot(413, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.03) for _ in xrange(len(X))])
errorbar(t,X+err,yerr=0.03, fmt='k.', label='3% error')
plot(t,X,'ro', ms=3, mec='r')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')
subplot(414, sharex=ax, sharey=ax)
err=array([random.gauss(0.,0.1) for _ in xrange(len(X))])
errorbar(t,X+err ,yerr=0.1 , fmt='k.', label='10% error')
plot(t,X,'ro', ms=3, mec='r')
xlabel('time, days')
ylabel('Rel. Mag.')
legend(loc=2)
ylim(-1.,1.)
#title('Rel. Mag. vs. Time\n for tau = 1e-2, 1e-1, 1e0, 1e1, 1e2\n from top to bottom\n X(0)=0, b=0, sigma=1, 1000 integration steps ')


show()