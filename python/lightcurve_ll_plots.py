# -*- coding: utf-8 -*-
from pylab import *
import random
import numpy as np
from drw_lightcurve import drw_lightcurve
from hubble_observation_schedule import hubble_obs_sched
#from qsr_ll import hojjati_ll
from qsr_ll import kelly_ll
import time

rcParams['font.size']=18
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.7
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.15
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.88
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 10, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1

#M_BH = 10**(7.78/8.) #in 10^8 M_sun
#R = 1./(1.1*M_BH)*121 # in 100 R_suns
z = 0.658 # redshift
b=-1.0
tau=121.
sigma=0.56 #for RX J1131 at 1000 A

X_0=b*tau #initial condition of light curve
print 'tau',tau
Nsteps=1000
Nsamples=1000

t=hubble_obs_sched(80)
#figure(1)
#subplot(211)
#plot(t, ones(len(t)), '.')
#xlabel('time, days')
#title('Hubble-like obs. schedule, 80 orbits')
#ylim(0.,2.)
#xlim(0.,max(t))
#subplot(212)
#plot(t, ones(len(t)), '.')
#xlabel('time, days')
#xlim(0.,1.)
#ylim(0.,2.)
#title('Hubble-like obs. schedule, 1 day zoom in')
X = drw_lightcurve(t, X_0, tau, sigma, b, z, Nsteps)
error=0.02
err=array([random.gauss(0.,error) for _ in xrange(len(X))])
X_data=X+err
#t*=(1+z) # convert back to time at the observer

figure()
subplot(111)
errorbar(t,X_data,yerr=error, fmt='k.', label='2% error')
plot(t,X, 'rs', ms=4, mec='r', label='Light Curve')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')
#show()

#SIGMA DEPENDENCE
figure()
tau_array=tau*10**arange(-5.,5.,1.)
#sigma_array=array([0.0056, 0.056,0.56, 5.6, 56.])
sigma_array=10**arange(-3,3,0.1)
Ntau=len(tau_array)
Nsig=len(sigma_array)
max_val=0.
for i in range(0,Ntau):
  ll=[]
  for j in range(0,Nsig):
    #val = hojjati_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    val = kelly_ll([sigma_array[j], tau_array[i], b], t, X_data, err)
    ll.append(val)
  max_val=max(max_val,max(ll))
  #print len(sigma_array), len(ll)
  if(i<=6): semilogx(sigma_array,ll, '+-',lw=2, label=r'$\tau$=%1.2e'%(tau_array[i]))
  if(i>6): semilogx(sigma_array,ll, '+--', lw=3, label=r'$\tau$=%1.2e'%(tau_array[i]))
  ylim(-10.*max_val,2.*max_val)
  xlabel(r'sigma, mag day$^{-1/2}$')
  ylabel('log likelihood')
  title('True tau=%1.2e'%(tau))
legend(loc=3)
#show()


#TAU DEPENDENCE
figure()
tau_array=10**arange(-5.,5.,0.1)
#sigma_array=array([0.0056, 0.056,0.56, 5.6, 56.])
sigma_array=sigma*10**arange(-2.,2.5,0.5)
Ntau=len(tau_array)
Nsig=len(sigma_array)
max_val=0.
for j in range(0,Nsig):
  ll=[]
  for i in range(0,Ntau):
    #val = hojjati_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    val = kelly_ll([sigma_array[j], tau_array[i], b], t, X_data, err)
    ll.append(val)
  max_val=max(max_val,max(ll))
  if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(sigma_array[j]))
  if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(sigma_array[j]))
  ylim(-10.*max_val,2.*max_val)
  xlabel('tau, days')
  ylabel('log likelihood')
  title('True sigma=%1.2e'%(sigma))
legend(loc=4)
#show()

tau_array=10**arange(-5.,5.,0.1)
sigma_array=10**arange(-3,3,0.1)
print 'len(X)',len(X)
Ntau=len(tau_array)
Nsig=len(sigma_array)
Xg, Yg=np.meshgrid(log10(array(sigma_array)), log10(array(tau_array)))
Z = Xg*Yg
figure()
max_val=0.
for j in range(0,Nsig):
  ll=[]
  for i in range(0,Ntau):
    #val = hojjati_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    val = kelly_ll([sigma_array[j], tau_array[i], b], t, X_data, err)
    ll.append(val)
    Z[i][j]=val
  max_val=max(max_val,max(ll))
  #if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(sigma_array[j]))
  #if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(sigma_array[j]))
  #ylim(-max_val,max_val)
  #xlabel('tau, days')
#legend(loc=3)


#Z=(Yg-0.1*Xg)
figure()
#max_val=max(Z)
for i in range(0,Ntau):
  for j in range(0,Nsig):
   if(Z[i][j]<-10.*max_val): Z[i][j]=-10*max_val
#delta = 0.025
#Z = bivariate_normal(10**X, 10**Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
dx=Xg[1]-Xg[0]
dy=Yg[1]-Yg[0]
pcolor(Yg, Xg, Z)
cbar=colorbar()
cbar.ax.set_ylabel(r'$2\ln(p(Y|\sigma,\tau)$)', rotation=270)
levels=arange(-max_val,max_val,0.1*max_val)
contour(Yg, Xg, Z, levels=levels, colors='k')
#pcolor(Yg-dy/2, Xg-dx/2., Z)
plot([log10(tau)],[log10(sigma)],'wo', mew=2, mfc='none', mec='w')
ylim(min(log10(sigma_array)),max(log10(sigma_array)))
xlim(min(log10(tau_array)),max(log10(tau_array)))
xlabel(r'$\log_{10}(\tau/$days$)$')
ylabel(r'$\log_{10}(\sigma/$mag day$^{-1/2})$')
#title('Hojjati Likelihood\nfor One Quasar Light Curve')
title('Kelly Likelihood\nfor One Quasar Light Curve')
#clabel(CS, inline=1, fontsize=10)
#title('Simplest default with labels')



show()