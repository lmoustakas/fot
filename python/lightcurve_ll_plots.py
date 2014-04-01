# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
from drw_lightcurve import *
from qsr_ll import hojjati_ll
from qsr_ll import kelly_ll
import time

rcParams['font.size']=18
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.7
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.1
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.88
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 10, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1

M_BH = 10**(7.78/8.) #in 10^8 M_sun
R = 1./(1.1*M_BH)*121 # in 100 R_suns
z = 0.658 # redshift

X_0=0. #initial condition of light curve
#tau=1.1 * M_BH * R
tau=121.
#tau=1.e-5
print 'tau',tau
b=0.0
sigma=0.56 #for RX J1131 at 1000 A
#sigma=1.e-3 #for RX J1131 at 1000 A
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

#t/=(1+z) #time at the rest frame of the quasar
for k in range(1,len(t)):
  #t.append(t[k-1]+1.+abs(random.gauss(0.,dt)))
  #t.append(t[k-1]+dt)
  X.append(X_int(t[k]-t[k-1], X[k-1], tau, b, sigma, Nsteps))
X=array(X)
error=0.015
err=array([random.gauss(0.,error) for _ in xrange(len(X))])
X_data=X+err
#t*=(1+z) # convert back to time at the observer

figure(2)
subplot(111)
errorbar(t,X_data,yerr=error, fmt='k.', label='1.5% error')
plot(t,X, 'rs', ms=4, mec='r', label='Light Curve')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')

t_avg=[]
for k in range(0,len(t)-1):
  t_avg.append(0.5*(t[k+1]+t[k]))
figure(3)
subplot(211)
plot(t_avg,diff(X_data), '.')
plot(t_avg,diff(X), '.')
subplot(212)
hist(diff(X_data), bins=25)
hist(diff(X), bins=25)

#M = np.identity(3)
#tau_array=arange(1.,300.,3.)
tau_array=10**arange(-5.,5.,0.1)
#sigma_array=arange(0.1,1.0,0.1)
#sigma_array=10**arange(-3,3,0.5)
sigma_array=10**arange(-3,3,0.1)
#Nsamp=len(X)
#Nsamp=60
print 'len(X)',len(X)
Ntau=len(tau_array)
Nsig=len(sigma_array)
#Xg, Yg=np.meshgrid(sigma_array, tau_array)
Xg, Yg=np.meshgrid(log10(array(sigma_array)), log10(array(tau_array)))
Z = Xg*Yg
#M = 0.*np.random.rand(Ntau*Nsig, Nsamp,Nsamp)
figure(4)
max_val=0.

for j in range(0,Nsig):
  ll=[]
  for i in range(0,Ntau):
    #val = hojjati_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    val = kelly_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    ll.append(val)
    Z[i][j]=val
  max_val=max(max_val,max(ll))
  if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(sigma_array[j]))
  if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(sigma_array[j]))
  ylim(-max_val,max_val)
  xlabel('tau, days')
legend(loc=3)

#Z=(Yg-0.1*Xg)
figure(5)
#max_val=max(Z)
for i in range(0,Ntau):
  for j in range(0,Nsig):
   if(Z[i][j]<-max_val): Z[i][j]=-max_val
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
ylabel(r'$\log_{10}(\sigma/$arb. u.$)$')
#title('Hojjati Likelihood\nfor One Quasar Light Curve')
title('Kelly Likelihood\nfor One Quasar Light Curve')
#clabel(CS, inline=1, fontsize=10)
#title('Simplest default with labels')

show()