# -*- coding: utf-8 -*-
from pylab import *
import random
import numpy as np
from drw_lightcurve import drw_lightcurve
from hubble_observation_schedule import hubble_obs_sched
#from qsr_ll import hojjati_ll
from qsr_ll import kelly_ll

rcParams['font.size']=18
rcParams['legend.fontsize']=14
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.7
rcParams['figure.subplot.wspace']=0.5
rcParams['figure.subplot.left']=0.15
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.93
rcParams['figure.subplot.bottom']=0.13
rcParams['figure.figsize'] = 8, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1

#M_BH = 10**(7.78/8.) #in 10^8 M_sun
#R = 1./(1.1*M_BH)*121 # in 100 R_suns
z = 0.658 # redshift
#avg_mag=19.5
avg_mag=19.
tau=121.
sigma=0.56 #for RX J1131 at 1000 A

X_0=avg_mag #initial condition of light curve
Nsteps=1000
Nsamples=1000
t=hubble_obs_sched(80)
#t=arange(0.,1825.,5.)
X = drw_lightcurve(t, X_0, tau, sigma, avg_mag, z, Nsteps)
error=0.02
err=array([random.gauss(0.,error) for _ in xrange(len(X))])
X_data=X+err

figure()
subplot(111)
errorbar(t,X_data,yerr=error, fmt='k.', label='2% error')
plot(t,X, 'rs', ms=4, mec='r', label='Light Curve')
legend(loc=2)
xlabel('time, days')
ylabel('Rel. Mag.')
#show()

figure()
subplot(311)
avg_array=arange(avg_mag-50.,avg_mag+50.,1.)
ll_avg_array=[]
for m in avg_array:
  ll_avg_array.append(kelly_ll([sigma, tau, m], t, X_data, error*ones(len(t))))
print len(avg_array), len(ll_avg_array)
plot(avg_array,ll_avg_array, 'k-')
plot([avg_mag,avg_mag],[max(ll_avg_array)-10.,max(ll_avg_array)+1.],'k--')
ylim(max(ll_avg_array)-5.,max(ll_avg_array)+1.)
xlabel('avg magnitude, mag')
ylabel('log likelihood')

subplot(312)
tau_array=10**arange(-1.,5.,0.01)
#tau_array=arange(tau-50.,tau+50.,0.1)
ll_tau_array=[]
for m in tau_array:
  ll_tau_array.append(kelly_ll([sigma, m, avg_mag], t, X_data, error*ones(len(t))))
semilogx(tau_array,ll_tau_array, 'k-')
plot([tau,tau],[max(ll_tau_array)-10.,max(ll_tau_array)+1.],'k--')
ylim(max(ll_tau_array)-5.,max(ll_tau_array)+1.)
xlabel(r'$\tau$, day')
ylabel('log likelihood')

subplot(313)
sig_array=10**arange(-1.,1.,0.005)
ll_sig_array=[]
for m in sig_array:
  ll_sig_array.append(kelly_ll([m, tau, avg_mag], t, X_data, error*ones(len(t))))
print len(sig_array), len(ll_sig_array)
semilogx(sig_array,ll_sig_array, 'k-')
plot([sigma,sigma],[max(ll_sig_array)-10.,max(ll_sig_array)+1.],'k--')
ylim(max(ll_sig_array)-5.,max(ll_sig_array)+1.)
xlabel(r'$\sigma$, mag day$^{-1/2}$')
ylabel('log likelihood')

show()

# tau-avg_mag likelihood map
print 'COMPUTING tau-avg_mag LIKELIHOOD MAP'
tau_array=10**arange(-5.,5.,0.1)
avg_mag_array=arange(avg_mag-4.,avg_mag+4.,0.1)
print 'len(X)',len(X)
Ntau=len(tau_array)
Navg_mag=len(avg_mag_array)
Xg, Yg=np.meshgrid(array(avg_mag_array), log10(array(tau_array)))
Z = Xg*Yg
max_val=0.
print 'avg_mag',avg_mag
#for j in range(0,Navg_mag):
  #print 'avg_mag=', avg_mag_array[j]
  #ll=[]
  #for i in range(0,Ntau):
    #val = kelly_ll([sigma, tau_array[i], avg_mag_array[j]], t, X_data, err)
    #print '\t',val
    #ll.append(val)
    #Z[i][j]=val
  #max_val=max(max_val,max(ll))
  #if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(avg_mag_array[j]))
  #if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(avg_mag_array[j]))
  #ylim(-max_val,max_val)
  #xlabel('tau, days')
#legend(loc=3)
show()
figure()
print 'avg_mag',avg_mag
for i in range(0,Ntau):
  for j in range(0,Navg_mag):
   if(Z[i][j]<-max_val): Z[i][j]=-max_val
dx=Xg[1]-Xg[0]
dy=Yg[1]-Yg[0]
pcolor(Yg, Xg, Z)
cbar=colorbar()
cbar.ax.set_ylabel(r'$2\ln(p(Y|\sigma,\tau, avg_mag)$)', rotation=270)
levels=arange(-max_val,max_val,0.1*max_val)
contour(Yg, Xg, Z, levels=levels, colors='k')
plot([log10(tau)],[avg_mag],'wo', mew=2, mfc='none', mec='w')
ylim(min(avg_mag_array),max(avg_mag_array))
xlim(min(log10(tau_array)),max(log10(tau_array)))
xlabel(r'$\log_{10}(\tau/$days$)$')
ylabel(r'avg_mag magnitude')
#title('Hojjati Likelihood\nfor One Quasar Light Curve')
title('Kelly Likelihood\nfor One Quasar Light Curve')
#clabel(CS, inline=1, fontsize=10)
#title('Simplest default with labels')

#show()





# sigma-b likelihood map
print 'COMPUTING sigma-avg_mag LIKELIHOOD MAP'
sigma_array=10**arange(-3.,3.,0.1)
avg_mag_array=arange(avg_mag-4.,avg_mag+4.,0.1)
print 'len(X)',len(X)
Nsigma=len(sigma_array)
Navg_mag=len(avg_mag_array)
Xg, Yg=np.meshgrid(array(avg_mag_array), log10(array(sigma_array)))
Z = Xg*Yg
max_val=0.
print 'avg_mag',avg_mag
for j in range(0,Navg_mag):
  print 'avg_mag=', avg_mag_array[j]
  ll=[]
  for i in range(0,Nsigma):
    val = kelly_ll([sigma_array[i], tau, avg_mag_array[j]], t, X_data, err)
    ll.append(val)
    Z[i][j]=val
  max_val=max(max_val,max(ll))
  #if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(sigma_array[j]))
  #if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(sigma_array[j]))
  #ylim(-max_val,max_val)
  #xlabel('tau, days')
#legend(loc=3)
figure()
print 'avg_mag',avg_mag
for i in range(0,Nsigma):
  for j in range(0,Navg_mag):
   if(Z[i][j]<-max_val): Z[i][j]=-max_val
dx=Xg[1]-Xg[0]
dy=Yg[1]-Yg[0]
pcolor(Yg, Xg, Z)
cbar=colorbar()
cbar.ax.set_ylabel(r'$2\ln(p(Y|\sigma,\tau, avg_mag)$)', rotation=270)
#levels=arange(-max_val,max_val,0.1*max_val)
#contour(Yg, Xg, Z, levels=levels, colors='k')
plot([log10(sigma)],[avg_mag],'wo', mew=2, mfc='none', mec='w')
ylim(min(avg_mag_array),max(avg_mag_array))
xlim(min(log10(sigma_array)),max(log10(sigma_array)))
#xlabel(r'$\log_{10}(\tau/$days$)$')
ylabel(r' avg_mag magitude')
xlabel(r'$\log_{10}(\sigma/$mag day$^{-1/2})$')
#title('Hojjati Likelihood\nfor One Quasar Light Curve')
title('Kelly Likelihood\nfor One Quasar Light Curve')
#clabel(CS, inline=1, fontsize=10)
#title('Simplest default with labels')
show()



print 'COMPUTING tau-sigma LIKELIHOOD MAP'

tau_array=10**arange(-5.,5.,0.1)
sigma_array=10**arange(-3,3,0.1)
print 'len(X)',len(X)
Ntau=len(tau_array)
Nsig=len(sigma_array)
Xg, Yg=np.meshgrid(log10(array(sigma_array)), log10(array(tau_array)))
Z = Xg*Yg
#figure()
max_val=0.
for j in range(0,Nsig):
  ll=[]
  for i in range(0,Ntau):
    #val = hojjati_ll(t, X_data, err, sigma_array[j], tau_array[i], 0.)
    val = kelly_ll([sigma_array[j], tau_array[i], avg_mag], t, X_data, err)
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
ylabel(r'$\log_{10}(\sigma/$mag day$^{-1/2})$')
#title('Hojjati Likelihood\nfor One Quasar Light Curve')
title('Kelly Likelihood\nfor One Quasar Light Curve')
#clabel(CS, inline=1, fontsize=10)
#title('Simplest default with labels')
show()


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
    val = kelly_ll([sigma_array[j], tau_array[i], avg_mag], t, X_data, err)
    ll.append(val)
  max_val=max(max_val,max(ll))
  if(j<=6): semilogx(tau_array,ll, '+-',lw=2, label='$\sigma$=%1.2e'%(sigma_array[j]))
  if(j>6): semilogx(tau_array,ll, '+--', lw=3, label='$\sigma$=%1.2e'%(sigma_array[j]))
  ylim(-max_val,2.*max_val)
  xlabel('tau, days')
  ylabel('log likelihood')
  title('True sigma=%1.2e'%(sigma))
legend(loc=4)
#show()

#SIGMA DEPENDENCE
#figure()
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
    val = kelly_ll([sigma_array[j], tau_array[i], avg_mag], t, X_data, err)
    ll.append(val)
  max_val=max(max_val,max(ll))
  #print len(sigma_array), len(ll)
  if(i<=6): semilogx(sigma_array,ll, '+-',lw=2, label=r'$\tau$=%1.2e'%(tau_array[i]))
  if(i>6): semilogx(sigma_array,ll, '+--', lw=3, label=r'$\tau$=%1.2e'%(tau_array[i]))
  ylim(-max_val,2.*max_val)
  xlabel(r'sigma, mag day$^{-1/2}$')
  ylabel('log likelihood')
  title('True tau=%1.2e'%(tau))
legend(loc=3)
#show()


show()