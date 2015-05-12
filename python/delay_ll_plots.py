# -*- coding: utf-8 -*-
from pylab import *
import numpy as np
import random
from delayed_lightcurve import delayed_lightcurve
from hubble_observation_schedule import hubble_obs_sched
from delay_ll import kelly_delay_ll

rcParams['font.size']=16
rcParams['legend.fontsize']=14
rcParams['legend.borderpad']=0.2
rcParams['legend.labelspacing']= 0.1  
rcParams['legend.handlelength']= 1.  
rcParams['legend.handleheight']= 0.1 
rcParams['legend.handletextpad']= 0.1 
rcParams['legend.borderaxespad']= 0.0 
rcParams['legend.columnspacing']= 1.  
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.5
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

delay=1.5
delta_mag=-1.

tau=121.
sigma=0.56
avg_mag=18.5
redshift = 0.658 # redshift
Nsteps=1000

t=hubble_obs_sched(80)
#t=arange(0.,1825.,5.)
lc1,lc2=delayed_lightcurve(t, delay, delta_mag, redshift,  tau, avg_mag, sigma, Nsteps)

error=avg_mag*.02
err1=array([random.gauss(0.,error) for _ in xrange(len(t))])
err2=array([random.gauss(0.,error) for _ in xrange(len(t))])
#X_data=X+err
#t*=(1+z) # convert back to time at the observer
#print lc1
figure()
ax=subplot(221)
plot(t,lc1, 'bo', ms=6,  label='light curve 1')
plot(t,lc2, 'r.', ms=3, label='light curve 2')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delayed Light Curves', fontsize=18)
#xlim(0., 6.)
subplot(223, sharex=ax, sharey=ax)
plot(t,lc1,       'bo', ms=6, label='light curve 1')
plot(t-delay,lc2-delta_mag, 'r.', ms=3, label='light curve 2')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delay Aligned Light Curves', fontsize=18)
#xlim(0., 6.)

lc1+=err1
lc2+=err2

subplot(222, sharex=ax, sharey=ax)
plot(t,lc1, 'bo', ms=6,  label='2% error')
plot(t,lc2, 'r.', ms=3, label='2% error')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delayed Light Curves\nwith Photometric Errors', fontsize=18)
#xlim(0., 6.)
subplot(224, sharex=ax, sharey=ax)
plot(t,lc1,       'bo', ms=6, label='2% error')
plot(t-delay,lc2-delta_mag, 'r.', ms=3, label='2% error')
legend(loc=1)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delay Aligned Light Curves\nwith Photometric Errors', fontsize=18)
#xlim(0., 6.)
#show()
'''
figure(3)
delay_array=arange(0.,3.,0.001)
ll=[]
for i in range(0,len(delay_array)):
  val = kelly_delay_ll(t, lc1, err1, lc2, err2, delay_array[i], sigma, tau, 0.)
  ll.append(val)

plot(delay_array,ll, '.-', ms=3)
max_ll=max(ll)
plot([delay, delay], [0.1*max_ll,1.1*max_ll], 'k--')
ylim(0.1*max_ll,1.1*max_ll)
xlabel('delay, days')
ylabel('Likelihood')
'''

e1=error*ones(len(t))
e2=error*ones(len(t))

figure()
subplot(321)
plot(t,lc1,       'bo', ms=6, label='2% error')
plot(t,lc2, 'r.', ms=3, label='2% error')
legend(loc=0)
xlabel('time, days')
ylabel('flux, mag')
title('Simulated Delayed Light Curves\nwith Photometric Errors', fontsize=18)
#xlim(0., 6.)
ylim(min(min(lc1),min(lc2))-0.1, max(max(lc1),max(lc2))+2.)

subplot(322)
delay_array=arange(delay-0.5,delay+0.5,0.005)
ll_delay_array=[]
for m in delay_array:
  ll_delay_array.append(kelly_delay_ll([m, delta_mag, sigma, tau, avg_mag], t, lc1, e1, lc2, e2))
print len(delay_array), len(ll_delay_array)
plot(delay_array,ll_delay_array, 'k-')
plot([delay,delay],[max(ll_delay_array)-10.,max(ll_delay_array)+1.],'k--')
ylim(max(ll_delay_array)-10.,max(ll_delay_array)+1.)
xlabel('delay, days')
ylabel('log likelihood')

subplot(323)
delta_mag_array=arange(delta_mag-0.5,delta_mag+0.5,0.005)
ll_delta_mag_array=[]
for m in delta_mag_array:
  ll_delta_mag_array.append(kelly_delay_ll([delay, m, sigma, tau, avg_mag], t, lc1, e1, lc2, e2))
print len(delta_mag_array), len(ll_delta_mag_array)
plot(delta_mag_array,ll_delta_mag_array, 'k-')
plot([delta_mag,delta_mag],[max(ll_delta_mag_array)-10.,max(ll_delta_mag_array)+1.],'k--')
ylim(max(ll_delta_mag_array)-10.,max(ll_delta_mag_array)+1.)
xlabel('delta_mag, mag')
ylabel('log likelihood')

subplot(324)
sig_array=arange(0.,1.,0.005)
ll_sig_array=[]
for m in sig_array:
  ll_sig_array.append(kelly_delay_ll([delay, delta_mag, m, tau, avg_mag], t, lc1, e1, lc2, e2))
print len(sig_array), len(ll_sig_array)
plot(sig_array,ll_sig_array, 'k-')
plot([sigma,sigma],[max(ll_sig_array)-10.,max(ll_sig_array)+1.],'k--')
ylim(max(ll_sig_array)-5.,max(ll_sig_array)+1.)
xlabel(r'$\sigma$, mag day$^{-1/2}$')
ylabel('log likelihood')

subplot(325)
#tau_array=10**arange(-1.,5.,0.01)
tau_array=arange(1.e-3,2*tau,tau/1000.)
ll_tau_array=[]
for m in tau_array:
  ll_tau_array.append(kelly_delay_ll([delay, delta_mag, sigma, m, avg_mag], t, lc1, e1, lc2, e2))
plot(tau_array,ll_tau_array, 'k-')
plot([tau,tau],[max(ll_tau_array)-10.,max(ll_tau_array)+1.],'k--')
ylim(max(ll_tau_array)-5.,max(ll_tau_array)+1.)
xlabel(r'$\tau$, day')
ylabel('log likelihood')

subplot(326)
#tau_array=10**arange(-1.,5.,0.01)
avg_mag_array=arange(avg_mag-50.,avg_mag+50.,1.)
ll_avg_mag_array=[]
for m in avg_mag_array:
  ll_avg_mag_array.append(kelly_delay_ll([delay, delta_mag, sigma, tau, m], t, lc1, e1, lc2, e2))
plot(avg_mag_array,ll_avg_mag_array, 'k-')
plot([avg_mag,avg_mag],[max(ll_avg_mag_array)-10.,max(ll_avg_mag_array)+1.],'k--')
ylim(max(ll_avg_mag_array)-5.,max(ll_avg_mag_array)+1.)
xlabel(r'avg_mag, mag')
ylabel('log likelihood')


show()






Ntrials=1000
diff_val=[]
for n in range(0,Ntrials):
  #print 'Trial',n
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
    t_sample.append(t_val)
    if(mod(t_val,90.)<t_obs*T_orbit):
      t.append(t_val)
    if(t_val>80.*90.): bail=1
  t=array(t)/(60.*24.)

  delta_mag=0.
  delay=1.2+random.gauss(0.,0.1)
  lc1,lc2=delayed_lightcurve(t, delay, delta_mag, tau, b, sigma, 1000)
  error=0.015
  err1=array([random.gauss(0.,error) for _ in xrange(len(t))])
  err2=array([random.gauss(0.,error) for _ in xrange(len(t))])
  lc1+=err1
  lc2+=err2

  delay_array=arange(delay-1.,delay+1.,0.001)
  #delay_array=arange(0.,3.,0.001)
  ll=[]
  for i in range(0,len(delay_array)):
    val = kelly_delay_ll(	)
    ll.append(val)

  k_max=argmax(ll)
  print delay_array[k_max]-delay
  sys.stdout.flush()

  diff_val.append(delay_array[k_max]-delay)

figure(4)
hist(diff_val)
diff_val=array(diff_val)
print cumsum(diff_val)[len(diff_val)-1]/len(diff_val)
print sqrt(cumsum(diff_val**2)[len(diff_val)-1]/len(diff_val)-(cumsum(diff_val)[len(diff_val)-1]/len(diff_val))**2)








show()
tau_array=10**arange(-5.,5.,0.5)
sigma_array=10**arange(-3,3,0.1)
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