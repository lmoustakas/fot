# -*- coding: utf-8 -*-
from pylab import *
import random
from drw_lightcurve import X_int
from drw_lightcurve import drw_lightcurve
#import time

rcParams['font.size']=12
rcParams['legend.fontsize']=12
rcParams['axes.formatter.limits']=-4,4
rcParams['figure.facecolor'] = 'white'
rcParams['figure.subplot.hspace']=0.9
rcParams['figure.subplot.wspace']=0.8
rcParams['figure.subplot.left']=0.15
rcParams['figure.subplot.right']=0.9
rcParams['figure.subplot.top']=0.9
rcParams['figure.subplot.bottom']=0.1
rcParams['figure.figsize'] = 8, 8 
rcParams['legend.numpoints'] = 1 
rcParams['lines.markersize'] = 10 
rcParams['lines.linewidth'] = 1
rcParams['lines.markeredgewidth'] = 1

# CHECK KELLY'S RELATIONSHIPS ON THE MEAN (Kelly, 2009, Eq. A4
print 'Checking generator mean and variance relationships (Kelly, 2009 Equations A4 and A5)'

def check_sigma_dependence():
  figure()
  print 'CHECKING sigma dependence'
  b=3.
  X_s=2.
  delta_t=0.014
  tau=121.
  sigma_array=10**arange(-5.,6.,1.)
  Nsteps=100
  Ntrials=100
  figure(1)
  for i in range(0,len(sigma_array)):
    s=0.
    s2=0.
    for j in range(0,Ntrials):
      val=X_int(delta_t,X_s, tau, b, sigma_array[i], Nsteps)
      s+=val
      s2+=val**2
    mean=s/Ntrials
    var=s2/Ntrials-mean**2
    print i,'of', len(sigma_array)
    print '\tmean %1.2e %1.2e'%(mean, exp(-delta_t/tau)*X_s+b*tau*(1.-exp(-delta_t/tau)))
    print '\tvar %1.2e %1.2e'%(var, (tau*sigma_array[i]**2)/2.*(1.-exp(-2.*delta_t/tau)))
    subplot(211)
    errorbar(sigma_array[i],[mean],yerr=[sqrt(var)/sqrt(Ntrials)], fmt='ko-', ms=3)
    subplot(212)
    errorbar(sigma_array[i],[var], yerr=[var*sqrt(2/(Ntrials-1))], fmt='ko-', ms=3)

  ax1=subplot(211)
  ax1.set_xscale('log')
  #ax1.set_yscale('log')
  sigma_array=10**arange(-5.,6.,0.1)
  plot(sigma_array,exp(-delta_t/tau)*X_s+b*tau*(1.-exp(-delta_t/tau))*ones(len(sigma_array)))
  xlabel(r'sigma, mag/day$^{1/2}$')
  ylabel(r'$E(X(t)|X(s))$')
  xlim(3.e-6, 3.e5)
  ylim(-1000.,1000.)
  ax1=subplot(212)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  plot(sigma_array,(tau*sigma_array**2)/2.*(1-exp(-2.*delta_t/tau)))
  xlabel(r'sigma, mag/day$^{1/2}$')
  ylabel(r'$Var(X(t)|X(s))$')
  xlim(3.e-6, 3.e5)
  #ylim(1.e-6,1.e3)
  #ylim(1.e-5,1.e2)

  show()


def check_tau_dependence():
  figure()
  print 'CHECKING tau dependence'
  b=3.
  sigma=0.56
  X_s=2.
  delta_t=1.e2
  tau_array=10**arange(-5.,6.,1.)
  Nsteps=1000
  Ntrials=1000
  figure(1)
  for i in range(0,len(tau_array)):
    s=0.
    s2=0.
    for j in range(0,Ntrials):
      val=X_int(delta_t,X_s, tau_array[i], b, sigma, Nsteps)
      s+=val
      s2+=val**2
    mean=s/Ntrials
    var=s2/Ntrials-mean**2
    print i,'of', len(tau_array)
    print '\ttau=%1.2e  delta_t=%1.2e  delta_t/tau = %1.2e'%(tau_array[i], delta_t, delta_t/tau_array[i])
    print '\tmean %1.2e %1.2e'%(mean, exp(-delta_t/tau_array[i])*X_s+b*tau_array[i]*(1.-exp(-delta_t/tau_array[i] )))
    print '\tvar %1.2e %1.2e'%(var, (tau_array[i]*sigma**2)/2.*(1.-exp(-2.*delta_t/tau_array[i])))
    subplot(211)
    errorbar(tau_array[i],[mean],yerr=[sqrt(var)/sqrt(Ntrials)], fmt='k.-')
    subplot(212)
    errorbar(tau_array[i],[var], yerr=[var*sqrt(2/(Ntrials-1))], fmt='k.-')

  ax1=subplot(211)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  tau_array=10**arange(-5.,6.,0.1)
  plot(tau_array,exp(-delta_t/tau_array)*X_s+b*tau_array*(1.-exp(-delta_t/tau_array)),'k--', lw=2, label='delta_t=100 days')
  #plot(tau_array,b*tau_array*(1.-exp(-delta_t/tau_array)))
  #plot(tau_array,exp(-delta_t/tau_array)*X_s)
  xlabel('tau, days')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  ax1=subplot(212)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  plot(tau_array,(tau_array*sigma**2)/2.*(1-exp(-2.*delta_t/tau_array)),'k--', lw=2, label='delta_t=100 days')
  #plot(tau_array,(delta_t*sigma**2)*ones(len(tau_array)))
  #plot(tau_array,(tau_array*sigma**2)/2.)
  xlabel('tau, days')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  #ylim(1.e-5,1.e2)




  b=3.
  sigma=0.56
  X_s=2.
  delta_t=1.
  tau_array=10**arange(-5.,6.,1.)
  #Nsteps=1000
  #Ntrials=1000
  figure(1)
  for i in range(0,len(tau_array)):
    s=0.
    s2=0.
    for j in range(0,Ntrials):
      val=X_int(delta_t,X_s, tau_array[i], b, sigma, Nsteps)
      s+=val
      s2+=val**2
    mean=s/Ntrials
    var=s2/Ntrials-mean**2
    print i,'of', len(tau_array)
    print '\ttau=%1.2e  delta_t=%1.2e  delta_t/tau = %1.2e'%(tau_array[i], delta_t, delta_t/tau_array[i])
    print '\tmean %1.2e %1.2e'%(mean, exp(-delta_t/tau_array[i])*X_s+b*tau_array[i]*(1.-exp(-delta_t/tau_array[i] )))
    print '\tvar %1.2e %1.2e'%(var, (tau_array[i]*sigma**2)/2.*(1.-exp(-2.*delta_t/tau_array[i])))
    subplot(211)
    errorbar(tau_array[i],[mean],yerr=[sqrt(var)/sqrt(Ntrials)], fmt='b.-')
    subplot(212)
    errorbar(tau_array[i],[var], yerr=[var*sqrt(2/(Ntrials-1))], fmt='b.-')

  ax1=subplot(211)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  tau_array=10**arange(-5.,6.,0.1)
  plot(tau_array,exp(-delta_t/tau_array)*X_s+b*tau_array*(1.-exp(-delta_t/tau_array)),'b--', lw=2, label='delta_t=1 day')
  #plot(tau_array,b*tau_array*(1.-exp(-delta_t/tau_array)))
  #plot(tau_array,exp(-delta_t/tau_array)*X_s)
  xlabel('tau, days')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  ax1=subplot(212)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  plot(tau_array,(tau_array*sigma**2)/2.*(1-exp(-2.*delta_t/tau_array)),'b--', lw=2, label='delta_t=1 day')
  #plot(tau_array,(delta_t*sigma**2)*ones(len(tau_array)))
  #plot(tau_array,(tau_array*sigma**2)/2.)
  xlabel('tau, days')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  #ylim(1.e-5,1.e2)


  b=3.
  sigma=0.56
  X_s=2.
  delta_t=0.01
  tau_array=10**arange(-5.,6.,1.)
  #Nsteps=1000
  #Ntrials=1000
  figure(1)
  for i in range(0,len(tau_array)):
    s=0.
    s2=0.
    for j in range(0,Ntrials):
      val=X_int(delta_t,X_s, tau_array[i], b, sigma, Nsteps)
      s+=val
      s2+=val**2
    mean=s/Ntrials
    var=s2/Ntrials-mean**2
    print i,'of', len(tau_array)
    print '\ttau=%1.2e  delta_t=%1.2e  delta_t/tau = %1.2e'%(tau_array[i], delta_t, delta_t/tau_array[i])
    print '\tmean %1.2e %1.2e'%(mean, exp(-delta_t/tau_array[i])*X_s+b*tau_array[i]*(1.-exp(-delta_t/tau_array[i] )))
    print '\tvar %1.2e %1.2e'%(var, (tau_array[i]*sigma**2)/2.*(1.-exp(-2.*delta_t/tau_array[i])))
    subplot(211)
    errorbar(tau_array[i],[mean],yerr=[sqrt(var)/sqrt(Ntrials)], fmt='r.-')
    subplot(212)
    errorbar(tau_array[i],[var], yerr=[var*sqrt(2/(Ntrials-1))], fmt='r.-')

  ax1=subplot(211)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  tau_array=10**arange(-5.,6.,0.1)
  plot(tau_array,exp(-delta_t/tau_array)*X_s+b*tau_array*(1.-exp(-delta_t/tau_array)),'r--', lw=2, label='delta_t=.01 day')
  #plot(tau_array,b*tau_array*(1.-exp(-delta_t/tau_array)))
  #plot(tau_array,exp(-delta_t/tau_array)*X_s)
  xlabel('tau, days')
  ylabel(r'$E(X(t)|X(s))$')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  legend(loc=4)

  ax1=subplot(212)
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  plot(tau_array,(tau_array*sigma**2)/2.*(1-exp(-2.*delta_t/tau_array)),'r--', lw=2, label='delta_t=.01 days')
  #plot(tau_array,(delta_t*sigma**2)*ones(len(tau_array)))
  #plot(tau_array,(tau_array*sigma**2)/2.)
  xlabel('tau, days')
  ylabel(r'$Var(X(t)|X(s))$')
  xlim(3.e-6, 3.e5)
  ylim(1.e-6,1.e3)
  #ylim(1.e-5,1.e2)
  legend(loc=2)

  show()


def check_redshift_dependence():
  b=3.
  #delta_t=1.
  #tau=1.
  delta_t=1.
  tau=100.
  X_0=b*tau
  sigma=0.56
  redshift_array=arange(0.,5.,1.)
  Nsteps=1000
  t=arange(1000,dtype=float64)*delta_t
  for i in range(0,len(redshift_array)):
    X=drw_lightcurve(t, X_0, tau, sigma, b, redshift_array[i], Nsteps)
    dX=diff(X)
    print ''
    #print 'theoretical mean=%1.2e'%(b*tau*(1.-exp(-delta_t/tau)))
    print 'theoretical var=%1.2e'%(sigma**2 * tau/2. *(1.-exp(-2.*delta_t/tau)))
    print 'z=%1.2f\tvar=%1.2e'%(redshift_array[i],cumsum(dX**2)[len(dX)-1]/len(dX)-(cumsum(dX)[len(dX)-1]/len(dX))**2 )
    figure(1)
    n=510+i
    subplot(n)
    plot(t,X)
    xlabel('Observer Time, days')
    ylabel('Flux, arb. u.')
    
    suptitle(r'b=%1.2f,$\tau$=%1.2f, $\sigma$=%1.2f, $\Delta t$=%1.2f'%(b,tau,sigma, delta_t), fontsize=14)
    title(r'z=%1.1f'%(redshift_array[i]))
    figure(2)
    n=511+i
    subplot(n)
    plot(t[0:len(t)-1],dX)
    figure(3)
    n=510+i
    subplot(n)
    plot(X[0:len(t)-1],dX,'.', ms=3)
    #plot(X[0:len(t)-1],(exp(-delta_t/tau)-1.)*X[0:len(t)-1]+b*tau*(1.-exp(-delta_t/tau)),'r.', ms=3)
    X=np.sort(X)
    #errorbar(X[0:len(t)-1],(exp(-delta_t/tau)-1.)*X[0:len(t)-1]+b*tau*(1.-exp(-delta_t/tau)),yerr=sigma**2*tau/2*(1-exp(-2.*delta_t/tau)), fmt='r.-')
    #fill_between(x, y-error, y+error)
    fill_between(X[0:len(t)-1], (exp(-delta_t/tau)-1.)*X[0:len(t)-1]+b*tau*(1.-exp(-delta_t/tau)) - sigma**2*tau/2*(1-exp(-2.*delta_t/tau)), (exp(-delta_t/tau)-1.)*X[0:len(t)-1]+b*tau*(1.-exp(-delta_t/tau)) + sigma**2*tau/2*(1-exp(-2.*delta_t/tau)), alpha=0.7, edgecolor='r', facecolor='r')
    #show()
  show()  
  

#SELECT WHICH DEPENDENCE YOU WANT TO TEST
check_sigma_dependence()
#check_tau_dependence()
#check_redshift_dependence()
