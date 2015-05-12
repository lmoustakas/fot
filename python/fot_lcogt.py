#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
Full Of Time
'''

import matplotlib
matplotlib.use('Agg') # Note, this MUST be before importing pylab or matplotlib.pyplot


if __name__ == "__main__":
    
    from fot_library import * 
    import argparse
    parser=argparse.ArgumentParser(description='fot_delay routine to calculate delay inference')
    parser.add_argument("-i","--datafile",	help="time delay challenge data file",type=str)
    parser.add_argument("-l","--image1",	default='A', help="Image 1 name (e.g. 'A')",type=str)
    parser.add_argument("-m","--image2",	default='B', help="Image 2 name (e.g. 'B')",type=str)
    parser.add_argument("-su","--systematic",	default=0., help="Additional systematic uncertainty (e.g. 0.1) [mag]",type=float)

    
    parser.add_argument("-dtp",   "--delay_prior",          default = 0., help="Delay prior for emcee",type=float)
    parser.add_argument("-dtpmin","--delay_prior_min",      default = -1.5, help="Delay prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dtpmax","--delay_prior_max",      default = +1.5, help="Delay prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-dmp",   "--delta_mag_prior",      default = 0., help="Difference in magniitude between images prior for emcee",type=float)
    parser.add_argument("-dmpmin","--delta_mag_prior_min",  default = -10., help="Difference in magniitude between images prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dmpmax","--delta_mag_prior_max",  default = +10., help="Difference in magniitude between images prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-sp",     "--sigma_prior",         default = 0.1, help="Quasar light curve sigma prior for emcee",type=float)
    parser.add_argument("-spmin",  "--sigma_prior_min",     default = 7.e-5, help="Quasar light curve sigma prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-spmax",  "--sigma_prior_max",     default = 7.e1, help="Quasar light curve sigma prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-tp",     "--tau_prior",           default = 100., help="Quasar light curve prior for emcee",type=float)
    parser.add_argument("-tpmin",  "--tau_prior_min",       default = 1., help="Quasar light curve prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-tpmax",  "--tau_prior_max",       default = 1.e5, help="Quasar light curve prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-mp",    "--avg_mag_prior",        default = 18., help="Quasar light curve magnitude prior for emcee",type=float)
    parser.add_argument("-mpmin", "--avg_mag_prior_min",    default = 0., help="Quasar light curve magnitude prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-mpmax", "--avg_mag_prior_max",    default = 100., help="Quasar light curve magnitude prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-z",      "--redshift",            default = 0., help="Quasar redshift",type=float)
    parser.add_argument("-o",      "--output_tag",          default = 'tmp', help="Output tag, in quotes",type=str)
    
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    # define the right column header names
    print args
    
    # read in the data
    #print rung, args.datafile
    import os
    fname = '/home/romerowo/lcolens/data/he0435-1223_lcogt_magnitudes.dat'

    if(os.path.isfile(fname)==False):
    	print fname
    	print 'File not found.'
    	exit()

    data=ascii.read(fname,data_start=6)
    print data
    #exit()
    t=data['mjd']
    m1=data['mag_A']
    em1=data['magerr_A']
    m2=data['mag_B']
    em2=data['magerr_B']
    m3=data['mag_C']
    em3=data['magerr_C']
    m4=data['mag_D']
    em4=data['magerr_D']
    
    '''
    t = array(t)
    m1=array(m1)
    m2=array(m2)
    em1=array(em1)
    em2=array(em2)
    em1 = np.add(em1,args.systematic)
    em2 = np.add(em2,args.systematic)
    
    t,m1,m2,em1,em2 = remove_outliers(t,m1,m2,em1,em2, display=False)
    '''
    sig1, log10_tau1 = emcee_lightcurve_estimator(t, m1, em1, args.output_tag+'_lc1', 
			  args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)

    sig2, log10_tau2 = emcee_lightcurve_estimator(t, m2, em2, args.output_tag+'_lc2', 
              args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
              args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)

    sig3, log10_tau3 = emcee_lightcurve_estimator(t, m3, em3, args.output_tag+'_lc3', 
              args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
              args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    
    sig4, log10_tau4 = emcee_lightcurve_estimator(t, m4, em4, args.output_tag+'_lc4', 
              args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
              args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
   
    '''
    '''
    tau1=10.**log10_tau1
    tau2=10.**log10_tau2
    tau3=10.**log10_tau3
    tau4=10.**log10_tau4
    sig=0.25*(sig1+sig2+sig3+sig4)
    tau = 10.**(0.25*(log10_tau1+log10_tau2+log10_tau3+log10_tau4)) #geometric mean
    print '#############################'
    print 'sigma_1\t\tsigma2\t\ttau1\t\ttau2'
    print '%1.2e\t%1.2e\t%1.2e\t%1.2e'%(sig1, sig2, tau1, tau2)
    print '#############################'
    emcee_delay_estimator(t, m1, em1, m2, em2 ,args.output_tag+'_dt12', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)

    emcee_delay_estimator(t, m1, em1, m3, em3 ,args.output_tag+'_dt13', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    emcee_delay_estimator(t, m1, em1, m4, em4 ,args.output_tag+'_dt14', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    emcee_delay_estimator(t, m2, em2, m3, em3 ,args.output_tag+'_dt23', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    emcee_delay_estimator(t, m2, em2, m4, em4 ,args.output_tag+'_dt24', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    emcee_delay_estimator(t, m3, em3, m4, em4 ,args.output_tag+'_dt34', 
              args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
              args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
              sig,            args.sigma_prior_min,    args.sigma_prior_max, 
              tau,                  min(tau1,tau2),      args.tau_prior_max, 
              args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    
    '''
    emcee_delay_estimator(t, m1, em1, m2, em2 ,args.output_tag+'_dt', 
			  args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
			  args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  args.sigma_prior,	args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,	args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max,
			  poly=3)
    '''
