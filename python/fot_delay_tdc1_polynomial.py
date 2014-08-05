#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
Full Of Time
'''



if __name__ == "__main__":
    
    from fot_library import * 
    import argparse
    parser=argparse.ArgumentParser(description='fot_delay routine to calculate delay inference')
    parser.add_argument("-i","--datafile",	help="time delay challenge data file",type=str)
    parser.add_argument("-l","--image1",	default='A', help="Image 1 name (e.g. 'A')",type=str)
    parser.add_argument("-m","--image2",	default='B', help="Image 2 name (e.g. 'B')",type=str)
    parser.add_argument("-su","--systematic",	default=0., help="Additional systematic uncertainty (e.g. 0.1) [mag]",type=float)

    
    parser.add_argument("-dtp",   "--delay_prior",          default = 0., help="Delay prior for emcee",type=float)
    parser.add_argument("-dtpmin","--delay_prior_min",      default = -1600., help="Delay prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dtpmax","--delay_prior_max",      default = +1600., help="Delay prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-dmp",   "--delta_mag_prior",      default = 0., help="Difference in magniitude between images prior for emcee",type=float)
    parser.add_argument("-dmpmin","--delta_mag_prior_min",  default = -10., help="Difference in magniitude between images prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dmpmax","--delta_mag_prior_max",  default = +10., help="Difference in magniitude between images prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-sp",     "--sigma_prior",         default = 1., help="Quasar light curve sigma prior for emcee",type=float)
    parser.add_argument("-spmin",  "--sigma_prior_min",     default = 7.e-5, help="Quasar light curve sigma prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-spmax",  "--sigma_prior_max",     default = 7.e1, help="Quasar light curve sigma prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-tp",     "--tau_prior",           default = 3500., help="Quasar light curve prior for emcee",type=float)
    parser.add_argument("-tpmin",  "--tau_prior_min",       default = 10., help="Quasar light curve prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-tpmax",  "--tau_prior_max",       default = 1.e5, help="Quasar light curve prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-mp",    "--avg_mag_prior",        default = -13., help="Quasar light curve magnitude prior for emcee",type=float)
    parser.add_argument("-mpmin", "--avg_mag_prior_min",    default = -30., help="Quasar light curve magnitude prior minimum acceptable value for emcee",type=float)
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
    if('rung' not in args.datafile or 'pair' not in args.datafile):
	print '%s is not a valid filename for fot_delay_tdc1.py'
	exit()
    rung = int(args.datafile.split('rung')[1].split('_')[0])
    #pair = int(args.datafile.split('pair')[1].split('.')[0])
    print rung, args.datafile
    import os
    fname = os.environ['FOTDIR']+'/data/tdc1/' + 'rung%d/'%rung + args.datafile

    if(os.path.isfile(fname)==False):
	print fname
	print 'File not found.'
	exit()

    data=ascii.read(fname,data_start=6)
    t_in=data['col1']
    f1=data['col2']
    ef1=data['col3']
    f2=data['col4']
    ef2=data['col5']
       
    t=[]
    m1=[]
    m2=[]
    em1=[]
    em2=[]
    for k in range(0,len(t_in)):
	if(f1[k]>0. and f2[k]>0.):
          t.append(t_in[k])
	  m1.append(22.5-2.5*log10(f1[k]))
	  em1.append(2.5*ef1[k]/f1[k]/log(10.))
	  m2.append(22.5-2.5*log10(f2[k]))
	  em2.append(2.5*ef2[k]/f2[k]/log(10.))
    t = array(t)
    m1=array(m1)
    m2=array(m2)
    em1=array(em1)
    em2=array(em2)
    em1 = np.add(em1,args.systematic)
    em2 = np.add(em2,args.systematic)
    
    t,m1,m2,em1,em2 = remove_outliers(t,m1,m2,em1,em2)
    '''
    t = t[::20]
    m1 = m1[::20]
    m2 = m2[::20]
    em1 = em1[::20]
    em2 = em2[::20]
    '''
    #N = len(t)
    #N_sub = int(0.05*N)
    #print 'number of points', N
    #print 'subsample number of points', N_sub
    '''
    from numpy.random import randint
    re_sample = 16
    offset = randint(0,re_sample)
    t   = (t[offset:])[::re_sample]
    m1  = (m1[offset:])[::re_sample]
    m2  = (m2[offset:])[::re_sample]
    em1 = (em1[offset:])[::re_sample]
    em2 = (em2[offset:])[::re_sample]
    '''
    #print t
    #exit()
    '''
    indices = randint(0,N,N_sub)
    t   = t[indices]
    m1  = m1[indices]
    m2  = m2[indices]
    em1 = em1[indices]
    em2 = em2[indices]
    '''
    # run the emcee delay estimator
    #emcee_delay_estimator(time, m[mag1],me[magerr1],m[mag2],me[magerr2],args.outputtag)
    '''
    sig1, log10_tau1 = emcee_lightcurve_estimator(t, m1, em1, args.output_tag+'_lc1', 
			  args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)

    sig2, log10_tau2 = emcee_lightcurve_estimator(t, m2, em2, args.output_tag+'_lc2', 
			  args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    tau1=10.**log10_tau1
    tau2=10.**log10_tau2
    sig=0.5*(sig1+sig2)
    tau = 10.**(0.5*(log10_tau1+log10_tau2)) #geometric mean
    print '#############################'
    print 'sigma_1\t\tsigma2\t\ttau1\t\ttau2'
    print '%1.2e\t%1.2e\t%1.2e\t%1.2e'%(sig1, sig2, tau1, tau2)
    print '#############################'
    emcee_delay_estimator(t, m1, em1, m2, em2 ,args.output_tag+'_dt', 
			  args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
			  args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  sig,    		args.sigma_prior_min,    args.sigma_prior_max, 
			  tau,                  min(tau1,tau2),      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    '''
    emcee_delay_estimator(t, m1, em1, m2, em2 ,args.output_tag+'_dt', 
			  args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
			  args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  args.sigma_prior,	args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,	args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max,
			  poly=5)
   
