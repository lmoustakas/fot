#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
Full Of Time
'''
if __name__ == "__main__":
    
    from fot_library import * 
    import argparse
    parser=argparse.ArgumentParser(description='fot_hubble_sim routine to simulat hubble observations and calculate delay inference')
    
    parser.add_argument("-dt",    "--delay",                help="Delay between images",type=float)
    parser.add_argument("-dtp",   "--delay_prior",          help="Delay prior for emcee",type=float)
    parser.add_argument("-dtpmin","--delay_prior_min",      help="Delay prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dtpmax","--delay_prior_max",      help="Delay prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-dm",    "--delta_mag",            help="Difference in magntiude between images",type=float)
    parser.add_argument("-dmp",   "--delta_mag_prior",      help="Difference in magniitude between images prior for emcee",type=float)
    parser.add_argument("-dmpmin","--delta_mag_prior_min",  help="Difference in magniitude between images prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dmpmax","--delta_mag_prior_max",  help="Difference in magniitude between images prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-s",      "--sigma",               help="Quasar light curve sigma parameter",type=float)
    parser.add_argument("-sp",     "--sigma_prior",         help="Quasar light curve sigma prior for emcee",type=float)
    parser.add_argument("-spmin",  "--sigma_prior_min",     help="Quasar light curve sigma prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-spmax",  "--sigma_prior_max",     help="Quasar light curve sigma prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-t",      "--tau",                 help="Quasar light curve tau parameter",type=float)
    parser.add_argument("-tp",     "--tau_prior",           help="Quasar light curve prior for emcee",type=float)
    parser.add_argument("-tpmin",  "--tau_prior_min",       help="Quasar light curve prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-tpmax",  "--tau_prior_max",       help="Quasar light curve prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-m",     "--avg_mag",              help="Quasar light curve average magnitude",type=float)
    parser.add_argument("-mp",    "--avg_mag_prior",        help="Quasar light curve magnitude prior for emcee",type=float)
    parser.add_argument("-mpmin", "--avg_mag_prior_min",    help="Quasar light curve magnitude prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-mpmax", "--avg_mag_prior_max",    help="Quasar light curve magnitude prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-z",      "--redshift",            help="Quasar redshift",type=float)
    parser.add_argument("-n",      "--numorbits",           help="Number of Hubble orbits to be used",type=int)
    parser.add_argument("-p",      "--photometric_error",   help="Photometric error uncertainty (e.g. 0.02) [mag]",type=float)
    parser.add_argument("-o",      "--output_tag",          help="Output tag, in quotes",type=str)
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    #PAESE ARGUMENTS
    args=parser.parse_args()
    #Nsteps is a simulation parameter. It is recommended to leave it at 1000.
    Nsteps=1000
    print ''
    print 'SIMULATION PARAMETERS' 
    print args
    print ''
    time_array = hubble_obs_sched(args.numorbits)
    mag1, mag2 = delayed_lightcurve(time_array, args.delay, args.delta_mag, args.redshift, args.tau, args.avg_mag, args.sigma, Nsteps)
    #mag1, mag2 = delayed_lightcurve(time_array, delay, delta_mag, redshift, tau, avg_mag, sigma, Nsteps)
    mag1_dat = mag1 + array([random.gauss(0.,args.photometric_error) for _ in xrange(len(time_array))])
    mag2_dat = mag2 + array([random.gauss(0.,args.photometric_error) for _ in xrange(len(time_array))])
    err1 = args.photometric_error*ones(len(time_array))
    err2 = args.photometric_error*ones(len(time_array))
    
    import numpy
    import datetime
    t0 = datetime.datetime.now()
    date_string = t0.strftime("_%Y_%m_%d_%H:%M")
    outputdir=os.environ['FOTDIR']+'/outputs/'
    numpy.savez(outputdir+'sim_data_%s%s'%(args.output_tag,date_string), time_array=time_array, mag1=mag1, mag2=mag2, mag1_dat=mag1_dat, mag2_dat=mag2_dat)
    
    emcee_delay_estimator(time_array, mag1_dat,err1,mag2_dat,err2,args.output_tag, 
			  args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
			  args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
    #show()

    
    # define the right column header names
    #mag1, mag2, magerr1, magerr2 = 'mag_'+args.image1, 'mag_'+args.image2, 'magerr_'+args.image1, 'magerr_'+args.image2
    # read in the data
    #time,m,me=read_cosmograil_data(args.datafile,[mag1,mag2],[magerr1,magerr2])
    # add optionally specified systematic uncertainty
    #me[magerr1] = np.add(me[magerr1],args.systematic)
    #me[magerr2] = np.add(me[magerr2],args.systematic)
    # run the emcee delay estimator
    #N=len(time)
    #print 'N=',N
    #n1=N-70
    #n2=N
    #emcee_delay_estimator(time[n1:n2], m[mag1][n1:n2],me[magerr1][n1:n2],m[mag2][n1:n2],me[magerr2][n1:n2],args.outputtag)
