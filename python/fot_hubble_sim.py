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
    
    parser.add_argument("-dt","--delay",help="Delay between images",type=float)
    parser.add_argument("-dm","--delta_mag",help="Difference in magntiude between images",type=float)
    parser.add_argument("-s","--sigma",help="Quasar light curve sigma parameter",type=float)
    parser.add_argument("-t","--tau",help="Quasar light curve tau parameter",type=float)
    parser.add_argument("-m","--avg_mag",help="Quasar light curve average magnitude",type=float)
    parser.add_argument("-z","--redshift",help="Quasar redshift",type=float)
    parser.add_argument("-n","--numorbits",help="Number of Hubble orbits to be used",type=int)
    parser.add_argument("-p","--photometricerror",help="Photometric error uncertainty (e.g. 0.02) [mag]",type=float)
    parser.add_argument("-o","--outputtag",help="Output tag, in quotes",type=str)
    parser.add_argument("-dtp","--delayprior",help="Delay prior for emcee [optional]",type=float)
    args=parser.parse_args()
    Nsteps=1000
    print ''
    print 'SIMULATION PARAMETERS' 
    print args
    print ''
    time_array = hubble_obs_sched(args.numorbits)
    mag1, mag2 = delayed_lightcurve(time_array, args.delay, args.delta_mag, args.redshift, args.tau, args.avg_mag, args.sigma, Nsteps)
    #mag1, mag2 = delayed_lightcurve(time_array, delay, delta_mag, redshift, tau, avg_mag, sigma, Nsteps)
    mag1_dat = mag1 + array([random.gauss(0.,args.photometricerror) for _ in xrange(len(time_array))])
    mag2_dat = mag2 + array([random.gauss(0.,args.photometricerror) for _ in xrange(len(time_array))])
    err1 = args.photometricerror*ones(len(time_array))
    err2 = args.photometricerror*ones(len(time_array))
    
    import numpy
    import datetime
    t0 = datetime.datetime.now()
    date_string = t0.strftime("_%Y_%m_%d_%H:%M")
    print numpy.arange(0.,1.,0.1)
    numpy.savez('sim_data_%s%s'%(args.outputtag,date_string), time_array=time_array, mag1=mag1, mag2=mag2, mag1_dat=mag1_dat, mag2_dat=mag2_dat)
    #numpy.savez('test', time_array, mag1, mag2, mag1_dat, mag2_dat, 'time_array', 'mag1', 'mag2', 'mag1_dat', 'mag2_dat')
    #errorbar(time_array, lc1, err1, fmt='b.')
    #errorbar(time_array, lc2, err2, fmt='r.')
    emcee_delay_estimator(time_array, mag1_dat,err1,mag2_dat,err2,args.outputtag, delay_prior = args.delayprior)
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
