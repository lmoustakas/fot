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
    parser.add_argument("-i","--datafile",help="COSMOGRAIL data file",type=str)
    parser.add_argument("-l","--image1",help="Image 1 name (e.g. 'A')",type=str)
    parser.add_argument("-m","--image2",help="Image 2 name (e.g. 'B')",type=str)
    parser.add_argument("-su","--systematic",help="Additional systematic uncertainty (e.g. 0.1) [mag]",type=float)

    
    parser.add_argument("-dtp",   "--delay_prior",          help="Delay prior for emcee",type=float)
    parser.add_argument("-dtpmin","--delay_prior_min",      help="Delay prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dtpmax","--delay_prior_max",      help="Delay prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-dmp",   "--delta_mag_prior",      help="Difference in magniitude between images prior for emcee",type=float)
    parser.add_argument("-dmpmin","--delta_mag_prior_min",  help="Difference in magniitude between images prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-dmpmax","--delta_mag_prior_max",  help="Difference in magniitude between images prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-sp",     "--sigma_prior",          help="Quasar light curve sigma prior for emcee",type=float)
    parser.add_argument("-spmin",  "--sigma_prior_min",      help="Quasar light curve sigma prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-spmax",  "--sigma_prior_max",      help="Quasar light curve sigma prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-tp",     "--tau_prior",            help="Quasar light curve prior for emcee",type=float)
    parser.add_argument("-tpmin",  "--tau_prior_min",        help="Quasar light curve prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-tpmax",  "--tau_prior_max",        help="Quasar light curve prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-mp",    "--avg_mag_prior",        help="Quasar light curve magnitude prior for emcee",type=float)
    parser.add_argument("-mpmin", "--avg_mag_prior_min",    help="Quasar light curve magnitude prior minimum acceptable value for emcee",type=float)
    parser.add_argument("-mpmax", "--avg_mag_prior_max",    help="Quasar light curve magnitude prior maximum acceptable value for emcee",type=float)
    parser.add_argument("-z",      "--redshift",             help="Quasar redshift",type=float)
    parser.add_argument("-o",      "--output_tag",           help="Output tag, in quotes",type=str)
    
    #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
    for i, arg in enumerate(sys.argv):
      if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    args=parser.parse_args()
    # define the right column header names
    print args
    
    mag1, mag2, magerr1, magerr2 = 'mag_'+args.image1, 'mag_'+args.image2, 'magerr_'+args.image1, 'magerr_'+args.image2
    # read in the data
    time,m,me=read_cosmograil_data(args.datafile,[mag1,mag2],[magerr1,magerr2])
    # add optionally specified systematic uncertainty
    me[magerr1] = np.add(me[magerr1],args.systematic)
    me[magerr2] = np.add(me[magerr2],args.systematic)
    # run the emcee delay estimator
    #emcee_delay_estimator(time, m[mag1],me[magerr1],m[mag2],me[magerr2],args.outputtag)

    emcee_delay_estimator(time, m[mag1],me[magerr1],m[mag2],me[magerr2],args.output_tag, 
			  args.delay_prior,    args.delay_prior_min,    args.delay_prior_max,
			  args.delta_mag_prior, args.delta_mag_prior_min, args.delta_mag_prior_max, 
			  args.sigma_prior,    args.sigma_prior_min,    args.sigma_prior_max, 
			  args.tau_prior,      args.tau_prior_min,      args.tau_prior_max, 
			  args.avg_mag_prior,   args.avg_mag_prior_min,   args.avg_mag_prior_max)
   