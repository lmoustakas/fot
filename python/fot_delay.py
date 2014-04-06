#!/usr/bin/env python
'''
2014 April 5
ARW & LAM, JPL/Caltech
Full Of Time
'''
if __name__ == "__main__":
    
    import fot_library
    import argparse
    parser=argparse.ArgumentParser(description='fot_delay routine to calculate delay inference')
    parser.add_argument("-i","--datafile",help="COSMOGRAIL data file",type=str)
    parser.add_argument("-l","--image1",help="Image 1 name (e.g. 'A')",type=str)
    parser.add_argument("-m","--image2",help="Image 2 name (e.g. 'B')",type=str)
    parser.add_argument("-o","--outputtag",help="Output tag, in quotes",type=str)
    parser.add_argument("-s","--systematic",help="Additional systematic uncertainty (e.g. 0.1) [mag]",type=float)
    args=parser.parse_args()
    infile=args.infile
    # read in the data
    mag1, mag2, magerr1, magerr2 = 'mag_'+args.image1, 'mag_'+args.image2, 'magerr_'+args.image1, 'magerr_'+args.image2
    time,m,me=read_cosmograil_data(args.datafile,[mag1,mag2],[magerr1,magerr2])
    # add optionally specified systematic uncertainty
    me=me+args.systematic
    emcee_delay_estimator(time, m['mag1'],me['magerr1'],m['mag2'],me['magerr2'],args.outputtag)
