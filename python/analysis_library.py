from fot_library import * 
from pylab import *
import numpy
import glob
import triangle
import matplotlib.pyplot as plt 
from astropy.io import ascii
from scipy.optimize import curve_fit

def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))

def get_data(fnm, season, img1, img2, pmu):
	fname = '/data2/fot/data/cosmograil/' +fnm+'_season%d'%(season) + '.rdb'
	data=ascii.read(fname,data_start=2)
    	return data['mhjd'], data['mag_'+img1], data['magerr_'+img1], data['mag_'+img2], data['magerr_'+img2]

def get_chain(filename):
	npzfile = numpy.load(filename)
	print 'npzfile.files', npzfile.files
	samples = npzfile['arr_0']
	return samples

def get_parameter_samples(samples):
	return samples[:,0], samples[:,1], samples[:,2], samples[:,3], samples[:,4]

def filter_samples(array, n_iterations, n_iteration_filter):
	array_index = range(0,len(array))
	return array[mod(array_index,n_iterations)>n_iteration_filter]
	
def hist_param(array, y, norm=True):
	p,d =  numpy.histogram(array, bins = len(y), range = [min(y),max(y)] )
	delta_d = d[1]-d[0]
	d = d[1:len(d)] - delta_d
	p = p.astype(float)
	if(norm): p/=max(p) 
	return p,d

def test(chain_file_name):
	samples = get_chain(chain_file_name)
	delay, d_mag, sigma, tau, avg_mag = get_parameter_samples(dirc, fnm, img1, img2, pmu)
	nwalkers     = 100
	n_iterations = 1000
	k = range(0,len(delay))
	ax = subplot(221)
	hist(delay, bins = 2000, range = [-200, 200.] )
	subplot(223)
	plot(delay)
	delay = delay[mod(k,1000)>200]
	subplot(222, sharex = ax)
	hist(delay, bins = 2000, range = [-200, 200.] )
	subplot(224)
	plot(delay)



