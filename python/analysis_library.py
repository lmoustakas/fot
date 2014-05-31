from fot_library import * 
from pylab import *
import numpy
import glob
import triangle
import matplotlib.pyplot as plt 
from astropy.io import ascii
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import beta


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

def read_sim_hubble_light_curve(filename):
	npzfile = numpy.load(filename)
	#print 'npzfile.files', npzfile.files
	mag1 = npzfile['mag1']
	mag2 = npzfile['mag2']
	time_array = npzfile['time_array']
	mag1_dat = npzfile['mag1_dat']
	mag2_dat = npzfile['mag2_dat']
	return mag1, mag2, time_array, mag1_dat, mag2_dat

def rms(array):
  s2=cumsum(array**2)[len(array)-1]/len(array)
  r=s2-mean(array)**2
  return sqrt(r)

def fit_val(array, numbins, minval, maxval):
	hist, bin_edges = numpy.histogram(array, bins=numbins, range=[minval,maxval], density=False)
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
	p0 = [20., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
	coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
	# Get the fitted curve
	hist_fit = gauss(bin_centres, *coeff)
	plot(bin_centres, hist, lw=3, drawstyle='steps-mid', label='Test data')
	#plot(bin_centres, cumsum(hist), lw=3, drawstyle='steps-mid', label='Test data')
	x=arange(minval,maxval,(maxval-minval)*1.e-3)
	plot(x, gauss(x,*coeff), 'r-', lw=2)
	xticks(rotation=45 )
	xlabel('rec.-true delays, hours')
	ylabel('Count (total of 100)')

	#print coeff
	print '\tFitted mean = %1.2e'%(coeff[1])
	print '\tFitted standard deviation = %1.2e'%(coeff[2])
	xlim(minval,maxval)
	return coeff[1],coeff[2]

def hist_ana(dir1, run1, run2, nbins, norbits, photom_err,ax):
	print 'N=%d orbits, pme = %1.2f mag'%(norbits, photom_err)
	run, del_val, err_up, err_low, ll = read_logs(dir1,run1, run2)
	m,s=fit_val((del_val-1.5)*24., nbins, -30.,30.)
	scs = len(del_val[abs(del_val-1.5)*24.<3.*s])
	print '\tSuccess rate %1.0f%s'%(scs,'%')
	text(0.97, 0.85,r'%d%s success'%(scs,'%'), ha='right', va='center', transform=ax.transAxes,fontsize=13)
	text(0.97, 0.65,r'$\mu$=%1.2f hours'%(m), ha='right', va='center', transform=ax.transAxes,fontsize=13)
	text(0.97, 0.45,r'$\sigma$=%1.2f hours'%(s), ha='right', va='center', transform=ax.transAxes,fontsize=13)
	title('%d orbits, %1.2f mag uncertainty'%(norbits, photom_err))
	xticks([-24, -12, 0, 12, 24])
	return m,s,scs

def binom_interval(success, total, confint=0.68):
	quantile = (1 - confint) / 2.
	lower = beta.ppf(quantile, success, total - success + 1)
	upper = beta.ppf(1 - quantile, success + 1, total - success)
	if(success==total):
		upper=1.0
	if(success==0):
		lower=0.
	return (lower*total, upper*total)

def read_logs(dirc,num_orbits, ph_uncrtnty):
	os.chdir(dirc)
	os.system('pwd')
	fnames = os.listdir('./')
	#print 'in read_logs: dirc, num_orbits, ph_uncrtnty', dirc, num_orbits, ph_uncrtnty
	#print fnames
	run=[]
	del_val=[]
	err_up=[]
	err_low=[]
	ll=[]
	for f in fnames:
	  if(len(f.split('.'))>0 and len(f.split('_'))>4):
  	    tag=f.split('.')[len(f.split('.'))-1]
	    norb = (f.split('_')[2])
	    pu = (f.split('_')[4])
	    if (tag=='log' and norb == num_orbits and ph_uncrtnty == pu and os.stat(f).st_size!=0):
	      print f
	      if(os.stat(f).st_size!=0):
	        #print int(f.split('.')[0].split('_r')[1])
	        run_num = int(f.split('.')[0].split('_r')[1])
	        #print r1,r2,f
	        lc=0
	        mc_table=0
	        for line in file(f):
		  lc+=1
		  if('MC_rec' in line): mc_table+=1
		  if('delay' in line and mc_table==2):
		    #print lc, line
		    run.append(run_num)
		    del_val.append(float(line.split()[1]))
		    err_up.append(float(line.split()[2]))
		    err_low.append(abs(float(line.split()[3])))
		  if('Likelihood' in line):
		    ll.append(float(line.split(':')[1]))
	del_val=array(del_val)
	err_up=array(err_up)
	err_low=array(err_low)
	ll=array(ll)

	return run, del_val, err_up, err_low, ll



def ana_ll2(dirc, ll_cut, sig_cut, numbins, minval, maxval, norbits, photom_err, ax):
	print 'N=%d orbits, pme = %1.2f mag'%(norbits, photom_err)
	#run, del_val, err_up, err_low, ll = read_logs(dirc,run1, run2)
	print 'in ana_ll2: norbits, photom_err', norbits, photom_err
	pmu = '%1.2f'%(photom_err)
	pmu = pmu.replace('.','p')
	num_orb = '%s'%(norbits)
	run, del_val, err_up, err_low, ll = read_logs(dirc, num_orb, pmu)
	print run, del_val, err_up, err_low, ll
	array=(del_val-1.5)*24.
	hist, bin_edges = numpy.histogram(array, bins=numbins, range=[minval,maxval], density=False)
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

	plot(bin_centres, hist, lw=3, drawstyle='steps-mid', label='Simulation', color='0.6')
	p0 = [20., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
	coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
	# Get the fitted curve
	hist_fit = gauss(bin_centres, *coeff)
	x=arange(minval,maxval,(maxval-minval)*1.e-3)
	#plot(x, gauss(x,*coeff), 'r-', lw=2)
	print '\tmu    = %1.2e'%(coeff[1])
	print '\tsigma = %1.2e'%(coeff[2])
	title('Delay Reconsruction')
	#text(0.97, 0.8,r'$\mu$=%1.2f hours'%(coeff[1]), ha='right', va='center', transform=ax.transAxes,fontsize=13)
	#text(0.97, 0.7,r'$\sigma$=%1.2f hours'%(coeff[2]), ha='right', va='center', transform=ax.transAxes,fontsize=13)


	#subplot(223)
	#errorbar(ll,array,yerr=(err_low*24.,err_up*24.),fmt='.', ms=3)
	#ylabel('reconstructed - true delay, hours')
	#xlabel('log likelihood')
	#plot([ll_cut,ll_cut],[axis()[2],axis()[3]], 'r--')
	#xticks(rotation=45 )


	if(ll_cut>0.):
		#CUTS
		print '\tscs   = %d'%(len(ll[ ll>ll_cut]))
		cut2_list = [k for k in range(0,len(ll)) if (ll[k] > ll_cut and (abs(del_val[k]-1.5)*24.>sig_cut*coeff[2] or err_low[k]*24.>sig_cut*coeff[2] or err_up[k]*24.>sig_cut*coeff[2]))]
		print '\tscs   = %d'%(len(ll)-len(cut2_list))
		cut_tot_list = [k for k in range(0,len(ll)) if (ll[k] < ll_cut or (abs(del_val[k]-1.5)*24.>sig_cut*coeff[2] or err_low[k]*24.>sig_cut*coeff[2] or err_up[k]*24.>sig_cut*coeff[2]))]
		accept_tot_list = [k for k in range(0,len(ll)) if (ll[k] > ll_cut and (abs(del_val[k]-1.5)*24.<sig_cut*coeff[2] and err_low[k]*24.<sig_cut*coeff[2] and err_up[k]*24.<sig_cut*coeff[2]))]
		print '\tscs   = %d'%(len(ll)-len(cut_tot_list))

		#print axis()
		#ax=subplot(222)
		hist, bin_edges = numpy.histogram(array[accept_tot_list], bins=numbins, range=[minval,maxval], density=False)
		plot(bin_centres, hist, lw=3, drawstyle='steps-mid', label='Successful')
		p0 = [20., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
		coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
		# Get the fitted curve
		hist_fit = gauss(bin_centres, *coeff)
		x=arange(minval,maxval,(maxval-minval)*1.e-3)
		plot(x, gauss(x,*coeff), 'r-', lw=2, label='Fit')
		legend(loc=2)
		print '\tmu    = %1.2e'%(coeff[1])
		print '\tsigma = %1.2e'%(coeff[2])
		title('Delay Reconsruction With Cuts')
		text(0.97, 0.90,r'%d%s success'%(len(ll[accept_tot_list]),'%'), ha='right', va='center', transform=ax.transAxes,fontsize=13)
		text(0.97, 0.75,r'$\mu$=%1.2f hours'%(coeff[1]), ha='right', va='center', transform=ax.transAxes,fontsize=13)
		text(0.97, 0.60,r'$\sigma$=%1.2f hours'%(coeff[2]), ha='right', va='center', transform=ax.transAxes,fontsize=13)

		#subplot(224)
		#errorbar(ll[ll>ll_cut],array[ll>ll_cut],yerr=(err_low[ll>ll_cut]*24.,err_up[ll>ll_cut]*24.),fmt='.', ms=3)
		#if(len(cut2_list)!=0):
			#errorbar(ll[cut2_list], array[cut2_list], yerr=(err_low[cut2_list]*24.,err_up[cut2_list]*24.),  fmt='r.', ecolor='r', ms=3)
		#xticks(rotation=45 )
		#xlim(ll_cut, max(ll))
		#ylabel('rec. - true delay, hours')
		#xlabel('log likelihood')

	xticks([-6, -3, 0, 3, 6], rotation=45 )
	xlim(-6.,6.)

	#xticks([-12, -9, -6, -3, 0, 3, 6, 9, 12], rotation=45 )
	xlabel('reconstructed - true delay, hours')
	title('%d Orbits\n%1.2f Photometric Uncertainty'%(norbits,photom_err), fontsize=16)
	#title('%1.2f Photometric Uncertainty\n%d orbits'%(photom_err,norbits), fontsize=16)

	return coeff[1], coeff[2], len(ll[accept_tot_list])

def ll_plot2():
	pdfp=PdfPages('systematics_ll2.pdf')
	#dir1 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_08_hubble_sim_runs/'
	#dir2 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_09_hubble_sim_runs/' 
	#dir3 = '/data2/fot/python/'
	dir1 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_15_hubble_sim_runs/'
	sig_cut = 3.0

	# 40 ORBITS
	fig=figure()
	ax=subplot(325)
	m_40_0p02, s_40_0p02, scs_40_0p02 = ana_ll2(dir1, 800., sig_cut, 250, -30., 30., 40, 0.02, ax)
	ax=subplot(326)
	m_40_0p05, s_40_0p05, scs_40_0p05 = ana_ll2(dir1,500., sig_cut, 100, -30., 30., 40, 0.05, ax)
	
	# 60 ORBITS
	ax=subplot(323)
	m_60_0p02, s_60_0p02, scs_60_0p02 = ana_ll2(dir1, 1200., sig_cut, 250, -30., 30., 60, 0.02,ax)
	ax=subplot(324)
	m_60_0p05, s_60_0p05, scs_60_0p05 = ana_ll2(dir1,  850., sig_cut, 100, -30., 30., 60, 0.05,ax)

	# 80 ORBITS
	ax=subplot(321)
	m_80_0p02, s_80_0p02, scs_80_0p02 = ana_ll2(dir1, 1700., sig_cut, 250, -30., 30., 80, 0.02,ax)
	ax=subplot(322)
	m_80_0p05, s_80_0p05, scs_80_0p05 = ana_ll2(dir1, 1100., sig_cut, 100, -30., 30., 80, 0.05,ax)


	subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.1, hspace=0.99)
	show()
	pdfp.savefig(fig)
#ll_plot2()
#exit()

def ana_ll(dirc, num_orbits, ph_unc, ll_cut, sig_cut, numbins, minval, maxval):
	#, norbits, photom_err
	print 'N=%s orbits, pme = %s mag'%(num_orbits, ph_unc.replace('p','.'))
	run, del_val, err_up, err_low, ll = read_logs(dirc, num_orbits, ph_unc)
	#print run, del_val, err_up, err_low, ll
	array=(del_val-1.5)*24.
	hist, bin_edges = numpy.histogram(array, bins=numbins, range=[minval,maxval], density=False)
	bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

	ax=subplot(221)
	plot(bin_centres, hist, lw=3, drawstyle='steps-mid', label='Test data')
	xticks([-24, -12, 0, 12, 24], rotation=45 )
	xlabel('reconstructed - true delay, hours')
	p0 = [20., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
	coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
	# Get the fitted curve
	hist_fit = gauss(bin_centres, *coeff)
	x=arange(minval,maxval,(maxval-minval)*1.e-3)
	plot(x, gauss(x,*coeff), 'r-', lw=2)
	print '\tmu    = %1.2e'%(coeff[1])
	print '\tsigma = %1.2e'%(coeff[2])
	title('Delay Reconsruction')
	text(0.97, 0.8,r'$\mu$=%1.2f hours'%(coeff[1]), ha='right', va='center', transform=ax.transAxes,fontsize=13)
	text(0.97, 0.7,r'$\sigma$=%1.2f hours'%(coeff[2]), ha='right', va='center', transform=ax.transAxes,fontsize=13)


	subplot(223)
	errorbar(ll,array,yerr=(err_low*24.,err_up*24.),fmt='.', ms=3)
	ylabel('reconstructed - true delay, hours')
	xlabel('log likelihood')
	plot([ll_cut,ll_cut],[axis()[2],axis()[3]], 'r--')
	xlim(0.9*max(ll_cut, min(ll)), 1.1*max(ll))

	xticks(rotation=45 )


	if(ll_cut>0.):
		#CUTS
		print '\tscs   = %d'%(len(ll[ ll>ll_cut]))
		cut2_list = [k for k in range(0,len(ll)) if (ll[k] > ll_cut and (abs(del_val[k]-1.5)*24.>sig_cut*coeff[2] or err_low[k]*24.>sig_cut*coeff[2] or err_up[k]*24.>sig_cut*coeff[2]))]
		print '\tscs   = %d'%(len(ll)-len(cut2_list))
		cut_tot_list = [k for k in range(0,len(ll)) if (ll[k] < ll_cut or (abs(del_val[k]-1.5)*24.>sig_cut*coeff[2] or err_low[k]*24.>sig_cut*coeff[2] or err_up[k]*24.>sig_cut*coeff[2]))]
		accept_tot_list = [k for k in range(0,len(ll)) if (ll[k] > ll_cut and (abs(del_val[k]-1.5)*24.<sig_cut*coeff[2] and err_low[k]*24.<sig_cut*coeff[2] and err_up[k]*24.<sig_cut*coeff[2]))]
		print '\tscs   = %d'%(len(ll)-len(cut_tot_list))

		#print axis()
		ax=subplot(222)
		hist, bin_edges = numpy.histogram(array[accept_tot_list], bins=numbins, range=[minval,maxval], density=False)
		plot(bin_centres, hist, lw=3, drawstyle='steps-mid', label='Test data')
		xticks([-24, -12, 0, 12, 24], rotation=45 )
		#xticks([-3, -2, -1, 0, 1, 2, 3], rotation=45 )
		xlabel('reconstructed - true delay, hours')
		p0 = [20., 0., 1.] # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
		coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
		# Get the fitted curve
		hist_fit = gauss(bin_centres, *coeff)
		x=arange(minval,maxval,(maxval-minval)*1.e-3)
		plot(x, gauss(x,*coeff), 'r-', lw=2)
		print '\tmu    = %1.2e'%(coeff[1])
		print '\tsigma = %1.2e'%(coeff[2])
		title('Delay Reconsruction With Cuts')
		text(0.97, 0.9,r'%d%s success'%(len(ll[accept_tot_list]),'%'), ha='right', va='center', transform=ax.transAxes,fontsize=13)
		text(0.97, 0.8,r'$\mu$=%1.2f hours'%(coeff[1]), ha='right', va='center', transform=ax.transAxes,fontsize=13)
		text(0.97, 0.7,r'$\sigma$=%1.2f hours'%(coeff[2]), ha='right', va='center', transform=ax.transAxes,fontsize=13)

		subplot(224)
		errorbar(ll[ll>ll_cut],array[ll>ll_cut],yerr=(err_low[ll>ll_cut]*24.,err_up[ll>ll_cut]*24.),fmt='.', ms=3)
		if(len(cut2_list)!=0):
			errorbar(ll[cut2_list], array[cut2_list], yerr=(err_low[cut2_list]*24.,err_up[cut2_list]*24.),  fmt='r.', ecolor='r', ms=3)
		xticks(rotation=45 )
		xlim(0.9*max(ll_cut, min(ll)), 1.1*max(ll))
		ylabel('reconstructed - true delay, hours')
		xlabel('log likelihood')
	suptitle('%s orbits, %s Photometric Uncertainty'%(num_orbits, ph_unc.replace('p','.')), fontsize=20)

	subplots_adjust(left=None, bottom=None, right=None, top=0.85, wspace=None, hspace=0.5)
	return coeff[1], coeff[2], len(ll[accept_tot_list])





def test(chain_file_name):
	samples = get_chain(chain_file_name)
	delay, d_mag, sigma, tau, avg_mag = get_parameter_samples(samples)
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



