from pylab import *
from analysis_library import *

def single_obs_view(fnm,chain_fnm):

	mag1, mag2, time_array, mag1_dat, mag2_dat = read_sim_hubble_light_curve(fnm)


	nwalkers     = 100
	n_iterations = 1000
	n_iteration_filter=600

	samples = get_chain(chain_fnm)
	print len(samples[:,[0,4,1,2,3]])
	figure()
	subplot(311)
	plot(samples[:,[0]],',')
	subplot(312)
	array_index = range(0,len(samples[:,[0]]))
	plot(samples[:,[0]][mod(array_index,1000)>900],',')
	plot([0,1000],[1.5,1.5],'k--')
	delay, d_mag, sigma, tau, avg_mag = get_parameter_samples(samples)

	subplot(521)
	plot(delay, ',')
	subplot(523)
	plot(d_mag, ',')
	subplot(525)
	plot(log10(sigma), ',')
	subplot(527)
	plot(log10(tau), ',')
	subplot(529)
	plot(avg_mag, ',')

	delay   = filter_samples(delay,   n_iterations, n_iteration_filter)
	d_mag   = filter_samples(d_mag,   n_iterations, n_iteration_filter)
	sigma   = filter_samples(sigma,   n_iterations, n_iteration_filter)
	tau     = filter_samples(tau,   n_iterations, n_iteration_filter)
	avg_mag = filter_samples(avg_mag,   n_iterations, n_iteration_filter)

	subplot(522)
	plot(delay, ',')
	subplot(524)
	plot(d_mag, ',')
	subplot(526)
	plot(log10(sigma), ',')
	subplot(528)
	plot(log10(tau), ',')
	subplot(5,2,10)
	plot(avg_mag, ',')

	figure(3)
	subplot(321)
	dt_min = min(delay)-1.
	dt_max = max(delay)+1.
	delta_dt = (dt_max-dt_min)/1000
	p_dt,  x_dt  = hist_param(delay, arange(dt_min,dt_max,delta_dt))
	plot(x_dt,p_dt,'-', lw=2, label='%s'%chain_fnm, drawstyle='steps')
	#xlim(1.4,2.4)
	dm_min = -5.
	dm_max=5.
	delta_dm = 0.001
	p_dm,  x_dm  = hist_param(d_mag, arange(dm_min,dm_max,delta_dm))
	subplot(324)
	plot(x_dm,p_dm,'-', lw=2, label='%s'%chain_fnm, drawstyle='steps')
	xlim(0.25,0.35)
	subplot(323)
	plot(delay,d_mag,',')

	subplot(322)
	plot(time_array, mag1_dat,'b.', mec='b',         ms=5)
	plot(time_array-1.5, mag2_dat-0.3, 'r.', mec='r', ms=5)
	subplot(325)
	plot(time_array, mag1_dat,'b.', mec='b',         ms=5)
	plot(time_array+2.9, mag2_dat-0.65, 'r.', mec='r', ms=5)

	subplot(326)
	plot(time_array, mag1_dat,'b.', mec='b',         ms=5)
	plot(time_array-2.995, mag2_dat-0.3605, 'r.', mec='r', ms=5)

	filtered_samples = samples[mod(array_index,1000)>200]
	#fig= triangle.corner(samples[:,[0,4,1,2,3]], labels=["delay", "avg_mag", "$\Delta$mag", "$\sigma$", r"$\tau$"])
	fig= triangle.corner(filtered_samples[:,[0,1,4,2,3]], labels=["delay", "$\Delta$m", r"$\langle m \rangle$", "$\sigma$", r"$\tau$"])

	a1 = axes([.55, .83, .4, .14], axisbg='none')
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3)
	errorbar(time_array, mag2_dat ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3 )
	title('Light Curves')
	#xlabel('time, days')
	ylabel('Magnitude, arb. u.')
	a1.set_xticklabels('')
	setp(a1, yticks=arange(18.8,19.5,0.2))
	a2 = axes([.55, .65, .4, .14], axisbg='none')
	dt = x_dt[argmax(p_dt)]
	dm = x_dm[argmax(p_dm)]
	errorbar(time_array, mag1_dat,fmt='b+', yerr=0.02*ones(len(time_array)), mec='b', ms=3, label='light curve 1')
	errorbar(time_array-dt, mag2_dat-dm ,fmt='rx', yerr=0.02*ones(len(time_array)), mec='r', ms=3, label='light curve 2' )
	legend(loc=1)
	#title('Light Curves')
	xlabel('time, days')
	ylabel('Magnitude, arb. u.')
	setp(a2, yticks=arange(18.8,19.21,0.1))
	show()
	exit()

	figure(4)
	subplot(331)
	sigma_min = min(log10(sigma))
	sigma_max=max(log10(sigma))
	delta_sigma = (sigma_max-sigma_min)/100
	p_sigma,  x_sigma  = hist_param(log10(sigma), arange(sigma_min,sigma_max,delta_sigma))
	plot(x_sigma,p_sigma,'-', lw=2, label='%s'%chain_fnm, drawstyle='steps')


	subplot(335)
	tau_min = min(log10(tau))
	tau_max=max(log10(tau))
	delta_tau = (tau_max-tau_min)/100
	p_tau,  x_tau  = hist_param(log10(tau), arange(tau_min,tau_max,delta_tau))
	plot(x_tau,p_tau,'-', lw=2, label='%s'%chain_fnm, drawstyle='steps')
	subplot(334)
	hist(tau)
	#subplot(323)
	#plot(time_array, mag1_dat,'b,')
	#plot(time_array-1.5, mag2_dat-0.3, 'r.')

	#nwalkers     = 100
	#n_iterations = 1000



	show()

def ll_plot():
	pdfp=PdfPages('systematics_ll.pdf')
	pdfp2=PdfPages('systematics_smy.pdf')
	pdfp3=PdfPages('ana_ll_example.pdf')
	#pdfp=PdfPages('systematics_ll_0p01.pdf')
	#dir1 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_08_hubble_sim_runs/'
	#dir2 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_09_hubble_sim_runs/' 
	#dir3 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_10_hubble_sim_runs/'
	dir1 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_15_hubble_sim_runs/'
	sig_cut = 3.0
	'''
	fig=figure()
	m_40_0p01, s_40_0p01, scs_40_0p01 = ana_ll(dir3, 8000, 9100, 1700., sig_cut, 50, -2., 2., 40, 0.01)
	pdfp.savefig(fig)
	fig=figure()
	m_80_0p01, s_80_0p01, scs_80_0p01 = ana_ll(dir3, 9000, 9100, 1800., sig_cut, 50, -2., 2., 80, 0.01)
	show()
	pdfp.savefig(fig)
	pdfp.close()
	exit()
	'''

	num_orbits_list=['40', '50', '60', '70', '80', '90', '100']
	pu='0p02'
	m_0p02=[]
	s_0p02=[]
	scs_0p02=[]
  	for norb in num_orbits_list:
  	  fig=figure()
	  m, s, scs = ana_ll(dir1, norb, pu, 0.01, sig_cut, 250, -30., 30.)
	  m_0p02.append(m)
	  s_0p02.append(s)
	  scs_0p02.append(scs)
	  pdfp.savefig(fig)
	scs_0p02_err_up=[]
	scs_0p02_err_lo=[]
	for k in range(0,len(scs_0p02)):
		l,u = binom_interval(scs_0p02[k], 100, confint=0.68)
		scs_0p02_err_lo.append(scs_0p02[k]-l)
		scs_0p02_err_up.append(u-scs_0p02[k])

	pu='0p05'
	m_0p05=[]
	s_0p05=[]
	scs_0p05=[]
  	for norb in num_orbits_list:
  	  fig=figure()
	  m, s, scs = ana_ll(dir1, norb, pu, 0.01, sig_cut, 100, -30., 30.)
	  m_0p05.append(m)
	  s_0p05.append(s)
	  scs_0p05.append(scs)
	  pdfp.savefig(fig)
	  if(norb=='60'):
		print 'SAVING EXAMPLE FIGURE'
		pdfp3.savefig(fig)
		pdfp3.close()
	scs_0p05_err_up=[]
	scs_0p05_err_lo=[]
	for k in range(0,len(scs_0p05)):
		l,u = binom_interval(scs_0p05[k], 100, confint=0.68)
		scs_0p05_err_lo.append(scs_0p05[k]-l)
		scs_0p05_err_up.append(u-scs_0p05[k])
	pu='0p10'
	m_0p10=[]
	s_0p10=[]
	scs_0p10=[]
  	for norb in num_orbits_list:
  	  fig=figure()
	  m, s, scs = ana_ll(dir1, norb, pu, 0.01, sig_cut, 50, -30., 30.)
	  m_0p10.append(m)
	  s_0p10.append(s)
	  scs_0p10.append(scs)
	  pdfp.savefig(fig)
	scs_0p10_err_up=[]
	scs_0p10_err_lo=[]
	for k in range(0,len(scs_0p10)):
		l,u = binom_interval(scs_0p10[k], 100, confint=0.68)
		scs_0p10_err_lo.append(scs_0p10[k]-l)
		scs_0p10_err_up.append(u-scs_0p10[k])
		
	fig=figure()
	orbs = map(int, num_orbits_list)
	subplot(211)
	errorbar(orbs, m_0p10, yerr=s_0p10, fmt='r_', label='0.10 mag', capsize=0, elinewidth=30, ms=0, ecolor='r', mew=3)
	errorbar(orbs, m_0p05, yerr=s_0p05, fmt='b_', label='0.05 mag', capsize=0, elinewidth=20, ms=0, ecolor = 'b', mew=3)
	errorbar(orbs, m_0p02, yerr=s_0p02, fmt='k_', label='0.02 mag', capsize=0, elinewidth=10, ms=0, ecolor='k', mew=3)
	#errorbar(orbs, m_0p10, yerr=s_0p10, fmt='r_', label='0.10 mag', capsize=0, elinewidth=30, ms=0, ecolor='0.7', mew=3)
	#errorbar(orbs, m_0p05, yerr=s_0p05, fmt='b_', label='0.05 mag', capsize=0, elinewidth=20, ms=0, ecolor = '0.4', mew=3)
	#errorbar(orbs, m_0p02, yerr=s_0p02, fmt='k_', label='0.02 mag', capsize=0, elinewidth=10, ms=0, ecolor='0.0', mew=3)
	xlim(35.,125.)
	ylim(-5.9,5.9)
	ylabel('Light Curve Delay\nResolution, hours')
	#xlabel('Number of Orbits')
	grid(True)
	legend(loc=1, title='Photometric\nUncertainty')
	tick_params(axis='x', labelbottom='off')
	yticks(arange(-5.,5.1,1.))
	subplot(212)
	#print binom_interval(50, 100, confint=0.68)
		
	errorbar(orbs,scs_0p10, yerr=(scs_0p10_err_lo, scs_0p10_err_up), fmt='r_', color='r', capsize=0, elinewidth=30, label='0.10 mag', ms=0)
	errorbar(orbs,scs_0p05, yerr=(scs_0p05_err_lo, scs_0p05_err_up), fmt='b_', color='b', capsize=0, elinewidth=20, label='0.05 mag', ms=0)
	errorbar(orbs,scs_0p02, yerr=(scs_0p02_err_lo, scs_0p02_err_up), fmt='k_', color='k', capsize=0, elinewidth=10, label='0.02 mag', ms=0)
	legend(loc=1, title='Photometric\nUncertainty')
	ylim(40.,100.)
	xlim(35.,125.)
	ylabel('Reconstruction\nSuccess Rate, %')
	xlabel('Number of Orbits')
	grid(True)
	#show()
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.1)
	pdfp.savefig(fig)
	pdfp2.savefig(fig)
	pdfp.close()
	pdfp2.close()
	exit()
#fnm       = '../outputs/sim_data_hubble_test_2014_05_28_16:46.npz'
#chain_fnm = '../outputs/chain_samples_hubble_test_2014_05_28_16:46.npz'

#fnm       = '../outputs/sim_data_hubble_test_2014_05_29_10:05.npz'
#chain_fnm = '../outputs/chain_samples_hubble_test_2014_05_29_10:05.npz'

#dir1 = '/data2/fot_archived_outputs_and_logs/logs/2014_04_15_hubble_sim_runs/'

fnm	  = '/data2/fot_archived_outputs_and_logs/outputs/2014_04_15_hubble_sim_runs/sim_data_hubble_sim_90_orbits_0p02_pu_r97_2014_04_16_23:21.npz'
chain_fnm = '/data2/fot_archived_outputs_and_logs/outputs/2014_04_15_hubble_sim_runs/chain_samples_hubble_sim_90_orbits_0p02_pu_r97_2014_04_16_23:21_2014_04_16_23:21.npz'

#single_obs_view(fnm,chain_fnm)

ll_plot()
#ll_plot2()
exit()

