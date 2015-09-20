import numpy as np
fout = open('fot_lcogt_runs.txt', 'w')
for delay in np.arange(-2.0,2.1, 0.1):
	tag = 'lcogt_dAB_dt_%1.1f'%(delay)
	outdir = '/disk4/romerowo/fot_lcogt_outputs/20150919/'
	com = './fot_lcogt_delay.py -i /nisushome/romerowo/lcolens/python/arw/he0435-1223_lcogt_magnitudes.dat -im1 A -im2 B -dtp %1.1f -dtpmin %1.1f -dtpmax %1.1f -sp 0.25 -spmin 0.01 -spmax 1.0 -tp 160. -tpmin 5. -tpmax 1000. -o %s -od %s > /disk4/romerowo/fot_lcogt_outputs/20150919/%s.log'%(delay, delay-0.1, delay+0.1, tag, outdir, tag)
	fout.write(com+'\n')
fout.close()
