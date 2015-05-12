from pylab import *

def break_into_seasons(fname):
    lc=0
    l1=[]
    l2=[]
    t_array=[]
    t_prev=0.
    prev_line=[]
    season=1
    fname_out = fname.split('.')[0]+'_season%d'%season+'.rdb'
    f_out = open(fname_out,'w')
    for line in file(fname):
      lc+=1
      if(lc==1): l1=line
      if(lc==2): l2=line
      if(lc<3): print line
      
      if(lc>=3):
	 t=line.split()[0]
	 t_array.append(float(t))
	 dt = float(t)-float(t_prev)
	 if(dt>60.):
	   f_out.close()
	   fname_out = fname.split('.')[0]+'_season%d'%season+'.rdb'
	   f_out = open(fname_out,'w')
	   f_out.write(l1)
	   f_out.write(l2)
	   if(season>0):
	    if(len(prev_line)!=0): print prev_line
	    print '\nseason', season, ' dt=',dt
	    print '------------------------------------------------------------------------------------------------------------------'
	    print line
	   season+=1
	 t_prev = t
	 prev_line = line
	 f_out.write(line)
    
    f_out.close()
    #plot(t_array,'.')
    #show()
#break_into_seasons('RXJ1131_Tewes2013.rdb')
break_into_seasons('J1206_Eulaers2013.rdb')
break_into_seasons('J1001_Rathnakumar2013.rdb')
break_into_seasons('HS2209_Eulaers2013.rdb')
break_into_seasons('HE0435_Courbin2011.rdb')



