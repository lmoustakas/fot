from pylab import *

def read_data(fname):
  jd=[]
  mag1=[]
  e1=[]
  mag2=[]
  e2=[]
  mag3=[]
  e3=[]
  mag4=[]
  e4=[]
  lc=0
  for line in file(fname):
    #print line
    lc+=1
    if(lc>2):
      line.split()[1]
      jd.append(float(line.split()[0]))
      mag1.append(float(line.split()[1]))
      e1.append(float(line.split()[2]))
      mag2.append(float(line.split()[3]))
      e2.append(float(line.split()[4]))
      mag3.append(float(line.split()[5]))
      e3.append(float(line.split()[6]))
      mag4.append(float(line.split()[7]))
      e4.append(float(line.split()[8]))
    
  return array(jd), array(mag1),array(e1),array(mag2),array(e2),array(mag3),array(e3),array(mag4),array(e4)

#jd, mag1, e1, mag2, e2, mag3, e3, mag4, e4 = read_data('../data/cosmograil/RXJ1131_Tewes2013.rdb') 
#print mag1
#errorbar(jd, mag1, yerr=e1, fmt='.')
#errorbar(jd, mag2, yerr=e2, fmt='.')
#errorbar(jd, mag3, yerr=e3, fmt='.')
#errorbar(jd, mag4, yerr=e4, fmt='.')
#show()