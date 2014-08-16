from analysis_library import *

fnm='chain_samples_tdc1_rung0_double_pair20_dt_2014_08_08_14:02.npz'

chain = get_chain('../outputs/'+fnm)

#print chain

#print len(chain[:,1])
#print len(chain[:,1])/100/4000
figure(figsize=(10,14))
for n in range(0,8):
  subplot(8,2,2*n+1)
  for k in range(0,100):
    plot((chain[k*4000:(k+1)*4000,n]))
  subplot(8,2,2*n+2)
  for k in range(0,100):
    plot((chain[k*4000+2000:(k+1)*4000,n]))
subplots_adjust(bottom=0.03, top=0.97)

print shape(chain[:,0]), (chain[:,0])
print shape(((chain[:,0].reshape((100, 4000)))[:,2000:4000])), ((chain[:,0].reshape((100, 4000)))[:,2000:4000])
#print shape(reshape(chain[:,0], (-1, 100))[2000:4000]), reshape(chain[:,0], (-1, 100))[2000:4000]
#print shape(ravel(reshape(chain[:,0], (-1, 100))[2000:4000,:]))
#print ravel(reshape(chain[:,0], (-1, 100))[2000:4000,:])
#print ravel(reshape(chain[:,0], (-1, 100))[2000:4000,:])
#print ravel(reshape(chain[:,0], (-1, 100))[2000:4000,:])
figure()
subplot(411)
plot((chain[:,0]))
subplot(412)
plot(ravel((chain[:,0].reshape((100, 4000)))[:,2000:4000]))
subplot(413)
print len(sort(ravel((chain[:,0].reshape((100, 4000)))[:,2000:4000])))
plot(sort(ravel((chain[:,0].reshape((100, 4000)))[:,2000:4000])), arange(200000, dtype=float64)/200000., ',')
subplot(414)
x = sort(ravel((chain[:,0].reshape((100, 4000)))[:,2000:4000]))
c = arange(200000, dtype=float64)/200000.
N=len(x)
print 'calculating minimum'
alpha = 0.68
'''
figure()
for n in range(0,int((1.-alpha)*200000))[::2000]:
  #print n, len(c[n:]), c[n], max(c[n:]-c[n]),  (c[n:]-c[n])[c[n:]-c[n]<0.68], (x[n:]-x[n])[c[n:]-c[n]<0.68]
  N = len((x[n:]-x[n])[c[n:]-c[n]<0.68])
  plot((x[n:]-x[n])[c[n:]-c[n]<alpha],(c[n:]-c[n])[c[n:]-c[n]<alpha],',')
  print n, max((x[n:]-x[n])[c[n:]-c[n]<alpha]), x[n], c[n], c[numpy.max(numpy.where(abs(c[n:]-c[n])<=alpha))]
  #,  len(c[n:]), c[n], max(c[n:]-c[n]),(x[n:]-x[n])[c[n:]-c[n]<0.68],(x[n:]-x[n])[c[n:]-c[n]<0.68]
  #print (numpy.where(abs(c[n:]-c[n]-0.001)<=alpha))
  #print numpy.max(numpy.where(abs(c[n:]-c[n]-0.001)<=alpha))
print 'DONE WITH LOOP'
'''
newlist = array([ max((x[n:]-x[n])[c[n:]-c[n]<alpha]) for n in range(0,int((1.-alpha)*200000))[::200]])
#newlist = [ x[numpy.max(numpy.where(abs(c[n:]-c[n])<=alpha))] for n in range(0,int((1.-alpha)*200000))[::20]]
#newlist = [ x[numpy.max(numpy.where(abs(c[n:]-c[n]-0.001)<=alpha))]-x[n] for n in range(0,int((1.-alpha)*200000))[::20]]
print 'done calculating minimum'
print numpy.where((newlist == min(newlist)))[0]
print x[numpy.where((newlist == min(newlist)))[0]*200][0],min(newlist)


figure(figsize=(8,8))
for n in range(0,8):
  subplot(3,3,n+1)
  hist(ravel((chain[:,n].reshape((100, 4000)))[:,2000:4000]), bins=50)
  if(n==0):
    plot([x[numpy.where((newlist == min(newlist)))[0]*200][0],x[numpy.where((newlist == min(newlist)))[0]*200][0]], [0,25000], 'r-')
    plot([x[numpy.where((newlist == min(newlist)))[0]*200][0]+min(newlist),x[numpy.where((newlist == min(newlist)))[0]*200][0]+min(newlist)], [0,25000], 'r-')
    m = x[numpy.where((newlist == min(newlist)))[0]*200][0]+min(newlist)/2.
    plot([m,m], [0,25000], 'r-')
subplots_adjust(bottom=0.03, top=0.97)
show()