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
subplot(311)
plot((chain[:,0]))
subplot(312)
plot(ravel((chain[:,0].reshape((100, 4000)))[:,2000:4000]))
subplot(313)
plot((chain[2000:,0]))
#show()

figure(figsize=(8,8))
for n in range(0,8):
  subplot(3,3,n+1)
  hist(ravel((chain[:,n].reshape((100, 4000)))[:,2000:4000]), bins=100)
subplots_adjust(bottom=0.03, top=0.97)

figure(figsize=(8,8))
for n in range(0,8):
  subplot(3,3,n+1)
  hist(ravel((chain[:,n].reshape((100, 4000)))[:,2000:4000]), bins=100)
subplots_adjust(bottom=0.03, top=0.97)
show()