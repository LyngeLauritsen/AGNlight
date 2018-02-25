import numpy as np
import scipy as sc
import scipy.optimize as opt
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from numpy.random import randint
import gc
#%%
runs = 100
#%%
def lognorm(x,mu,sigma):
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)

#fast_leave = 1193. #Minimum days to leave dust cloud
x_list = np.linspace(1,5000,5000)

mean = 200.
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma
log = []
for i in range(len(x_list)):
    log.append(lognorm(x_list[i],mu,sigma))
    
print np.nansum(log)*5000/5000.

#plt.figure()
#plt.plot(x_list,log)
#plt.show()
'''
plt.figure()
plt.plot(x_list,log)
plt.xlim([0,1000])
plt.show()
'''
#%%
data = np.loadtxt('NOVEMBER/Kelly-NGC3783K')
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')
#data[:,0] = data[:,0] - min(data[:,0])
#error[:,0] = float(int(error[:,0]))
#%%
print data[:,0]
print data[:,1]
#print error[:,0]
print error[:,1]
#%%
cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]),1)
print np.shape(cont_days)
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = np.nanmean(data[:,1])
data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0
#print log[:]
chi1 = 1e10000
#plt.figure()
for i in range(runs):
    print i
    for j in range(len(cont_days)):
        #print i, j/float(len(cont_days))
        data_comp[:,2] = 0
        #place1 = randint(0,len(cont_days))
        #place2 = randint(0,len(cont_days))
        #place3 = randint(0,len(cont_days))
        #place4 = randint(0,len(cont_days))
        #place5 = randint(0,len(cont_days))
        change = (-1)**randint(0,2)*cont[j,1]*random.random()*0.05
        #change1 = (-1)**randint(0,2)*cont[place1,1]*random.random()*0.05
        #change2 = (-1)**randint(0,2)*cont[place2,1]*random.random()*0.05
        #change3 = (-1)**randint(0,2)*cont[place3,1]*random.random()*0.05
        #change4 = (-1)**randint(0,2)*cont[place4,1]*random.random()*0.05
        #change5 = (-1)**randint(0,2)*cont[place5,1]*random.random()*0.05
        cont[j,1] += change
        #cont[place1,1] += change1
        #cont[place2,1] += change2
        #cont[place3,1] += change3
        #cont[place4,1] += change4
        #cont[place5,1] += change5        
        h = 0
        #print 'j = ', j
        for h in range(len(data_comp[:,0])):
            for k in range(len(cont_days)):
                if cont_days[k] < data_comp[h,0]:
                    data_comp[h,2] += cont[k,1]*log[abs(int((cont[k,0]-data_comp[h,0])))] #int((cont[k,0]-data_comp[h,0])/10.)
                    #print 'yes'
            #print 'k = ', k
        #print data_comp[:,2]
        chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2)
        #print data_comp[:,2]
        #print chi2, chi1
        if chi2 >= chi1:
            #print 'YES'
            cont[j,1] -= change
            #cont[place1,1] -= change1
            #cont[place2,1] -= change2
            #cont[place3,1] -= change3
            #cont[place4,1] -= change4
            #cont[place5,1] -= change5
        elif chi2 < chi1:
            #print chi2 - chi1
            chi1 = chi2
        #plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
        #plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
        #plt.scatter(cont[:,0],cont[:,1],color='g',s=0.1)
        #plt.ylim([4e-15,4e-14])
        #plt.yscale('log')
        #plt.show()
        gc.collect()
    gc.collect()

print cont[:,1]

np.savetxt('continuum3',cont)
np.savetxt('data_comp3',data_comp)
