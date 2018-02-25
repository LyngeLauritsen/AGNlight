import numpy as np
import scipy as sc
import scipy.optimize as opt
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from numpy.random import randint
import math
from scipy.optimize import minimize as mini
from scipy.optimize import least_squares as ls
from scipy.optimize import leastsq
import gc

runs = 30

def lognorm(x,mu,sigma):
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)

x_list = np.linspace(0.01,1500,1500)

mean = 140.
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma
log = []
for i in range(len(x_list)):
    log.append(lognorm(x_list[i],mu,sigma))
    
print np.nansum(log)*1500/5000.

data = np.loadtxt('NOVEMBER/Kelly-NGC3783K')
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')

print data[:,0]
print data[:,1]

print error[:,1]
step = 10
cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]),step)
print np.shape(cont_days)
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = np.nanmean(data[:,1])
data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0

chi1 = 1e10000
max_jump = 0.01

for i in range(runs):
    #if i == int(0.4*runs):
    #    max_jump = 0.01
    #if i == int(0.8*runs):
    #    max_jump = 0.005
    for i1 in range(runs):
        print i, i1/float(runs)
        for j in range(len(cont_days)-2):
            #print i, j/float(len(cont_days))
            data_comp[:,2] = 0
            change = (-1)**randint(0,2)*cont[j,1]*random.random()*max_jump
            cont[j,1] += change
            h = 0
            for h in range(len(data_comp[:,0])):
                for k in range(len(cont_days)-1):
                    if cont_days[k] < data_comp[h,0]:
                        data_comp[h,2] += step*np.mean((cont[k,1],cont[k+1,1]))*log[abs(int((cont[k,0]-data_comp[h,0])))]
            #print data_comp[:,2]
            chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2)
            chi_test = random.random()
            #jump = abs((cont[j-2,1]+cont[j-1,1]+cont[j+1,1]+cont[j+2,1])/4. - cont[j,1])
            if chi2 <= chi1: # or chi_test < 0.02: # and jump < 0.05*cont[j,1]
                chi1 = chi2
            else:
                cont[j,1] -= change
            gc.collect()
    plt.figure()
    plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    plt.scatter(cont[:,0],cont[:,1],color='g')
    plt.ylim([4e-15,1.2e-14])
    plt.show(block=False)
    gc.collect()

print cont[:,1]

np.savetxt('continuum3',cont)
np.savetxt('data_comp3',data_comp)
