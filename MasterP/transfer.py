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
from multiprocessing import Process

'''The Continuum and Observed Light Curves are loaded'''
CLC = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1')
OLC = np.loadtxt('NOVEMBER/NOV-NGC3783-K')
OLC_error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')

#print OLC_error[:,1]
#print OLC[:,1]

step = 5.

mean = 140. #Transfer function
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma
normalize = (sigma*np.sqrt(2*np.pi))

'''Starting parameters are made for the transfer functions'''
mean_hot = 130. #Transfer function
sigma_hot = 0.8
mu_hot = np.log(mean_hot) - 1/2.*sigma_hot
normalize_hot = sigma_hot*np.sqrt(2*np.pi)

mean_cold = 200. #Transfer function
sigma_cold = 0.8
mu_cold = np.log(mean_cold) - 1/2.*sigma_cold
normalize_cold = sigma_cold*np.sqrt(2*np.pi)

vars = [mean_hot,sigma_hot,mu_hot,normalize_hot,mean_cold,sigma_cold,mu_cold,normalize_cold]

delay_days = np.linspace(0.01,1500,1500) #X-axis of transfer function. 1500 days are chosen as no observed light curves has this length.

#print delay_days

def create_lognorm(x,mu,sigma,normalize):
    '''Defines the transfer function'''
    sigma_array = np.zeros((np.shape(delay_days))) + float(sigma)
    mu_array = np.zeros((np.shape(delay_days))) + float(mu)
    exp_term = np.zeros((np.shape(delay_days)))
    exp_term = -((np.log(x)-mu_array)**2/(2*sigma_array**2))
    front = 1/(x*normalize)
    return front*np.exp(exp_term)

lognorm_hot = create_lognorm(delay_days,mu_hot,sigma_hot,normalize_hot) #creating the transfer function for testing.
lognorm_cold = create_lognorm(delay_days,mu_cold,sigma_cold,normalize_cold) #creating the transfer function for testing.
lognorm = create_lognorm(delay_days,mu,sigma,normalize) #creating the transfer function for testing.

lognorm_sum = lognorm_hot + lognorm_cold

plt.figure()
plt.scatter(delay_days,lognorm_hot,color='b')
plt.scatter(delay_days,lognorm_cold,color='yellow')
plt.scatter(delay_days,lognorm,color='g')
plt.scatter(delay_days,lognorm_sum,color='r')
plt.show(block=False)


def transfer(cont_days,cont,data_comp,x,mu,sigma,normalize):
    '''Transfer function array'''
    log_norm_transfer = create_lognorm(x,mu,sigma,normalize)
    transfer = np.zeros((len(cont_days),len(data_comp[:,0]))) #Array with the transferred light in any combination of the continuum and the observed.
    for k in range(len(cont_days)-1):
        '''Creating the array for the amount of transmitted light'''
        for h in range(len(data_comp[:,0])):
            if cont_days[k] < data_comp[h,0]:
                transfer[k,h] = log_norm_transfer[abs(int((cont[k,0]-data_comp[h,0])))]
    return transfer.transpose()

def model_data(cont,data_comp,x,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold):
    '''The simulated observed light curve as a result of the continuum light curve and the transfer function'''

    '''The creation of the transfer function arrays'''
    transfer_array_hot = transfer(cont[:,0],cont,data_comp,x,mu_hot,sigma_hot,normalize_hot)
    transfer_array_cold = transfer(cont[:,0],cont,data_comp,x,mu_cold,sigma_cold,normalize_cold)

    cont_flux = np.array([cont[:,1],]*len(data_comp[:,0]))

    '''The simulated observed light curves for the two transfer functions'''
    sim_OLC_hot = np.sum(step*cont_flux*transfer_array_hot,axis=1)
    sim_OLC_cold = np.sum(step*cont_flux*transfer_array_cold,axis=1)

    '''The total simulated observed light curves'''
    sim_OLC = sim_OLC_hot + sim_OLC_cold

    return sim_OLC

model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)

chi = np.nansum((OLC[:,1] - model_qual)**2)

def residual(vars,delay_days,OLC,OLC_error,CLC):
    '''Starting parameters are made for the transfer functions'''
    mean_hot = vars[0] #Transfer function
    sigma_hot = vars[1]
    mu_hot = vars[2]
    normalize_hot = vars[3]

    mean_cold = vars[4] #Transfer function
    sigma_cold = vars[5]
    mu_cold = vars[6]
    normalize_cold = vars[7]

    model = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)

    res = ((OLC[:,1] - model)/OLC_error[:,1])
    return res[np.isfinite(res)]

out = leastsq(residual, vars, args=(delay_days,OLC,OLC_error,CLC))

print vars
print out

vars = out[0]

mean_hot = vars[0] #Transfer function
sigma_hot = vars[1]
mu_hot = vars[2]
normalize_hot = vars[3]

mean_cold = vars[4] #Transfer function
sigma_cold = vars[5]
mu_cold = vars[6]
normalize_cold = vars[7]

lognorm_hot = create_lognorm(delay_days,mu_hot,sigma_hot,normalize_hot) #creating the transfer function for testing.
lognorm_cold = create_lognorm(delay_days,mu_cold,sigma_cold,normalize_cold) #creating the transfer function for testing.
lognorm = create_lognorm(delay_days,mu,sigma,normalize) #creating the transfer function for testing.

lognorm_sum = lognorm_hot + lognorm_cold

plt.figure()
plt.scatter(delay_days,lognorm_hot,color='b')
plt.scatter(delay_days,lognorm_cold,color='yellow')
plt.scatter(delay_days,lognorm,color='g')
plt.scatter(delay_days,lognorm_sum,color='r')
plt.show(block=False)

plt.figure()
'''Testing transfer functions'''
plt.scatter(OLC[:,0],OLC[:,1],color='r')
plt.scatter(OLC[:,0],model_qual,color='g')
plt.ylim([4e-15,1.2e-14])
plt.show()


'''
print chi

for i in range(3000):
    for j in range(100):
        mean_hot1 = mean_hot + (-1)**randint(0,2)*0.1
        model_qual1 = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual1)**2)
        print chi1, chi
        if chi1 < chi:
            chi = chi1
            mean_hot = mean_hot1
            print j

        mean_cold1 = mean_cold + (-1)**randint(0,2)*0.1
        model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual)**2)
        if chi1 < chi:
            chi = chi1
            mean_cold = mean_cold1
            print j

        sigma_hot1 = sigma_hot + (-1)**randint(0,2)*0.1
        model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual)**2)
        if chi1 < chi:
            chi = chi1
            sigma_hot = sigma_hot1
            print j

        sigma_cold1 = sigma_cold + (-1)**randint(0,2)*0.1
        model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual)**2)
        if chi1 < chi:
            chi = chi1
            sigma_cold = sigma_cold1
            print j

        normalize_hot1 = normalize_hot + (-1)**randint(0,2)*0.1
        model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual)**2)
        if chi1 < chi:
            chi = chi1
            normalize_hot = normalize_hot1
            print j

        normalize_cold1 = normalize_cold + (-1)**randint(0,2)*0.1
        model_qual = model_data(CLC,OLC,delay_days,mu_hot,mu_cold,sigma_hot,sigma_cold,normalize_hot,normalize_cold)
        chi1 = np.nansum((OLC[:,1] - model_qual)**2)
        if chi1 < chi:
            chi = chi1
            normalize_cold = normalize_cold1
            print j

    lognorm_hot = create_lognorm(delay_days,mu_hot,sigma_hot,normalize_hot) #creating the transfer function for testing.
    lognorm_cold = create_lognorm(delay_days,mu_cold,sigma_cold,normalize_cold) #creating the transfer function for testing.
'''
