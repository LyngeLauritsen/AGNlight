import numpy as np
from numpy import *
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
import colorednoise as cn
from numpy.fft import ifft, fftfreq


runs = 3000 #Number of runs (arbitrary due to the gradual updates, would probably take close to a month to run)
runs1 = 100
days_before = 700. #How far to extend the lightcurve back (days)
days_after = 100. #How far to extend the lightcurve forward (days)
step = 5. #Timestep (days)
slopeaim = -2.5 #Slope that is targeted
slopeallow = -0.05
slope = 0
mean = 140. #Transfer function
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma
kelly_repeat = 20
dd_change = 1.1

# There are signal processing package in scipy that defines filters, prob. can find lognormal in there.
def lognorm(x,mu,sigma):
    '''Defines the transfer function'''
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)

def colour(alpha,length):
    return cn.powerlaw_psd_gaussian(alpha,length)

x_list = np.linspace(0.01,1500,1500) #X-axis of transfer function



log_norm_transfer = []
for i in range(len(x_list)):
    '''Creates the list used to identify the degree of transfered light'''
    log_norm_transfer.append(lognorm(x_list[i],mu,sigma))
#Suggest:
# log = lognorm(np.arange(len(x_list)))

print np.nansum(log_norm_transfer)*1500/1500.

data = np.loadtxt('NOVEMBER/Kelly-NGC3783K') #'NOVEMBER/NOV-NGC3783-K')
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt') #Load data

#plt.figure()
#plt.scatter(data[:,0],data[:,1])
#plt.ylim([6e-15,1.1e-14])

plt.show()

dd_ddt_data = []
for j in range(len(data[:,0])-2):
    '''Purely a test to see the largest double derivative of the observed light curve'''
    dd_ddt_data.append(abs(((data[j+1,1] - data[j,1])/(data[j+1,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+1,0] - data[j-1,0])/2.))

#
#    dd_ddt_data.append(abs(((data[j+4,1] - data[j,1])/(data[j+4,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+4
print np.nanmax(dd_ddt_data)

#Have a look at the numpy/scipy Runge-Kutta implementation
def num_deriv(data):
    dd_flux = abs(np.diff(data[:,1])/np.diff(data[:,0]))
    dd_time = data[:-1,0] + np.diff(data[:,0])/2.
    dd_data = np.array((dd_time,dd_flux)).transpose()
    dd_ddt_flux = abs(np.diff(dd_data[:,1])/np.diff(dd_data[:,0]))/2.
    dd_ddt_time = dd_data[:-1,0] + np.diff(dd_data[:,0])/2.0
    dd_ddt_data = np.array((dd_ddt_time,dd_ddt_flux)).transpose()
    return dd_ddt_data

dd_ddt_data1 = num_deriv(data)
print np.nanmax(dd_ddt_data1)


dd_ddt_allowed = np.nanmax(dd_ddt_data)*dd_change #Identifies the maximal double derivative allowed in the continuum light curve

'''cont denotes the model for the continuum light curve, whereas data_comp denotes the array giving indication of the
observed light curve and how the continuum light curve fitted through the transfer function compares'''
cont_days = np.arange(min(data[:,0])-days_before,max(data[:,0]+days_after),step) #Defines the spacing of the light curve
#cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]+100),step) #Defines the spacing of the light curve
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = colour(-2.5,len(cont[:,1]))*1e-18
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1')
day = cont_days[0]

data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0
#data_comp = np.loadtxt('CONTINUUM/NGC3783-data_comp-slope-Kelly-2_5-1')
#print data[:,2]

def transfer(cont_days,data_comp,log_norm_transfer):
    transfer = zeros((len(cont_days),len(data_comp[:,0]))) #Array with the transferred light in any combination of the continuum and the observed.
    for k in range(len(cont_days)-1):
        '''Creating the array for the amount of transmitted light'''
        for h in range(len(data_comp[:,0])):
            if cont_days[k] < data_comp[h,0]:
                transfer[k,h] = log_norm_transfer[abs(int((cont[k,0]-data_comp[h,0])))]
    return transfer.transpose()

transfer_array = transfer(cont_days,data_comp,log_norm_transfer)

flux1 = cont[:,1] #Series of items for the Kelly 2009
time1 = cont[:,0]
sigma1 = cont[:,2] + np.nanmean(error[:,1])
time = np.insert(time1,0,0)
sigma = np.insert(sigma1,0,np.mean(sigma1))
flux = np.insert(flux1,0,np.mean(flux1))

print shape(flux1)

freq = []

for i in range(int(len(cont_days)/2.)):
    '''The frequency determination for the PSD'''
    freq.append((i+1)/(60.*60.*24*float(cont_days[len(cont_days)-1] - cont_days[0])))

def PSD(freq,transfer_array,cont):
    '''Code designed to determine the Power Spectral Density Slope'''
    cont1 = cont[:,1] - np.nanmean(cont[:,1]) #Removing zero frequency

    cont_zero_freq = transpose(np.array([cont1,]*len(freq)))
    days_cont = transpose(np.array([cont[:,0],]*len(freq)))
    freq_array = np.array([freq,]*len(cont[:,0]))

    F_N_v = np.sum(cont_zero_freq*np.cos(2*np.pi*freq_array*60.*60.*24*days_cont),axis=0)**2 \
    + np.sum(cont_zero_freq*np.sin(2*np.pi*freq_array*60.*60.*24*days_cont),axis=0)**2

    #print len(cont[:,0]), F_N_v
    P_v = 2*60.*60.*24*float(cont[:,0][len(cont[:,0])-1] - cont[0,0]) \
    /(np.nanmean(cont[:,1])**2*len(cont[:,0])**2)*F_N_v #Finding PSD

    slope = np.polyfit(np.log10(freq),np.log10(P_v),1)[0] #Finding PSD slope

    return slope, P_v

def model_data(cont,data_comp,transfer_array):
    '''The simulated observed light curve as a result of the continuum light curve and the transfer function'''

    cont_flux = array([cont[:,1],]*len(data_comp[:,0]))

    return sum(step*cont_flux*transfer_array,axis=1)

model3 = cont[:,1]
tau = 800. #The Kelly relaxation time
sigma_tot = np.nanmean(error[:,1]) #Kelly sigma
chi1 = 1e100 #To ensure that the first change is always accepted to get it started
max_jump = 1 #0.01 #Maximal allowed jump relative to the value at point

slopeolder1 = 0 #Temporary constant

F_N_v = np.zeros((len(freq)))
P_v = 0
print F_N_v

for i in range(runs):
    '''The MCMC'''
    for i1 in range(runs1):
        '''Doing 2 loops to allow regular progress reports'''
        cont_real = cont
        cont_real[:,1] = ifft(cont[:,1]).real
        plt.figure()
        plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
        plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
        plt.scatter(cont_real[:,0],cont_real[:,1],color='g')
        plt.ylim([4e-30,1.2e5])
        plt.yscale('log')
        plt.show()
        
        #print i, i1/float(runs1), slopeolder1
        for j in range(1):
            '''Looping over the course of the CLC'''
            cont_save = cont
            F_N_v = np.zeros((len(freq))) #Creating the base array for the PSD
            data_comp[:,2] = 0 #Resetting the simulated OLC

            colour_change = colour(-slope,len(cont[:,1]))
            #colour_change = colour_change/np.std(colour_change)
            
            change = colour_change #Finding the magnitude and direction of change
            print np.nanmean(change), np.std(change)
            cont[:,1] = (cont[:,1] * (change))/abs(np.mean(change)) #/np.sqrt(2.) #Implementing change
            #print change
            #print cont[:,1]
            h = 0
            b = np.nanmean(cont[:,1]) #Finding mean value of new Kelly 2009
            flux1 = cont[:,1] #Updating global variable for Kelly 2009
            param = [b,tau,sigma_tot] #See above

            cont_real = cont
            cont_real[:,1] = ifft(cont[:,1]).real
            #cont_real[:,1] = cont_real[:,1] #/np.std(cont_real[:,1])

            slope, P_v = PSD(freq,transfer_array,cont_real) #Finding PSD slope

            data_comp[:,2] = model_data(cont_real,data_comp,transfer_array) #The data after running through the transfer function

            chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2) #Finding residuals squared summed

            slopechange = (slope - slopeaim)**2 - (slopeolder1 - slopeaim)**2 #Is the new PSD slope a better fit.
            

            if chi2 <= chi1 \
               and slopechange < 0.:
                '''The first instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 1, slope,j

            elif chi2 <= chi1 \
               and abs(slope) < abs(slopeaim - slopeallow):
                '''The first instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 1.1, slope,j

            elif chi2 <= chi1 \
                 and slopeaim - slopeallow < slope < slopeaim + slopeallow:
                '''The second instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 2, slope,j

            elif random.random < 0.05 \
                 and abs(slope) < abs(slopeaim + slopeallow): # and (dd2_ddt1 - dd2_ddt) < 0 and (dd3_ddt1 - dd3_ddt) < 0: # \
                 #and random.random() < 0.1:
                '''The third instance of acceptance, an attempt to encourage a smoother curve.'''
                chi1 = chi2
                slopeolder1 = slope
                print 3, slope,j

            else:
                '''Rejection'''
                cont = cont_save

            #if chi2 < chi1:
            #    print chi2, chi1

            #model3 = model2
            #print cont_real
            gc.collect()
    
    P_show = log(P_v)
    #np.savetxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1',cont) #Updating the files
    #np.savetxt('CONTINUUM/NGC3783-data_comp-slope-Kelly-2_5-1',data_comp)
    plt.figure()
    plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    plt.scatter(cont_real[:,0],cont_real[:,1],color='g')
    #plt.ylim([4e-15,1.2e-14])
    plt.ylim([4e-30,1.2e5])
    plt.yscale('log')
    plt.show(block=False)

print cont[:,1]
