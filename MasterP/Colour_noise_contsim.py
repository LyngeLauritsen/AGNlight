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


#ifft(s).real
#y / std(y)

'''Usefull Physical Constants for the program'''
h_c = 6.626*10**(-34)   #Planck's Constant
c_c = 3.*10**8          #Speed of light
k_B = 1.38*10**(-23)    #Boltzmann Constant
conv = 10**(-9)         #Convergence from nm to m

''' Usefull parameters for observational bands'''
K_cen = 2190*10**(-9)   #Center of K band in m
H_cen = 1630*10**(-9)   #Center of H band in m
J_cen = 1220*10**(-9)   #Center of J band in m
g_cen = 474.7*10**(-9)   #Center of g band in m
r_cen = 621.4*10**(-9)   #Center of r band in m
i_cen = 762.8*10**(-9)   #Center of i band in m
z_cen = 906.8*10**(-9)   #Center of z band in m

K_width = 390*10**(-9)  #Width of K band in m
H_width = 307*10**(-9)  #Width of H band in m
J_width = 213*10**(-9)  #Width of J band in m
g_width = 99.6*10**(-9)  #Width of g band in m
r_width = 95.7*10**(-9)  #Width of r band in m
i_width = 105.6*10**(-9)  #Width of i band in m
z_width = 122.7*10**(-9)  #Width of z band in m

runs = 100 #3000 #Number of runs (arbitrary due to the gradual updates, would probably take close to a month to run)
runs1 = 100 #1000 #1000
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

'''Changeable constants'''
T = 3000 # K
wl_slope_T = 0.5
w_T = 0.5
wl_slope_power = 0.5
wl_intercept_power = 0.5
w_power = 0.5
w_intercept_power = 0.5



# There are signal processing package in scip y that defines filters, prob. can find lognormal in there.
def lognorm(x,mu,sigma):
    '''Defines the transfer function'''
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)

def planck(wl,temperature):
    front = 2*h_c*c_c**2/wl**5
    exp = h_c*c_c/(wl*k_B*temperature)
    return front*1/(np.exp(exp)-1)

def colour(alpha,length):
    return cn.powerlaw_psd_gaussian(alpha,length)

def lag_equation(wl,wl_slope,wl_intercept):
    return wl_intercept*wl**wl_slope

def w_equation(w,w_slope,w_intercept):
    return w_intercept*w**w_slope


x_list = np.linspace(0.01,1500,1500) #X-axis of transfer function



log_norm_transfer = []
for i in range(len(x_list)):
    '''Creates the list used to identify the degree of transfered light'''
    log_norm_transfer.append(lognorm(x_list[i],mu,sigma))
#Suggest:
# log = lognorm(np.arange(len(x_list)))

print np.nansum(log_norm_transfer)*1500/1500.

data_K = np.loadtxt('NOVEMBER/NOV-NGC3783-K') #'NOVEMBER/Kelly-NGC3783K')
error_K = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt') #Load data

data_H = np.loadtxt('NOVEMBER/NOV-NGC3783-H') #'NOVEMBER/Kelly-NGC3783K')
error_H = np.loadtxt('NOVEMBER/NGC3783_NOISE_H.txt') #Load data

data_J = np.loadtxt('NOVEMBER/NOV-NGC3783-J') #'NOVEMBER/Kelly-NGC3783K')
error_J = np.loadtxt('NOVEMBER/NGC3783_NOISE_J.txt') #Load data

data_g = np.loadtxt('NOVEMBER/NOV-NGC3783-g') #'NOVEMBER/Kelly-NGC3783K')
error_g = np.loadtxt('NOVEMBER/NGC3783_NOISE_g.txt') #Load data

data_r = np.loadtxt('NOVEMBER/NOV-NGC3783-r') #'NOVEMBER/Kelly-NGC3783K')
error_r = np.loadtxt('NOVEMBER/NGC3783_NOISE_r.txt') #Load data

data_i = np.loadtxt('NOVEMBER/NOV-NGC3783-i') #'NOVEMBER/Kelly-NGC3783K')
error_i = np.loadtxt('NOVEMBER/NGC3783_NOISE_i.txt') #Load data

data_z = np.loadtxt('NOVEMBER/NOV-NGC3783-z') #'NOVEMBER/Kelly-NGC3783K')
error_z = np.loadtxt('NOVEMBER/NGC3783_NOISE_z.txt') #Load data

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

'''


def linear(a_c,b_c,lambda):
    return a*lambda + b_c
'''


dd_ddt_allowed = np.nanmax(dd_ddt_data)*dd_change #Identifies the maximal double derivative allowed in the continuum light curve

'''cont denotes the model for the continuum light curve, whereas data_comp denotes the array giving indication of the
observed light curve and how the continuum light curve fitted through the transfer function compares'''
cont_days = np.arange(min(data[:,0])-days_before,max(data[:,0]+days_after),step) #Defines the spacing of the light curve
#cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]+100),step) #Defines the spacing of the light curve
cont = np.zeros((len(cont_days),3))

colour_start = colour(-slope,len(cont[:,1]))
colour_start = colour_start/np.std(colour_start)

cont[:,0] = cont_days
cont[:,1] = colour_start
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1')
day = cont_days[0]

data_comp_K = np.zeros((len(data_K[:,1]),3))
data_comp_K[:,0] = data_K[:,0]
data_comp_K[:,1] = data_K[:,1]
data_comp_K[:,2] = 0

data_comp_H = np.zeros((len(data_H[:,1]),3))
data_comp_H[:,0] = data_H[:,0]
data_comp_H[:,1] = data_H[:,1]
data_comp_H[:,2] = 0

data_comp_J = np.zeros((len(data_J[:,1]),3))
data_comp_J[:,0] = data_J[:,0]
data_comp_J[:,1] = data_J[:,1]
data_comp_J[:,2] = 0

data_comp_g = np.zeros((len(data_g[:,1]),3))
data_comp_g[:,0] = data_g[:,0]
data_comp_g[:,1] = data_g[:,1]
data_comp_g[:,2] = 0

data_comp_r = np.zeros((len(data_r[:,1]),3))
data_comp_r[:,0] = data_r[:,0]
data_comp_r[:,1] = data_r[:,1]
data_comp_r[:,2] = 0

data_comp_i = np.zeros((len(data_i[:,1]),3))
data_comp_i[:,0] = data_i[:,0]
data_comp_i[:,1] = data_i[:,1]
data_comp_i[:,2] = 0

data_comp_z = np.zeros((len(data_z[:,1]),3))
data_comp_z[:,0] = data_z[:,0]
data_comp_z[:,1] = data_z[:,1]
data_comp_z[:,2] = 0
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
max_jump = 0.000001 #Maximal allowed jump relative to the value at point

plt.figure()
plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
plt.scatter(cont[:,0],ifft(cont[:,1])*np.nanmean(data[:,1]).real,color='yellow')
#plt.scatter(cont[:,0],cont[:,1],color='g')
plt.ylim([1.2e-25,1e1])
#plt.ylim([4e-15,1.2e-14])
plt.yscale('log')
plt.show(block=False)

chi1 = 1e1000 #np.nansum((data_comp[:,1] - data_comp[:,2])**2)

slopeolder1 = 0 #Temporary constant

F_N_v = np.zeros((len(freq)))
P_v = 0
print F_N_v

for i in range(runs):
    '''The MCMC'''
    for i1 in range(runs1):
        '''Doing 2 loops to allow regular progress reports'''
        #print i, i1/float(runs1), slopeolder1
        '''Changing the parameters'''
        T += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) # K
        wl_slope_T += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)
        w_T += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)
        wl_slope_power += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)
        wl_intercept_power += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)
        w_power += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)
        w_intercept_power += (-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1)

        '''Looping over the course of the CLC'''
        F_N_v = np.zeros((len(freq))) #Creating the base array for the PSD
        #print np.shape(colour(-slope,len(cont[:,1])))
        colour_change = colour(-slope,len(cont[:,1]))
        change = colour_change*np.nansum(cont[:,1])*max_jump #Finding the magnitude and direction of change
        cont[:,1] += change #Implementing change
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

        data_comp_K[:,2] = model_data(cont_real,data_comp_K,transfer_array) #The data after running through the transfer function
        data_comp_H[:,2] = model_data(cont_real,data_comp_H,transfer_array)
        data_comp_J[:,2] = model_data(cont_real,data_comp_J,transfer_array)
        data_comp_g[:,2] = model_data(cont_real,data_comp_g,transfer_array)
        data_comp_r[:,2] = model_data(cont_real,data_comp_r,transfer_array)
        data_comp_i[:,2] = model_data(cont_real,data_comp_i,transfer_array)
        data_comp_z[:,2] = model_data(cont_real,data_comp_z,transfer_array)

        chi2_K = np.nansum((data_comp_K[:,1] - data_comp_K[:,2])**2) #Finding residuals squared summed
        chi2_H = np.nansum((data_comp_H[:,1] - data_comp_H[:,2])**2)
        chi2_J = np.nansum((data_comp_J[:,1] - data_comp_J[:,2])**2)
        chi2_g = np.nansum((data_comp_g[:,1] - data_comp_g[:,2])**2)
        chi2_r = np.nansum((data_comp_r[:,1] - data_comp_r[:,2])**2)
        chi2_i = np.nansum((data_comp_i[:,1] - data_comp_i[:,2])**2)
        chi2_z = np.nansum((data_comp_z[:,1] - data_comp_z[:,2])**2)

        chi2 = chi2_K + chi2_H + chi2_J + chi2_g + chi2_r + chi2_i + chi2_z

        print chi2 - chi1
        slopechange = (slope - slopeaim)**2 - (slopeolder1 - slopeaim)**2 #Is the new PSD slope a better fit.

        MCMC = random.random()
        #if (chi2 <= chi1 \
        #    and slopechange < 0.):
        #    '''The first instance of acceptance'''
        #    chi1 = chi2
        #    slopeolder1 = slope
            #print 1, slope

        if chi2 <= chi1 \
           and abs(slope) < abs(slopeaim - slopeallow): #abs(slopeolder1) <
            '''The first instance of acceptance'''
            chi1 = chi2
            slopeolder1 = slope
            #print 1.1, slope

        elif chi2 <= chi1 \
             and slopeaim - slopeallow < slope < slopeaim + slopeallow:
            '''The second instance of acceptance'''
            chi1 = chi2
            slopeolder1 = slope
            #print 2, slope

        elif 0.1 > MCMC \
             and abs(slope) < abs(slopeaim + slopeallow): # and (dd2_ddt1 - dd2_ddt) < 0 and (dd3_ddt1 - dd3_ddt) < 0: # \
             #and random.random() < 0.1:
            '''The third instance of acceptance, an attempt to encourage a smoother curve.'''
            chi1 = chi2
            slopeolder1 = slope
            #print MCMC
            #print 3, slope

        else:
            '''Rejection'''
            cont[:,1] -= change

        #if chi2 < chi1:
        #    print chi2, chi1

        #model3 = model2
        gc.collect()
    P_show = log(P_v)
    np.savetxt('CONTINUUM/NGC3783-continuum-slope-colour-2_5',cont_real) #Updating the files
    np.savetxt('CONTINUUM/NGC3783-data_comp-slope-colour-2_5',data_comp)
    plt.figure()
    plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    plt.scatter(cont[:,0],ifft(cont[:,1]).real,color='yellow')
    plt.scatter(cont[:,0],cont[:,1],color='g')
    plt.ylim([1.2e-25,1e1])
    #plt.ylim([4e-15,1.2e-14])
    plt.yscale('log')
    plt.show(block=False)
    #plt.figure()
    #plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    #plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    #plt.scatter(cont[:,0],cont[:,1],color='g')
    #plt.ylim([4e-15,1.2e-14])
    #plt.show(block=False)
    print slopeolder1
    #plt.figure()
    #plt.scatter(freq,P_show,color='b')
    #plt.ylim([1e-21,1e-31])
    #plt.xlim([1e-8,3e-7])
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.show(block=False)
    gc.collect()

#print cont[:,1]
