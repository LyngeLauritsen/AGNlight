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
#import colorednoise as cn
from numpy import concatenate, real, std, abs, min
from numpy.fft import ifft, fftfreq, fft, rfft, irfft
from numpy.random import normal
import time

def powerlaw_psd_gaussian(exponent, samples, fmin=0):
    """Gaussian (1/f)**beta noise.
    """

    # frequencies (we asume a sample rate of one)
    f = fftfreq(samples)

    # scaling factor for all frequencies
    ## though the fft for real signals is symmetric,
    ## the array with the results is not - take neg. half!
    s_scale = abs(concatenate([f[f<0], [f[-1]]]))
    ## low frequency cutoff?!?
    if fmin:
        ix = sum(s_scale>fmin)
        if ix < len(f):
            s_scale[ix:] = s_scale[ix]
    s_scale = s_scale**(-exponent/2.)

    # scale random power + phase
    sr = s_scale * normal(size=len(s_scale))
    si = s_scale * normal(size=len(s_scale))
    if not (samples % 2): si[0] = si[0].real

    s = sr + 1J * si
    # this is complicated... because for odd sample numbers,ss
    ## there is one less positive freq than for even sample numbers
    s = concatenate([s[1-(samples % 2):][::-1], s[:-1].conj()])

    # time series
    y = ifft(s).real
    y*= 10**np.random.normal(-13.7,1.5,1)[0]
    y += abs(min(y))

    #y = rfft(y)
    #y = irfft(y)

    #y1 = y / std(y)
    #print y1

    #y1 = y1

    #y2 = ifft(y1)

    ret = y

    #print 'HEY', ret
    #print 'HEY HO'

    return ret

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
runs1 = 2500 #1000 #1000
days_before = 700. #How far to extend the lightcurve back (days)
days_after = 10. #How far to extend the lightcurve forward (days)
step = 5. #Timestep (days)
slopeaim = -2.7 #Slope that is targeted
slopeallow = -0.5
slope = 0
#mean = 140. #Transfer function
#sigma = 0.8
#mu = np.log(mean) - 1/2.*sigma
kelly_repeat = 20
dd_change = 1.1

'''Changeable constants'''
T = 3000. # K
lag_thermal = 3.
width_thermal = 1.
lag_slope_power = 0.
lag_intercept_power = 30. #At 0 m
width_slope_power = 0.
width_intercept_power = 30. #At 0 m
A_T = 0.5
N_s_power = np.array([.5,.5,.5,.5,.5,.5,.5])
print N_s_power
scale = 12.5

'''Fitting the constants'''
mu_power_K = 3.
mu_power_H = 3.
mu_power_J = 3.
mu_power_g = 3.
mu_power_r = 3.
mu_power_i = 3.
mu_power_z = 3.

sigma_power_K = 1.
sigma_power_H = 1.
sigma_power_J = 1.
sigma_power_g = 1.
sigma_power_r = 1.
sigma_power_i = 1.
sigma_power_z = 1.

def create_lognorm(x,mu,sigma):
    '''Defines the transfer function'''
    delay_days = x
    sigma_array = np.zeros((len(delay_days)))
    sigma_array.fill(float(sigma))

    mu_array = np.zeros((len(delay_days)))
    mu_array.fill(float(mu))
    #print (2*sigma_array**2)
    exp_term = np.zeros((np.shape(delay_days)))
    exp_term = -((np.log(x)-mu_array)**2/(2*sigma_array**2))
    front = 1/(x*sigma_array*np.sqrt(2*np.pi))
    #print front*np.exp(exp_term)
    #print np.shape(front*np.exp(exp_term))
    #print front
    result = front*np.exp(exp_term)
    #print np.sum(result)
    return result #/max(result)

def planck(wl,temperature):
    front = 2*h_c*c_c**2/wl**5
    exp = h_c*c_c/(wl*k_B*temperature)
    return front*1/(np.exp(exp)-1)

def colour(alpha,length):
    return powerlaw_psd_gaussian(alpha,length)

def lag_equation(wl,lag_slope,lag_intercept):
    return lag_intercept*wl**lag_slope

def width_equation(w,width_slope,width_intercept):
    return width_intercept*w**width_slope

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

cont = np.loadtxt('NOVEMBER/NGC3783_CONT_TRY')
#plt.figure()
#plt.scatter(data[:,0],data[:,1])
#plt.ylim([6e-15,1.1e-14])

#plt.show()

dd_ddt_data = []
#for j in range(len(data[:,0])-2):
#    '''Purely a test to see the largest double derivative of the observed light curve'''
#    dd_ddt_data.append(abs(((data[j+1,1] - data[j,1])/(data[j+1,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+1,0] - data[j-1,0])/2.))

#
#    dd_ddt_data.append(abs(((data[j+4,1] - data[j,1])/(data[j+4,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+4


 #Identifies the maximal double derivative allowed in the continuum light curve

'''cont denotes the model for the continuum light curve, whereas data_comp denotes the array giving indication of the
observed light curve and how the continuum light curve fitted through the transfer function compares'''
cont_days = cont[:,0] #np.arange(min(data_K[:,0])-days_before,max(data_K[:,0])+days_after,step) #Defines the spacing of the light curve
cont_fft = rfft(cont[:,1])
cont_ifft = irfft(cont_fft)

print cont[:,1][0],cont[:,1][1],cont[:,1][2],cont[:,1][3],cont[:,1][4]
print cont_ifft[0],cont_ifft[1],cont_ifft[2],cont_ifft[3],cont_ifft[4]

time.sleep(5)

T = np.array([T]*len(cont_days)) #List of temperatures

'''Arrays defining the band filters'''
K_cen = np.array([K_cen]*len(cont_days))
H_cen = np.array([H_cen]*len(cont_days))
J_cen = np.array([J_cen]*len(cont_days))
g_cen = np.array([g_cen]*len(cont_days))
r_cen = np.array([r_cen]*len(cont_days))
i_cen = np.array([i_cen]*len(cont_days))
z_cen = np.array([z_cen]*len(cont_days))

K_width = np.array([K_width]*len(cont_days))
H_width = np.array([H_width]*len(cont_days))
J_width = np.array([J_width]*len(cont_days))
g_width = np.array([g_width]*len(cont_days))
r_width = np.array([r_width]*len(cont_days))
i_width = np.array([i_width]*len(cont_days))
z_width = np.array([z_width]*len(cont_days))

#cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]+100),step) #Defines the spacing of the light curve
#cont = np.empty((len(cont_days),3),dtype=complex)

colour_start = colour(np.random.normal(2.7,0.15,1)[0],len(cont[:,1])) #colour(-slope,len(cont[:,1]))
#print colour_start[0]
#colour_start = colour_start*1e-17

#cont[:,0] = np.real(cont_days)
#!cont[:,1] = colour_start #colour(2.7,len(cont[:,1])) #np.random.randint(-1e5,1e5,len(cont_days))*1.e-22 #colour(-slope,len(cont[:,1]))*1e-17 #colour_start
#print 'Cont 1 =', cont[0,1]
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1')
#day = cont_days[0]

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
data_comp_g[:,1] = data_g[:,2]
data_comp_g[:,2] = 0

data_comp_r = np.zeros((len(data_r[:,1]),3))
data_comp_r[:,0] = data_r[:,0]
data_comp_r[:,1] = data_r[:,2]
data_comp_r[:,2] = 0

data_comp_i = np.zeros((len(data_i[:,1]),3))
data_comp_i[:,0] = data_i[:,0]
data_comp_i[:,1] = data_i[:,2]
data_comp_i[:,2] = 0

data_comp_z = np.zeros((len(data_z[:,1]),3))
data_comp_z[:,0] = data_z[:,0]
data_comp_z[:,1] = data_z[:,2]
data_comp_z[:,2] = 0
#data_comp = np.loadtxt('CONTINUUM/NGC3783-data_comp-slope-Kelly-2_5-1')
#print data[:,2]

def transfer(cont_days,data_comp,log_norm_transfer):
    '''Transfer function array'''
    cont_days = np.real(cont_days)
    transfer = np.zeros((len(cont_days),len(data_comp[:,0]))) #Array with the transferred light in any combination of the continuum and the observed.
    #print data_comp[:,0]
    #print cont[:,0]
    for k in range(len(cont_days)-1):
        '''Creating the array for the amount of transmitted light'''
        mask = data_comp[:,0] > cont_days[k]
        z = data_comp[:,0][~mask]
        z.fill(nan)
        #print z
        #print np.shape(log_norm_transfer), np.sum(log_norm_transfer)
        #print np.shape(cont_days)
        #plt.figure()
        #plt.plot(cont_days,log_norm_transfer[:])
        #plt.show()
        #print np.concatenate([z,data_comp[:,0][mask]])
        #print cont[k,0]
        #print abs(cont[k,0] - np.concatenate(([z,data_comp[:,0][mask]]))).astype('int')
        #print np.shape(transfer), np.shape(z), np.shape(data_comp[:,0][mask])
        #print np.shape(np.concatenate(([z,data_comp[:,0][mask]]))), np.shape(log_norm_transfer)
        delay = abs(int(np.real(cont[k,0])) - np.concatenate(([z,data_comp[:,0][mask]])))
        #print delay
        delay[np.isnan(delay)] = 0
        #print delay
        transfer[k,:] = log_norm_transfer[delay.astype('int')]
    return transfer.transpose()

def model_data(cont,data_comp,transfer_array):
    '''The simulated observed light curve as a result of the continuum light curve and the transfer function'''

    '''The creation of the transfer function arrays'''
    #transfer_array_thermal = transfer(cont[:,0],cont,data_comp,x,mu_thermal,sigma_thermal,normalize_thermal)
    #transfer_array_power = transfer(cont[:,0],cont,data_comp,x,mu_power,sigma_power,normalize_power)

    cont_flux = np.array([cont[:,1],]*len(data_comp[:,0]))

    '''The simulated observed light curves for the two transfer functions'''
    sim_OLC = np.sum(step*cont_flux*transfer_array,axis=1)

    return sim_OLC

#transfer_array = transfer(cont_days,data_comp,log_norm_transfer)

#flux1 = cont[:,1] #Series of items for the Kelly 2009
#time1 = cont[:,0]
#sigma1 = cont[:,2] + np.nanmean(error[:,1])
#ime = np.insert(time1,0,0)
#sigma = np.insert(sigma1,0,np.mean(sigma1))
#flux = np.insert(flux1,0,np.mean(flux1))

#print shape(flux1)

freq = []

for i in range(int(len(cont_days)/2.)):
    '''The frequency determination for the PSD'''
    freq.append((i+1)/(60.*60.*24*float(cont_days[len(cont_days)-1] - cont_days[0])))

def PSD(freq,cont):
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

#def model_data(cont,data_comp,transfer_array):
#    '''The simulated observed light curve as a result of the continuum light curve and the transfer function'''
#
#    cont_flux = array([cont[:,1],]*len(data_comp[:,0]))
#
#    return sum(step*cont_flux*transfer_array,axis=1)

#model3 = cont[:,1]
tau = 800. #The Kelly relaxation time
#sigma_tot = np.nanmean(error[:,1]) #Kelly sigma
chi1 = 1e100 #To ensure that the first change is always accepted to get it started
max_jump = 0.000001 #Maximal allowed jump relative to the value at point

#plt.figure()
#plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
#plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
#plt.scatter(cont[:,0],ifft(cont[:,1])*np.nanmean(data[:,1]).real,color='yellow')
#plt.scatter(cont[:,0],cont[:,1],color='g')
#plt.ylim([1.2e-25,1e1])
#plt.ylim([4e-15,1.2e-14])
#plt.yscale('log')
#plt.show(block=False)

chi1 = 1e1000 #np.nansum((data_comp[:,1] - data_comp[:,2])**2)

'''Determining the quality of the various fits'''
chi1_K = np.nansum((data_comp_K[:,1] - data_comp_K[:,2])**2) #Finding residuals squared summed
chi1_H = np.nansum((data_comp_H[:,1] - data_comp_H[:,2])**2)
chi1_J = np.nansum((data_comp_J[:,1] - data_comp_J[:,2])**2)
chi1_g = np.nansum((data_comp_g[:,1] - data_comp_g[:,2])**2)
chi1_r = np.nansum((data_comp_r[:,1] - data_comp_r[:,2])**2)
chi1_i = np.nansum((data_comp_i[:,1] - data_comp_i[:,2])**2)
chi1_z = np.nansum((data_comp_z[:,1] - data_comp_z[:,2])**2)


slopeolder1 = 0 #Temporary constant

F_N_v = np.zeros((len(freq)))
P_v = 0
#print F_N_v

cont_real = np.real(cont)
cont_real_save = cont
cont_save = cont

data_comp_K_save = data_comp_K
data_comp_H_save = data_comp_H
data_comp_J_save = data_comp_J
data_comp_g_save = data_comp_g
data_comp_r_save = data_comp_r
data_comp_i_save = data_comp_i
data_comp_z_save = data_comp_z
scale_save = scale

x_list = np.linspace(0.01,len(cont_days)*step,len(cont_days)*step)

jump = 1.*10**(-3)
slope_jump = 0.005
jump_change = 0
slope_jump_change = 0

try1 = 0

T_direction = (-1)**np.random.randint(2,size=1) #np.array([(-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) for i in range(len(cont_days))]) # K
#print np.random.randint(2)
lag_thermal_direction = (-1)**np.random.randint(2)
'''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
width_thermal_direction = (-1)**np.random.randint(2)#*np.random.randint(5,size=1)
lag_slope_power_direction = (-1)**np.random.randint(2)
lag_intercept_power_direction = (-1)**np.random.randint(2)#*np.random.randint(5,size=1)
width_slope_power_direction = (-1)**np.random.randint(2)
width_intercept_power_direction = (-1)**np.random.randint(2)#*np.random.randint(5,size=1)
A_T_direction = (-1)**np.random.randint(2)
#print N_s_power, np.array([(-1)**np.random.randint(2)*random.random()*0.0001 for i in range(7)]) #.transpose()
N_s_power_direction = (-1)**np.random.randint(2,size=7) #.transpose()
scale_direction = (-1)**np.random.randint(2)

mu_power_K_save = np.copy(mu_power_K)
mu_power_H_save = np.copy(mu_power_H)
mu_power_J_save = np.copy(mu_power_J)
mu_power_g_save = np.copy(mu_power_g)
mu_power_r_save = np.copy(mu_power_r)
mu_power_i_save = np.copy(mu_power_i)
mu_power_z_save = np.copy(mu_power_z)

sigma_power_K_save = np.copy(sigma_power_K)
sigma_power_H_save = np.copy(sigma_power_H)
sigma_power_J_save = np.copy(sigma_power_J)
sigma_power_g_save = np.copy(sigma_power_g)
sigma_power_r_save = np.copy(sigma_power_r)
sigma_power_i_save = np.copy(sigma_power_i)
sigma_power_z_save = np.copy(sigma_power_z)

cont_fft_save = np.copy(cont_fft)

mu_power_K_direction = (-1)**np.random.randint(2,size=1)
mu_power_H_direction = (-1)**np.random.randint(2,size=1)
mu_power_J_direction = (-1)**np.random.randint(2,size=1)
mu_power_g_direction = (-1)**np.random.randint(2,size=1)
mu_power_r_direction = (-1)**np.random.randint(2,size=1)
mu_power_i_direction = (-1)**np.random.randint(2,size=1)
mu_power_z_direction = (-1)**np.random.randint(2,size=1)

sigma_power_K_direction = (-1)**np.random.randint(2,size=1)
sigma_power_H_direction = (-1)**np.random.randint(2,size=1)
sigma_power_J_direction = (-1)**np.random.randint(2,size=1)
sigma_power_g_direction = (-1)**np.random.randint(2,size=1)
sigma_power_r_direction = (-1)**np.random.randint(2,size=1)
sigma_power_i_direction = (-1)**np.random.randint(2,size=1)
sigma_power_z_direction = (-1)**np.random.randint(2,size=1)

TCK = []
SK = []
MK = []

#print np.real(ifft(cont[:,1]).real), PSD(freq,np.real((ifft(cont[:,1]).real - np.min(ifft(cont[:,1]).real)))*10**(-scale))

for i in range(runs):
    '''The MCMC'''
    start_time = time.time()
    for i1 in range(runs1):
        '''Doing 2 loops to allow regular progress reports'''

        jump_change += 1
        slope_jump_change += 1
        #print jump
        if jump_change > 50:
            #jump *= 0.7
            jump = 10**np.random.normal(-4,2,1)[0]
            jump_change = 0
            print 'jump =', jump
        if slope_jump_change > 300000:
            slope_jump *= 0.95

        '''Saving original arrays'''
        #print 'cont', cont[0,1]
        #print 'cont save', cont_save[0,1]
        cont_real_save = np.copy(cont_real)
        cont_save = np.copy(cont)
        #print cont_real[:,1]
        #time.sleep(10)

        data_comp_K_save = np.copy(data_comp_K)
        data_comp_H_save = np.copy(data_comp_H)
        data_comp_J_save = np.copy(data_comp_J)
        data_comp_g_save = np.copy(data_comp_g)
        data_comp_r_save = np.copy(data_comp_r)
        data_comp_i_save = np.copy(data_comp_i)
        data_comp_z_save = np.copy(data_comp_z)

        T_save = np.copy(T) # K
        lag_thermal_save = np.copy(lag_thermal)
        width_thermal_save = np.copy(width_thermal)
        lag_slope_power_save = np.copy(lag_slope_power)
        lag_intercept_power_save = np.copy(lag_intercept_power)
        width_slope_power_save = np.copy(width_slope_power)
        width_intercept_power_save = np.copy(width_intercept_power)
        A_T_save = np.copy(A_T)
        N_s_power_save = np.copy(N_s_power)
        scale_save = np.copy(scale)

        cont_fft_save = np.copy(cont_fft)

        mu_power_K_save = np.copy(mu_power_K)
        mu_power_H_save = np.copy(mu_power_H)
        mu_power_J_save = np.copy(mu_power_J)
        mu_power_g_save = np.copy(mu_power_g)
        mu_power_r_save = np.copy(mu_power_r)
        mu_power_i_save = np.copy(mu_power_i)
        mu_power_z_save = np.copy(mu_power_z)

        sigma_power_K_save = np.copy(sigma_power_K)
        sigma_power_H_save = np.copy(sigma_power_H)
        sigma_power_J_save = np.copy(sigma_power_J)
        sigma_power_g_save = np.copy(sigma_power_g)
        sigma_power_r_save = np.copy(sigma_power_r)
        sigma_power_i_save = np.copy(sigma_power_i)
        sigma_power_z_save = np.copy(sigma_power_z)


        if try1 == 0:
            #print i, i1/float(runs1), slopeolder1
            '''Changing the parameters'''
            change_size = random.randint(2,40)
            change = colour(np.random.normal(2.8,0.25,1)[0],change_size) #len(cont[:,1]))) rfft
            origin = np.copy(irfft(cont_fft))
            place = random.randint(0,len(cont_real[:,1])-change_size)
            #print np.shape(origin),np.shape(change), place
            origin[place:place+change_size] = change
            origin2 = rfft(origin)
            weight = random.random()*0.005 #0.05
            cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(origin2))
            #print try1
        if try1 == 1.1:
            '''Thermal component'''
            T_direction = (-1)**np.random.randint(2,size=1) #np.array([(-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) for i in range(len(cont_days))]) # K
            lag_thermal_direction = (-1)**np.random.randint(2)
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal_direction = (-1)**np.random.randint(2)#*np.random.randint(5,size=1)

            T += T_direction*np.random.random()*5
            lag_thermal += lag_thermal_direction*random.random()*jump*5
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal += width_thermal_direction*random.random()*jump*5#*np.random.randint(5,size=1)

        elif try1 == 1.2:
            '''Power component'''
            mu_power_K_direction = (-1)**np.random.randint(2)
            mu_power_H_direction = (-1)**np.random.randint(2)
            mu_power_J_direction = (-1)**np.random.randint(2)
            mu_power_g_direction = (-1)**np.random.randint(2)
            mu_power_r_direction = (-1)**np.random.randint(2)
            mu_power_i_direction = (-1)**np.random.randint(2)
            mu_power_z_direction = (-1)**np.random.randint(2)
            
            sigma_power_K_direction = (-1)**np.random.randint(2)
            sigma_power_H_direction = (-1)**np.random.randint(2)
            sigma_power_J_direction = (-1)**np.random.randint(2)
            sigma_power_g_direction = (-1)**np.random.randint(2)
            sigma_power_r_direction = (-1)**np.random.randint(2)
            sigma_power_i_direction = (-1)**np.random.randint(2)
            sigma_power_z_direction = (-1)**np.random.randint(2)

            mu_power_K += mu_power_K_direction*random.random()*jump
            mu_power_H += mu_power_H_direction*random.random()*jump
            mu_power_J += mu_power_J_direction*random.random()*jump
            mu_power_g += mu_power_g_direction*random.random()*jump
            mu_power_r += mu_power_r_direction*random.random()*jump
            mu_power_i += mu_power_i_direction*random.random()*jump
            mu_power_z += mu_power_z_direction*random.random()*jump

            sigma_power_K += sigma_power_K_direction*random.random()*jump
            sigma_power_H += sigma_power_H_direction*random.random()*jump
            sigma_power_J += sigma_power_J_direction*random.random()*jump
            sigma_power_g += sigma_power_g_direction*random.random()*jump
            sigma_power_r += sigma_power_r_direction*random.random()*jump
            sigma_power_i += sigma_power_i_direction*random.random()*jump
            sigma_power_z += sigma_power_z_direction*random.random()*jump

        elif try1 == 1.3:
            N_s_power_direction = (-1)**np.random.randint(2,size=7)
            A_T_direction = (-1)**np.random.randint(2)
            scale_direction = (-1)**np.random.randint(2)

            N_s_power += N_s_power_direction*np.random.rand(7)*0.5*jump #.transpose()
            A_T += A_T_direction*random.random()*jump
            scale += scale_direction*random.random()/1000.

        elif try1 == 2.1:
            T += (-1)**np.random.randint(2,size=1)*np.random.random()*5 #np.array([(-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) for i in range(len(cont_days))]) # K
            #print np.random.randint(2)
            lag_thermal += lag_thermal_direction*random.random()*jump*5
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal += width_thermal_direction*random.random()*jump*5#*np.random.randint(5,size=1)

        elif try1 == 2.2:
            mu_power_K += mu_power_K_direction*random.random()*jump
            mu_power_H += mu_power_H_direction*random.random()*jump
            mu_power_J += mu_power_J_direction*random.random()*jump
            mu_power_g += mu_power_g_direction*random.random()*jump
            mu_power_r += mu_power_r_direction*random.random()*jump
            mu_power_i += mu_power_i_direction*random.random()*jump
            mu_power_z += mu_power_z_direction*random.random()*jump

            sigma_power_K += sigma_power_K_direction*random.random()*jump
            sigma_power_H += sigma_power_H_direction*random.random()*jump
            sigma_power_J += sigma_power_J_direction*random.random()*jump
            sigma_power_g += sigma_power_g_direction*random.random()*jump
            sigma_power_r += sigma_power_r_direction*random.random()*jump
            sigma_power_i += sigma_power_i_direction*random.random()*jump
            sigma_power_z += sigma_power_z_direction*random.random()*jump
            
        elif try1 == 2.3:
            A_T += A_T_direction*random.random()*jump
            N_s_power += N_s_power_direction*np.random.rand(7)*jump #.transpose()
            scale += scale_direction*random.random()/1000.

        N_s_power = abs(N_s_power)
        mu_power_K = abs(mu_power_K)
        mu_power_H = abs(mu_power_H)
        mu_power_J = abs(mu_power_J)
        mu_power_g = abs(mu_power_g)
        mu_power_r = abs(mu_power_r)
        mu_power_i = abs(mu_power_i)
        mu_power_z = abs(mu_power_z)

        sigma_power_K = abs(sigma_power_K)
        sigma_power_H = abs(sigma_power_H)
        sigma_power_J = abs(sigma_power_J)
        sigma_power_g = abs(sigma_power_g)
        sigma_power_r = abs(sigma_power_r)
        sigma_power_i = abs(sigma_power_i)
        sigma_power_z = abs(sigma_power_z)

        #print K_width

        '''Looping over the course of the CLC'''
        F_N_v = np.zeros((len(freq))) #Creating the base array for the PSD
        #print np.shape(colour(-slope,len(cont[:,1])))
        #colour_change = colour(-slope,len(cont[:,1]))
        #change = colour_change
        #print change[0]
        #print cont[:,1][0]
        '''LOOK AT THE FFT FOR POSSIBLE PHASE SHIFT'''


        #print 'cont', cont[0,1]
        #print 'cont save', cont_save[0,1]
        #print cont[:,1]
        #print change
        cont_real[:,1] = irfft(cont_fft)#*10**(-scale) #Implementing change
        #print change[1]
        #print 'real', cont_real[0,1]
        #h = 0
        #b = np.nanmean(cont[:,1]) #Finding mean value of new Kelly 2009
        #flux1 = cont[:,1] #Updating global variable for Kelly 2009
        #param = [b,tau,sigma_tot] #See above

        #cont_real = cont
        #cont_real[:,1] = ifft(cont[:,1]).real
        #cont_real[:,1] = cont_real[:,1] #/np.std(cont_real[:,1])

        slope, P_v = PSD(freq,cont_real) #Finding PSD slope
        #print slope

        '''Power transfer functions'''
        transfer_power_K = (1-A_T)*N_s_power[0]*create_lognorm(x_list,mu_power_K,sigma_power_K)
        transfer_power_H = (1-A_T)*N_s_power[1]*create_lognorm(x_list,mu_power_H,sigma_power_H)
        transfer_power_J = (1-A_T)*N_s_power[2]*create_lognorm(x_list,mu_power_J,sigma_power_J)
        transfer_power_g = (1-A_T)*N_s_power[3]*create_lognorm(x_list,mu_power_g,sigma_power_g)
        transfer_power_r = (1-A_T)*N_s_power[4]*create_lognorm(x_list,mu_power_r,sigma_power_r)
        transfer_power_i = (1-A_T)*N_s_power[5]*create_lognorm(x_list,mu_power_i,sigma_power_i)
        transfer_power_z = (1-A_T)*N_s_power[6]*create_lognorm(x_list,mu_power_z,sigma_power_z)

        #print transfer_power_K
        #print np.shape(transfer_power_K)

        '''Black Boby radiation'''
        BB_K = planck(K_cen,T)*K_width
        BB_H = planck(H_cen,T)*H_width
        BB_J = planck(J_cen,T)*J_width
        BB_g = planck(g_cen,T)*g_width
        BB_r = planck(r_cen,T)*r_width
        BB_i = planck(i_cen,T)*i_width
        BB_z = planck(z_cen,T)*z_width

        '''Black Body Transfer function'''
        transfer_thermal = A_T*create_lognorm(x_list,lag_thermal,width_thermal)

        transfer_K = transfer_power_K + transfer_thermal
        transfer_H = transfer_power_H + transfer_thermal
        transfer_J = transfer_power_J + transfer_thermal
        transfer_g = transfer_power_g + transfer_thermal
        transfer_r = transfer_power_r + transfer_thermal
        transfer_i = transfer_power_i + transfer_thermal
        transfer_z = transfer_power_z + transfer_thermal

        #print 'x_list', len(x_list) #, x_list

        '''The creation of the different transfer arrays'''
        transfer_array_K = transfer(cont_days,data_comp_K,transfer_K)
        transfer_array_H = transfer(cont_days,data_comp_H,transfer_H)
        transfer_array_J = transfer(cont_days,data_comp_J,transfer_J)
        transfer_array_g = transfer(cont_days,data_comp_g,transfer_g)
        transfer_array_r = transfer(cont_days,data_comp_r,transfer_r)
        transfer_array_i = transfer(cont_days,data_comp_i,transfer_i)
        transfer_array_z = transfer(cont_days,data_comp_z,transfer_z)
        
        '''The fitting of the various models'''
        data_comp_K[:,2] = model_data(cont_real,data_comp_K,transfer_array_K) #The data after running through the transfer function
        data_comp_H[:,2] = model_data(cont_real,data_comp_H,transfer_array_H)
        data_comp_J[:,2] = model_data(cont_real,data_comp_J,transfer_array_J)
        data_comp_g[:,2] = model_data(cont_real,data_comp_g,transfer_array_g)
        data_comp_r[:,2] = model_data(cont_real,data_comp_r,transfer_array_r)
        data_comp_i[:,2] = model_data(cont_real,data_comp_i,transfer_array_i)
        data_comp_z[:,2] = model_data(cont_real,data_comp_z,transfer_array_z)

        '''Determining the quality of the various fits'''
        chi2_K = np.nansum((data_comp_K[:,1] - data_comp_K[:,2])**2) #Finding residuals squared summed
        chi2_H = np.nansum((data_comp_H[:,1] - data_comp_H[:,2])**2)
        chi2_J = np.nansum((data_comp_J[:,1] - data_comp_J[:,2])**2)
        chi2_g = np.nansum((data_comp_g[:,1] - data_comp_g[:,2])**2)
        chi2_r = np.nansum((data_comp_r[:,1] - data_comp_r[:,2])**2)
        chi2_i = np.nansum((data_comp_i[:,1] - data_comp_i[:,2])**2)
        chi2_z = np.nansum((data_comp_z[:,1] - data_comp_z[:,2])**2)

        #print data_comp_K[:,2]

        chi2 = chi2_K + chi2_H + chi2_J + chi2_g + chi2_r + chi2_i + chi2_z

        #if try1 == 1.2 and chi2 - chi1 == 0.0:
        #    jump = 1*10**(-5)
        #    slope_jump = 0.1
        #    print 'reset'

        #print chi2 - chi1
        slopechange = (slope - slopeaim)**2 - (slopeolder1 - slopeaim)**2 #Is the new PSD slope a better fit.

        MCMC = random.random()
        #if (chi2 <= chi1 \
        #    and slopechange < 0.):
        #    '''The first instance of acceptance'''
        #    chi1 = chi2
        #    slopeolder1 = slope
            #print 1, slope
        #print chi2 - chi1

        if chi2 < chi1 \
           and (abs(slopeaim - slopeallow) < abs(slope) < abs(slopeaim + slopeallow) or slopechange <= 0) \
           and (data_comp_K[:,2] >= 0.).all() == True \
           and (data_comp_H[:,2] >= 0.).all() == True \
           and (data_comp_J[:,2] >= 0.).all() == True \
           and (data_comp_g[:,2] >= 0.).all() == True \
           and (data_comp_r[:,2] >= 0.).all() == True \
           and (data_comp_i[:,2] >= 0.).all() == True \
           and (data_comp_z[:,2] >= 0.).all() == True: #abs(slopeolder1) <
            '''The first instance of acceptance'''
            print float(i1)/float(runs1)*100, '%', try1, chi2 - chi1, slope

            chi1 = chi2
            slopeolder1 = slope
            if try1 == 0:
                try1 = 0
            elif try1 == 1.1:
                try1 = 2.1
            elif try1 == 1.2:
                try1 = 2.2
            elif try1 == 2.1 and random.random() < 0.01:
                try1 = 1.2
                #jump_change = 0
                #jump_change = 0
            elif try1 == 1.3:
                try1 = 2.3
            elif try1 == 2.2 and random.random() < 0.005:
                try1 = 1.3
            elif try1 == 2.3 and random.random() < 0.01:
                try1 = 0
                #jump_change = 0
            #elif try1 == 2.1 or try1 == 2.2 or try1 == 2.3:
            #    try1 = 0
            #    jump_change = 0
            #elif try1 == 2.2 and random.random() < 0.1:
            #    try1 = 0
            #    jump_change = 0


        elif 0.001 > MCMC \
           and (data_comp_K[:,2] >= 0.).all() == True \
           and (data_comp_H[:,2] >= 0.).all() == True \
           and (data_comp_J[:,2] >= 0.).all() == True \
           and (data_comp_g[:,2] >= 0.).all() == True \
           and (data_comp_r[:,2] >= 0.).all() == True \
           and (data_comp_i[:,2] >= 0.).all() == True \
           and (data_comp_z[:,2] >= 0.).all() == True: # \
             #and abs(slope) < abs(slopeaim + slopeallow): # and (dd2_ddt1 - dd2_ddt) < 0 and (dd3_ddt1 - dd3_ddt) < 0: # \
             #and random.random() < 0.1:
            '''The third instance of acceptance, an attempt to encourage a smoother curve.'''

            print 3, try1, chi2 - chi1, slope

            chi1 = chi2
            slopeolder1 = slope
            #print MCMC
            try1 = 0


        else:
            '''Rejection'''
            #print min(data_comp_K[:,2]),min(data_comp_H[:,2]),min(data_comp_J[:,2]),min(data_comp_g[:,2]), \
            #      min(data_comp_r[:,2]),min(data_comp_i[:,2]),min(data_comp_z[:,2])

            cont_real[:,1] = cont_real_save[:,1]
            cont[:,1] = cont_save[:,1]

            data_comp_K = data_comp_K_save
            data_comp_H = data_comp_H_save
            data_comp_J = data_comp_J_save
            data_comp_g = data_comp_g_save
            data_comp_r = data_comp_r_save
            data_comp_i = data_comp_i_save
            data_comp_z = data_comp_z_save

            T = T_save # K
            lag_thermal = lag_thermal_save
            width_thermal = width_thermal_save
            lag_slope_power = lag_slope_power_save
            lag_intercept_power = lag_intercept_power_save
            width_slope_power = width_slope_power_save
            width_intercept_power = width_intercept_power_save
            A_T = A_T_save
            N_s_power = N_s_power_save
            scale = scale_save

            mu_power_K = np.copy(mu_power_K_save)
            mu_power_H = np.copy(mu_power_H_save)
            mu_power_J = np.copy(mu_power_J_save)
            mu_power_g = np.copy(mu_power_g_save)
            mu_power_r = np.copy(mu_power_r_save)
            mu_power_i = np.copy(mu_power_i_save)
            mu_power_z = np.copy(mu_power_z_save)

            cont_fft = np.copy(cont_fft_save)
            
            sigma_power_K = np.copy(sigma_power_K_save)
            sigma_power_H = np.copy(sigma_power_H_save)
            sigma_power_J = np.copy(sigma_power_J_save)
            sigma_power_g = np.copy(sigma_power_g_save)
            sigma_power_r = np.copy(sigma_power_r_save)
            sigma_power_i = np.copy(sigma_power_i_save)
            sigma_power_z = np.copy(sigma_power_z_save)

            #print 'FAIL', try1, chi2 - chi1, slope

            if try1 == 0 and random.random() < 0.005:
                try1 = 1.1
            elif (try1 == 1.1 or try1 == 2.1): # and random.random() < 0.2: # and random.random() < 0.3:
                try1 = 1.2
            elif (try1 == 1.2 or try1 == 2.2) and random.random() < 0.2: # and random.random() < 0.5:
                try1 = 1.3
            elif (try1 == 1.3 or try1 == 2.3): # and random.random() < 0.2: # and random.random() < 0.5:
                try1 = 0



        #if chi2 < chi1:
        #    print chi2, chi1

        #model3 = model2
        gc.collect()

    print 'lag_thermal =', lag_thermal
    print 'width_thermal =', width_thermal
    print 'lag_slope_power =', lag_slope_power
    print 'lag_intercept_power =', lag_intercept_power
    print 'width_slope_power =', width_slope_power
    print 'width_intercept_power =', width_intercept_power
    print 'A_T =', A_T
    print 'N_s_power =', N_s_power
    print 'Temperature =', T[0]
    P_show = log(P_v)
    np.savetxt('CONTINUUM/NGC3783-continuum-slope-colour-2_5',cont) #Updating the files
    np.savetxt('CONTINUUM/NGC3783-continuum-T-slope-colour-2_5',np.array([cont_days,T]))
    np.savetxt('CONTINUUM/NGC3783-data_comp_K-slope-colour-2_5',data_comp_K)
    np.savetxt('CONTINUUM/NGC3783-data_comp_H-slope-colour-2_5',data_comp_H)
    np.savetxt('CONTINUUM/NGC3783-data_comp_J-slope-colour-2_5',data_comp_J)
    np.savetxt('CONTINUUM/NGC3783-data_comp_g-slope-colour-2_5',data_comp_g)
    np.savetxt('CONTINUUM/NGC3783-data_comp_r-slope-colour-2_5',data_comp_r)
    np.savetxt('CONTINUUM/NGC3783-data_comp_i-slope-colour-2_5',data_comp_i)
    np.savetxt('CONTINUUM/NGC3783-data_comp_z-slope-colour-2_5',data_comp_z)

    plt.figure()

    plt.plot(data_comp_K[:,0],data_comp_K[:,1],color='b')
    plt.scatter(data_comp_K[:,0],data_comp_K[:,2],color='b',s=3)

    plt.plot(data_comp_H[:,0],data_comp_H[:,1],color='r')
    plt.scatter(data_comp_H[:,0],data_comp_H[:,2],color='r',s=3)

    plt.plot(data_comp_J[:,0],data_comp_J[:,1],color='g')
    plt.scatter(data_comp_J[:,0],data_comp_J[:,2],color='g',s=3)

    plt.plot(data_comp_g[:,0],data_comp_g[:,1],color='purple')
    plt.scatter(data_comp_g[:,0],data_comp_g[:,2],color='purple',s=3)

    plt.plot(data_comp_r[:,0],data_comp_r[:,1],color='yellow')
    plt.scatter(data_comp_r[:,0],data_comp_r[:,2],color='yellow',s=3)

    plt.plot(data_comp_i[:,0],data_comp_i[:,1],color='black')
    plt.scatter(data_comp_i[:,0],data_comp_i[:,2],color='black',s=3)

    plt.plot(data_comp_z[:,0],data_comp_z[:,1],color='orange')
    plt.scatter(data_comp_z[:,0],data_comp_z[:,2],color='orange',s=3)

    plt.scatter(np.real(cont[:,0]),cont_real[:,1],color='brown',s=6)

    print 'cont =', cont[:,1][0]

    print 'cont_real =', cont_real[:,1][0]

    plt.ylim([1.2e-15,0.4e-13])
    print 'H', data_comp_H[:,1][10],data_comp_H[:,2][10]
    print 'J', data_comp_J[:,1][10],data_comp_J[:,2][10]
    print 'K', data_comp_K[:,1][10],data_comp_K[:,2][10]
    print 'g', data_comp_g[:,1][10],data_comp_g[:,2][10]
    print 'r', data_comp_r[:,1][10],data_comp_r[:,2][10]
    print 'i', data_comp_i[:,1][10],data_comp_i[:,2][10]
    print 'z', data_comp_z[:,1][10],data_comp_z[:,2][10]
    #plt.ylim([4e-15,1.2e-14])
    #plt.yscale('log')
    plt.show(block=False)

    plt.figure()

    plt.plot(data_comp_K[:,0],data_comp_K[:,1],color='b')
    plt.scatter(data_comp_K[:,0],data_comp_K[:,2],color='b')

    plt.plot(data_comp_H[:,0],data_comp_H[:,1],color='r')
    plt.scatter(data_comp_H[:,0],data_comp_H[:,2],color='r')

    plt.plot(data_comp_J[:,0],data_comp_J[:,1],color='g')
    plt.scatter(data_comp_J[:,0],data_comp_J[:,2],color='g')

    plt.plot(data_comp_g[:,0],data_comp_g[:,1],color='purple')
    plt.scatter(data_comp_g[:,0],data_comp_g[:,2],color='purple')

    plt.plot(data_comp_r[:,0],data_comp_r[:,1],color='yellow')
    plt.scatter(data_comp_r[:,0],data_comp_r[:,2],color='yellow')

    plt.plot(data_comp_i[:,0],data_comp_i[:,1],color='black')
    plt.scatter(data_comp_i[:,0],data_comp_i[:,2],color='black')

    plt.plot(data_comp_z[:,0],data_comp_z[:,1],color='orange')
    plt.scatter(data_comp_z[:,0],data_comp_z[:,2],color='orange')

    plt.scatter(np.real(cont[:,0]),cont_real[:,1],color='brown')

    print 'cont =', cont[:,1][0]

    print 'cont_real =', cont_real[:,1][0]

    plt.ylim([1.2e-18,0.4e-09])
    print 'H', data_comp_H[:,1][10],data_comp_H[:,2][10]
    print 'J', data_comp_J[:,1][10],data_comp_J[:,2][10]
    print 'K', data_comp_K[:,1][10],data_comp_K[:,2][10]
    print 'g', data_comp_g[:,1][10],data_comp_g[:,2][10]
    print 'r', data_comp_r[:,1][10],data_comp_r[:,2][10]
    print 'i', data_comp_i[:,1][10],data_comp_i[:,2][10]
    print 'z', data_comp_z[:,1][10],data_comp_z[:,2][10]
    #plt.ylim([4e-15,1.2e-14])
    plt.yscale('log')
    plt.show(block=False)

    plt.figure()
    #plt.plot(x_list,transfer_thermal/A_T,label='Thermal Transfer')
    plt.plot([1,2,3,4,5,6,7],N_s_power,label='N_s_power')
    plt.plot([1,2,3,4,5,6,7],[mu_power_K,mu_power_H,mu_power_J,mu_power_g,\
                              mu_power_r,mu_power_i,mu_power_z],label='Lag Power')
    plt.plot([1,2,3,4,5,6,7],[sigma_power_K,sigma_power_H,sigma_power_J,sigma_power_g,\
                              sigma_power_r,sigma_power_i,sigma_power_z],label='Sigma Power')
    plt.xticks([1,2,3,4,5,6,7],['K','H','J','g','r','i','z'])

    plt.legend()
    #plt.xlim([0,400])
    #plt.ylim([0,np.max(transfer_power_K)])
    plt.show(block=False)

    print slopeolder1
    end_time = time.time()
    print 'Time =', end_time - start_time
    print 'transfer_thermal', np.sum(transfer_thermal)
    print 'transfer_power', np.sum(transfer_power_K),np.sum(transfer_power_H),np.sum(transfer_power_J), \
          np.sum(transfer_power_g),np.sum(transfer_power_r),np.sum(transfer_power_i),np.sum(transfer_power_z)
    TCK.append(np.max(transfer_power_K))
    MK.append(mu_power_K)
    SK.append(sigma_power_K)
    print TCK
    print MK
    print SK
    print np.sum(transfer_power_K/N_s_power[0]), np.sum(transfer_thermal/A_T)

#print cont[:,1]
