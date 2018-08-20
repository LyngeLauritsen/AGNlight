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
from scipy.optimize import fminbound
#from scipy.optimize import least_squares as ls
from scipy.optimize import leastsq
import gc
from multiprocessing import Process
#import colorednoise as cn
from numpy import concatenate, real, std, abs, min
from numpy.fft import ifft, fftfreq, fft, rfft, irfft
from numpy.random import normal
from scipy.optimize import curve_fit
import time

def powerlaw_psd_gaussian(exponent, samples, mean, fmin=0):
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
    y = ifft(s).real
    y = y / std(y)
    y = y*10**np.random.normal(-14.5,1.5,1)[0]
    y += mean #abs(min(y))
    #y = y*10**np.random.normal(-13.7,1.5,1)[0]
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
K_cen = 2190#*10**(-9)   #Center of K band in m
H_cen = 1630#*10**(-9)   #Center of H band in m
J_cen = 1220#*10**(-9)   #Center of J band in m
g_cen = 474.7#*10**(-9)   #Center of g band in m
r_cen = 621.4#*10**(-9)   #Center of r band in m
i_cen = 762.8#*10**(-9)   #Center of i band in m
z_cen = 906.8#*10**(-9)   #Center of z band in m

K_width = 390#*10**(-9)  #Width of K band in m
H_width = 307#*10**(-9)  #Width of H band in m
J_width = 213#*10**(-9)  #Width of J band in m
g_width = 99.6#*10**(-9)  #Width of g band in m
r_width = 95.7#*10**(-9)  #Width of r band in m
i_width = 105.6#*10**(-9)  #Width of i band in m
z_width = 122.7#*10**(-9)  #Width of z band in m

runs = 1000 #3000 #Number of runs (arbitrary due to the gradual updates, would probably take close to a month to run)
runs1 = 10000 #1000 #1000
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
T = 1200. # K
lag_thermal = 3.6
width_thermal = .5
lag_slope_power = -0.0
lag_intercept_power = 3.
width_slope_power = -0.0
width_intercept_power = .5
A_T = 0.0005
N_s_power = np.array([.5,.5,.5,.5,.5,.5])
print N_s_power
scale = 12.5
N_s_a = 0.5
N_s_b = 0.

'''Fitting the constants'''
mu_power_K = 2.3
mu_power_H = 2.3
mu_power_J = 2.3
mu_power_g = 2.3
mu_power_r = 2.3
mu_power_i = 2.3
mu_power_z = 2.3

sigma_power_K = .5
sigma_power_H = .5
sigma_power_J = .5
sigma_power_g = .5
sigma_power_r = .5
sigma_power_i = .5
sigma_power_z = .5

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

def power (x,a,b):
    return a*x**b

def planck(wl,temperature):
    front = 2*h_c*c_c**2/wl**5
    exp = h_c*c_c/(wl*k_B*temperature)
    return front*1/(np.exp(exp)-1)

def colour(alpha,length,mean):
    return powerlaw_psd_gaussian(alpha,length,mean)

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

data_g = np.loadtxt('NOVEMBER/NOV-NGC3783-g') #'NOVEMBER/NOV-NGC3783-g_short')
error_g = np.loadtxt('NOVEMBER/NGC3783_NOISE_g.txt') #Load data

data_r = np.loadtxt('NOVEMBER/NOV-NGC3783-r') #'NOVEMBER/NOV-NGC3783-r_short')
error_r = np.loadtxt('NOVEMBER/NGC3783_NOISE_r.txt') #Load data

data_i = np.loadtxt('NOVEMBER/NOV-NGC3783-i') #'NOVEMBER/NOV-NGC3783-i_short')
error_i = np.loadtxt('NOVEMBER/NGC3783_NOISE_i.txt') #Load data

data_z = np.loadtxt('NOVEMBER/NOV-NGC3783-z') #'NOVEMBER/NOV-NGC3783-z_short')
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

print np.shape(cont), np.shape(cont_ifft)

time1 = cont_days
time1 = np.insert(time1,len(time1),np.max(time1)+5)
flux1 = cont_ifft
sigma1 = cont[:,2]
flux1 = np.insert(flux1,len(flux1),flux1[len(flux1)-1])
sigma1 = np.insert(sigma1,len(sigma1),sigma1[len(sigma1)-1])
#print sigma2

time2 = np.insert(time1,0,0)
sigma = np.insert(sigma1,0,np.mean(sigma1))
flux = np.insert(flux1,0,np.mean(flux1))

b = 8.*10**(-2)
tau = 300.
sigma_tot = 1*10**(-16)

param = [b,tau,sigma_tot]

print cont[:,1][0],cont[:,1][1],cont[:,1][2],cont[:,1][3],cont[:,1][4]
print cont_ifft[0],cont_ifft[1],cont_ifft[2],cont_ifft[3],cont_ifft[4]

#time.sleep(5)

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

#colour_start = colour(np.random.normal(2.7,0.15,1)[0],len(cont[:,1])) #colour(-slope,len(cont[:,1]))
#print colour_start[0]
#colour_start = colour_start*1e-17

#cont[:,0] = np.real(cont_days)
#!cont[:,1] = colour_start #colour(2.7,len(cont[:,1])) #np.random.randint(-1e5,1e5,len(cont_days))*1.e-22 #colour(-slope,len(cont[:,1]))*1e-17 #colour_start
#print 'Cont 1 =', cont[0,1]
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-Kelly-2_5-1')
#day = cont_days[0]

data_comp_K = np.zeros((len(data_K[:,1]),4))
data_comp_K[:,0] = data_K[:,0]
data_comp_K[:,1] = data_K[:,1]/K_width[0]
data_comp_K[:,2] = 0
data_comp_K[:,3] = error_K[:,1]/K_width[0]

mask = ~np.isnan(data_comp_K[:,1])
data_comp_K = data_comp_K[mask,:]
mask = ~np.isnan(data_comp_K[:,3])
data_comp_K = data_comp_K[mask,:]
data_comp_K = data_comp_K[np.ma.masked_where(0<data_comp_K[:,1],data_comp_K[:,1]).mask,:]


data_comp_H = np.zeros((len(data_H[:,1]),4))
data_comp_H[:,0] = data_H[:,0]
data_comp_H[:,1] = data_H[:,1]/H_width[0]
data_comp_H[:,2] = 0
data_comp_H[:,3] = error_H[:,1]/H_width[0]
data_comp_H = data_comp_H[np.ma.masked_where(0<data_comp_H[:,1],data_comp_H[:,1]).mask,:]

mask = ~np.isnan(data_comp_H[:,1])
data_comp_H = data_comp_H[mask,:]
mask = ~np.isnan(data_comp_H[:,3])
data_comp_H = data_comp_H[mask,:]

data_comp_J = np.zeros((len(data_J[:,1]),4))
data_comp_J[:,0] = data_J[:,0]
data_comp_J[:,1] = data_J[:,1]/J_width[0]
data_comp_J[:,2] = 0
data_comp_J[:,3] = error_J[:,1]/J_width[0]
data_comp_J = data_comp_J[np.ma.masked_where(0<data_comp_J[:,1],data_comp_J[:,1]).mask,:]

mask = ~np.isnan(data_comp_J[:,1])
data_comp_J = data_comp_J[mask,:]
mask = ~np.isnan(data_comp_J[:,3])
data_comp_J = data_comp_J[mask,:]

data_comp_g = np.zeros((len(data_g[:,1]),4))
data_comp_g[:,0] = data_g[:,0]
data_comp_g[:,1] = data_g[:,2]/g_width[0]
data_comp_g[:,2] = 0
data_comp_g[:,3] = error_g[:,2]/g_width[0]
data_comp_g = data_comp_g[np.ma.masked_where(0<data_comp_g[:,1],data_comp_g[:,1]).mask,:]

mask = ~np.isnan(data_comp_g[:,1])
data_comp_g = data_comp_g[mask,:]
mask = ~np.isnan(data_comp_g[:,3])
data_comp_g = data_comp_g[mask,:]

data_comp_r = np.zeros((len(data_r[:,1]),4))
data_comp_r[:,0] = data_r[:,0]
data_comp_r[:,1] = data_r[:,2]/r_width[0]
data_comp_r[:,2] = 0
data_comp_r[:,3] = error_r[:,2]/r_width[0]
data_comp_r = data_comp_r[np.ma.masked_where(0<data_comp_r[:,1],data_comp_r[:,1]).mask,:]
data_comp_r = data_comp_r[np.ma.masked_where(3e-14>data_comp_r[:,3],data_comp_r[:,3]).mask,:]

mask = ~np.isnan(data_comp_r[:,1])
data_comp_r = data_comp_r[mask,:]
mask = ~np.isnan(data_comp_r[:,3])
data_comp_r = data_comp_r[mask,:]

data_comp_i = np.zeros((len(data_i[:,1]),4))
data_comp_i[:,0] = data_i[:,0]
data_comp_i[:,1] = data_i[:,2]/i_width[0]
data_comp_i[:,2] = 0
data_comp_i[:,3] = error_i[:,2]/i_width[0]
data_comp_i = data_comp_i[np.ma.masked_where(0<data_comp_i[:,1],data_comp_i[:,1]).mask,:]

#print data_comp_i[:,3]
mask = ~np.isnan(data_comp_i[:,1])
data_comp_i = data_comp_i[mask,:]
mask = ~np.isnan(data_comp_i[:,3])
data_comp_i = data_comp_i[mask,:]
#print data_comp_i[:,3]

data_comp_z = np.zeros((len(data_z[:,1]),4))
data_comp_z[:,0] = data_z[:,0]
data_comp_z[:,1] = data_z[:,2]/z_width[0]
data_comp_z[:,2] = 0
data_comp_z[:,3] = error_z[:,2]/z_width[0]
data_comp_z = data_comp_z[np.ma.masked_where(0<data_comp_z[:,1],data_comp_z[:,1]).mask,:]

mask = ~np.isnan(data_comp_z[:,1])
data_comp_z = data_comp_z[mask,:]
mask = ~np.isnan(data_comp_z[:,3])
data_comp_z = data_comp_z[mask,:]

print np.mean(data_comp_K[:,3])
print np.mean(data_comp_H[:,3])
print np.mean(data_comp_J[:,3])
print np.mean(data_comp_g[:,3])
print np.mean(data_comp_r[:,3])
print np.mean(data_comp_i[:,3])
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
        delay = abs(int(np.real(cont[k,0])) - np.concatenate(([z,data_comp[:,0][mask]])))
        delay[np.isnan(delay)] = 0
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

def determine_constants(cont_days,mu_p,sigma_p,lag_thermal,width_thermal,N_s,BB,A_T,T,data_comp,width):
    transfer_power = (1-abs(A_T))*abs(N_s)*create_lognorm(x_list,mu_p,sigma_p)
    transfer_thermal = abs(A_T)*create_lognorm(x_list,lag_thermal,width_thermal)
    
    transfer1 = transfer_power + transfer_thermal
    transfer_array = transfer(cont_days,data_comp,transfer1)

    data_comp_special = np.copy(data_comp)
    data_comp_special[:,1] = data_comp_special[:,1]/width
    
    result = model_data(cont_real,data_comp,transfer_array)
    return result

'''x denotes the observed fluxes, and t the observation times, 
sigma is the error variances and p is the probability function'''

def a_i(t2,t1,tau):
    tau = float(tau)
    return np.exp(-(t2-t1)/tau)

def omega_0(tau,sigma_tot):
    return 0.5*tau*sigma_tot**2

def omega_i(omega_0,a_i,omega_before,sigma_before):
    return omega_0*(1-a_i**2) + a_i**2*omega_before*(1 - omega_before/(omega_before + sigma_before**2))

def x_hat_i(a_i,x_hat_before,omega_before,sigma_before,x_star_before):
    return a_i*x_hat_before + ((a_i*omega_before)/((omega_before + sigma_before**2)))*(x_star_before - x_hat_before)

def x_star_i(x_i,b,tau):
    return x_i - b*tau

def parameters(param):
    x_hat = [0]
    omega = [0.5*param[1]*param[2]**2]
    x_star = [-param[0]*param[1]]
    a = []
    for i in range(len(flux1)):
        a.append(a_i(time2[i+1],time2[i],param[1]))
        x_hat.append(x_hat_i(a[i],x_hat[i],omega[i],sigma[i],x_star[i]))
        omega.append(omega_i(omega_0(param[1],param[2]),a[i],omega[i],sigma[i]))
        x_star.append(x_star_i(flux[i+1],param[0],param[1]))
    x_hat[0] = np.mean(x_hat)
    x_hat[1] = np.mean(x_hat)
    return x_hat,omega,x_star,a

def prob(param):
    x_hat,omega,x_star,a = parameters(param)
    #print x_hat
    #print x_star
    #print omega
    
    probability = 1
    for i in range(len(a)):
        #print omega[i+1],sigma[i+1]**2 #,x_hat[i+1],x_star[i+1],a[i]
        #print (x_hat[i+1] - x_star[i+1])
        part = np.log(1/np.sqrt(2*np.pi*(omega[i+1] + sigma[i+1]**2)))\
        -(1/2.)*(((x_hat[i+1] - x_star[i+1])**2.)/(omega[i+1] + sigma[i+1]**2.))
        #print part,sigma1[i]
        #print ((x_hat[i+1] - x_star[i+1])**2.)/(omega[i+1] + sigma[i+1]**2.)
        #print 1/np.sqrt(2*np.pi*(omega[i+1] + sigma[i+1]**2))
        probability += part
        #print part, probability
    return -probability

def dX(tau,sigma,dt,b,epsilon,X):
    dt = dt
    sigma = sigma
    return (-(1/tau)*X*dt + sigma*np.sqrt(dt)*epsilon + b*dt)


res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
time_model = np.arange(time1[:1],time1[-1:],5)

flux_test = flux1[0]

def flux_model(param):
    model = np.zeros((1,len(time_model)))
    for j in range(len(model[:,0])):
        model[j,0] = flux1[0]
        #flux_model = [flux1[0]]
        #print j
        for i in range(len(time_model)-1):
            #print i, model[j,i]
            idx_model = np.copy(time1)
            idx_model[idx_model <= time_model[i]] = 0
            idx_time = (np.abs(idx_model - time_model[i])).argmin()
            dt = abs(time1[idx_time]-time_model[i])
            epsilon = np.random.normal(0,1,1)
            dX1 = (flux1[idx_time])#*(-1) + model[j,0]
            change = dX(tau,sigma_tot,dt,b,epsilon,dX1) #flux1[i+1] - flux1[i]
            #print change,flux1[i+1]
            model[j,i+1] = change + flux1[idx_time] #model[j,i] # - flux_test)
    model2 = []
    time_model1 = []
    std = []
    for i in range(len(model[0,:])):
        model2.append(np.mean(model[:,i]))
        std.append(np.std(model[:,i]))
        time_model1.append(time_model[i])
    flux_model = np.array((time_model1,model2,std,res.x))
    return flux_model,time_model1, model2

def flux_model1(param):
    model = np.zeros((1,len(time_model)))
    for j in range(len(model[:,0])):
        model[j,0] = flux1[0]
        #flux_model = [flux1[0]]
        #print j
        for i in range(len(time_model)-1):
            #print i, model[j,i]
            idx_model = np.copy(time1)
            idx_model[idx_model >= time_model[i]] = 0
            idx_time = (np.abs(idx_model - time_model[i])).argmin()
            dt = abs(time1[idx_time]-time_model[i])
            epsilon = np.random.normal(0,1,1)
            dX1 = (flux1[idx_time])#*(-1) + model[j,0]
            change = dX(tau,sigma_tot,dt,b,epsilon,dX1) #flux1[i+1] - flux1[i]
            #print change,flux1[i+1]
            model[j,i+1] = change + flux1[idx_time] #model[j,i] # - flux_test)
    model2 = []
    time_model1 = []
    std = []
    for i in range(len(model[0,:])):
        model2.append(np.mean(model[:,i]))
        std.append(np.std(model[:,i]))
        time_model1.append(time_model[i])
    flux_model = np.array((time_model1,model2,std,res.x))
    return flux_model,time_model1, model2

#model, time_model, model2 = flux_model(param)

change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
weight = 0.1 #random.random()/float(random.randint(50,500)) #0.05
cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))

change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
weight = 0.1 #random.random()/float(random.randint(50,500)) #0.05
cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))

change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
weight = 0.1 #random.random()/float(random.randint(50,500)) #0.05
cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))

x_list = np.linspace(0.01,len(cont_days)*step,len(cont_days)*step)

cont_real = np.real(cont)
cont_real[:,1] = irfft(cont_fft)

'''Black Boby radiation'''
BB_K = planck(K_cen,T)*K_width
BB_H = planck(H_cen,T)*H_width
BB_J = planck(J_cen,T)*J_width
BB_g = planck(g_cen,T)*g_width
BB_r = planck(r_cen,T)*r_width
BB_i = planck(i_cen,T)*i_width
#BB_z = planck(z_cen,T)*z_width

param_A = [N_s_a,N_s_b] #,A_T]
def prob_A(param_A):
    N_s_K = abs(power(K_cen[0],param_A[0],param_A[1]))
    N_s_H = abs(power(H_cen[0],param_A[0],param_A[1]))
    N_s_J = abs(power(J_cen[0],param_A[0],param_A[1]))
    N_s_g = abs(power(g_cen[0],param_A[0],param_A[1]))
    N_s_r = abs(power(r_cen[0],param_A[0],param_A[1]))
    N_s_i = abs(power(i_cen[0],param_A[0],param_A[1]))
    N_1 = determine_constants(cont_days,mu_power_K,sigma_power_K,lag_thermal,width_thermal,N_s_K,BB_K,A_T,T,data_comp_K,K_width[0])
    N_2 = determine_constants(cont_days,mu_power_H,sigma_power_H,lag_thermal,width_thermal,N_s_H,BB_H,A_T,T,data_comp_H,H_width[0])
    N_3 = determine_constants(cont_days,mu_power_J,sigma_power_J,lag_thermal,width_thermal,N_s_J,BB_H,A_T,T,data_comp_J,J_width[0])
    N_4 = determine_constants(cont_days,mu_power_g,sigma_power_g,lag_thermal,width_thermal,N_s_g,BB_H,A_T,T,data_comp_g,g_width[0])
    N_5 = determine_constants(cont_days,mu_power_r,sigma_power_r,lag_thermal,width_thermal,N_s_r,BB_H,A_T,T,data_comp_r,r_width[0])
    N_6 = determine_constants(cont_days,mu_power_i,sigma_power_i,lag_thermal,width_thermal,N_s_i,BB_H,A_T,T,data_comp_i,i_width[0])
    #print np.mean(N_1), N_s_K
    return np.sum((data_comp_K[:,1] - N_1)**2/data_comp_K[:,3]**2) + np.sum((data_comp_H[:,1] - N_2)**2/data_comp_H[:,3]**2) + \
           np.sum((data_comp_J[:,1] - N_3)**2/data_comp_J[:,3]**2) + np.sum((data_comp_g[:,1] - N_4)**2/data_comp_g[:,3]**2) + \
           np.sum((data_comp_r[:,1] - N_5)**2/data_comp_r[:,3]**2) + np.sum((data_comp_i[:,1] - N_6)**2/data_comp_i[:,3]**2)

res_A_T = mini(prob_A,param_A,method='Nelder-Mead',tol=1e-18)

N_s_a = res_A_T.x[0]
N_s_b = res_A_T.x[1]
#A_T = abs(res_A_T.x[2])

N_s_power[0] = abs(power(K_cen[0],N_s_a,N_s_b))
N_s_power[1] = abs(power(H_cen[0],N_s_a,N_s_b))
N_s_power[2] = abs(power(J_cen[0],N_s_a,N_s_b))
N_s_power[3] = abs(power(g_cen[0],N_s_a,N_s_b))
N_s_power[4] = abs(power(r_cen[0],N_s_a,N_s_b))
N_s_power[5] = abs(power(i_cen[0],N_s_a,N_s_b))

#N_s_power[6] = abs(popt_z[0])

print 'Constants =', N_s_power, A_T

#print N_s_power[0], abs(res_K.x[0])
#print N_s_power[1], abs(res_H.x[0])
#print N_s_power[2], abs(res_J.x[0])
#print N_s_power[3], abs(res_g.x[0])
#print N_s_power[4], abs(res_r.x[0])
#print N_s_power[5], abs(res_i.x[0])

N_1 = determine_constants(cont_days,mu_power_K,sigma_power_K,lag_thermal,width_thermal,N_s_power[0],BB_K,A_T,T,data_comp_K,K_width[0])
N_2 = determine_constants(cont_days,mu_power_H,sigma_power_H,lag_thermal,width_thermal,N_s_power[1],BB_H,A_T,T,data_comp_H,H_width[0])
N_3 = determine_constants(cont_days,mu_power_J,sigma_power_J,lag_thermal,width_thermal,N_s_power[2],BB_H,A_T,T,data_comp_J,J_width[0])
N_4 = determine_constants(cont_days,mu_power_g,sigma_power_g,lag_thermal,width_thermal,N_s_power[3],BB_H,A_T,T,data_comp_g,g_width[0])
N_5 = determine_constants(cont_days,mu_power_r,sigma_power_r,lag_thermal,width_thermal,N_s_power[4],BB_H,A_T,T,data_comp_r,r_width[0])
N_6 = determine_constants(cont_days,mu_power_i,sigma_power_i,lag_thermal,width_thermal,N_s_power[5],BB_H,A_T,T,data_comp_i,i_width[0])

print np.mean(data_comp_K[:,1]), np.mean(N_1)
print np.mean(data_comp_H[:,1]), np.mean(N_2)
print np.mean(data_comp_J[:,1]), np.mean(N_3)
print np.mean(data_comp_g[:,1]), np.mean(N_4)
print np.mean(data_comp_r[:,1]), np.mean(N_5)
print np.mean(data_comp_i[:,1]), np.mean(N_6)

'''Power transfer functions'''
transfer_power_K = (1-A_T)*N_s_power[0]*create_lognorm(x_list,mu_power_K,sigma_power_K)
transfer_power_H = (1-A_T)*N_s_power[1]*create_lognorm(x_list,mu_power_H,sigma_power_H)
transfer_power_J = (1-A_T)*N_s_power[2]*create_lognorm(x_list,mu_power_J,sigma_power_J)
transfer_power_g = (1-A_T)*N_s_power[3]*create_lognorm(x_list,mu_power_g,sigma_power_g)
transfer_power_r = (1-A_T)*N_s_power[4]*create_lognorm(x_list,mu_power_r,sigma_power_r)
transfer_power_i = (1-A_T)*N_s_power[5]*create_lognorm(x_list,mu_power_i,sigma_power_i)
#transfer_power_z = (1-A_T)*N_s_power[6]*create_lognorm(x_list,mu_power_z,sigma_power_z)

'''Black Body Transfer function'''
transfer_thermal = A_T*create_lognorm(x_list,lag_thermal,width_thermal)

transfer_K = transfer_power_K + transfer_thermal
transfer_H = transfer_power_H + transfer_thermal
transfer_J = transfer_power_J + transfer_thermal
transfer_g = transfer_power_g + transfer_thermal
transfer_r = transfer_power_r + transfer_thermal
transfer_i = transfer_power_i + transfer_thermal
#transfer_z = transfer_power_z + transfer_thermal


'''The creation of the different transfer arrays'''
transfer_array_K = transfer(cont_days,data_comp_K,transfer_K)
transfer_array_H = transfer(cont_days,data_comp_H,transfer_H)
transfer_array_J = transfer(cont_days,data_comp_J,transfer_J)
transfer_array_g = transfer(cont_days,data_comp_g,transfer_g)
transfer_array_r = transfer(cont_days,data_comp_r,transfer_r)
transfer_array_i = transfer(cont_days,data_comp_i,transfer_i)
#transfer_array_z = transfer(cont_days,data_comp_z,transfer_z)

data_comp_K[:,2] = model_data(cont_real,data_comp_K,transfer_array_K) #The data after running through the transfer function
data_comp_H[:,2] = model_data(cont_real,data_comp_H,transfer_array_H)
data_comp_J[:,2] = model_data(cont_real,data_comp_J,transfer_array_J)
data_comp_g[:,2] = model_data(cont_real,data_comp_g,transfer_array_g)
data_comp_r[:,2] = model_data(cont_real,data_comp_r,transfer_array_r)
data_comp_i[:,2] = model_data(cont_real,data_comp_i,transfer_array_i)

print np.mean(data_comp_K[:,1]), np.mean(data_comp_K[:,2])

print 'K', np.mean(data_comp_K[:,1]), np.mean(data_comp_K[:,2])
print 'H', np.mean(data_comp_H[:,1]), np.mean(data_comp_H[:,2])
print 'J', np.mean(data_comp_J[:,1]), np.mean(data_comp_J[:,2])
print 'g', np.mean(data_comp_g[:,1]), np.mean(data_comp_g[:,2])
print 'r', np.mean(data_comp_r[:,1]), np.mean(data_comp_r[:,2])
print 'i', np.mean(data_comp_i[:,1]), np.mean(data_comp_i[:,2])

#time.sleep(15)

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
chi1_K = np.nansum((data_comp_K[:,1] - data_comp_K[:,2])**2/data_comp_K[:,3]**2) #Finding residuals squared summed
chi1_H = np.nansum((data_comp_H[:,1] - data_comp_H[:,2])**2/data_comp_H[:,3]**2)
chi1_J = np.nansum((data_comp_J[:,1] - data_comp_J[:,2])**2/data_comp_J[:,3]**2)
chi1_g = np.nansum((data_comp_g[:,1] - data_comp_g[:,2])**2/data_comp_g[:,3]**2)
chi1_r = np.nansum((data_comp_r[:,1] - data_comp_r[:,2])**2/data_comp_r[:,3]**2)
chi1_i = np.nansum((data_comp_i[:,1] - data_comp_i[:,2])**2/data_comp_i[:,3]**2)
#chi1_z = np.nansum((data_comp_z[:,1] - data_comp_z[:,2])**2)

#chi1 = chi1_K + chi1_H + chi1_J + chi1_g + chi1_r + chi1_i + chi1_z


slopeolder1 = 0 #Temporary constant

F_N_v = np.zeros((len(freq)))
P_v = 0
#print F_N_v

#cont_real = np.real(cont)
cont_real_save = cont
cont_save = cont

data_comp_K_save = data_comp_K
data_comp_H_save = data_comp_H
data_comp_J_save = data_comp_J
data_comp_g_save = data_comp_g
data_comp_r_save = data_comp_r
data_comp_i_save = data_comp_i
#data_comp_z_save = data_comp_z
scale_save = scale



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
#mu_power_z_save = np.copy(mu_power_z)

sigma_power_K_save = np.copy(sigma_power_K)
sigma_power_H_save = np.copy(sigma_power_H)
sigma_power_J_save = np.copy(sigma_power_J)
sigma_power_g_save = np.copy(sigma_power_g)
sigma_power_r_save = np.copy(sigma_power_r)
sigma_power_i_save = np.copy(sigma_power_i)
#sigma_power_z_save = np.copy(sigma_power_z)

cont_fft_save = np.copy(cont_fft)

mu_power_K_direction = (-1)**np.random.randint(2,size=1)
mu_power_H_direction = (-1)**np.random.randint(2,size=1)
mu_power_J_direction = (-1)**np.random.randint(2,size=1)
mu_power_g_direction = (-1)**np.random.randint(2,size=1)
mu_power_r_direction = (-1)**np.random.randint(2,size=1)
mu_power_i_direction = (-1)**np.random.randint(2,size=1)
#mu_power_z_direction = (-1)**np.random.randint(2,size=1)

sigma_power_K_direction = (-1)**np.random.randint(2,size=1)
sigma_power_H_direction = (-1)**np.random.randint(2,size=1)
sigma_power_J_direction = (-1)**np.random.randint(2,size=1)
sigma_power_g_direction = (-1)**np.random.randint(2,size=1)
sigma_power_r_direction = (-1)**np.random.randint(2,size=1)
sigma_power_i_direction = (-1)**np.random.randint(2,size=1)
#sigma_power_z_direction = (-1)**np.random.randint(2,size=1)

TCK = []
SK = []
MK = []

zero_val = 0
one_val = 0
two_val = 0

thermal_list = []
K_list = np.zeros((runs,3+len(data_comp_K[:,0])))
H_list = np.zeros((runs,3+len(data_comp_H[:,0])))
J_list = np.zeros((runs,3+len(data_comp_J[:,0])))
g_list = np.zeros((runs,3+len(data_comp_g[:,0])))
r_list = np.zeros((runs,3+len(data_comp_r[:,0])))
i_list = np.zeros((runs,3+len(data_comp_i[:,0])))
z_list = np.zeros((runs,3+len(data_comp_z[:,0])))
transfer_list = []
continuous_list = np.zeros((runs,len(np.real(cont[:,0]))+len(cont_real[:,1])))

'''
plt.figure()
plt.scatter(data_comp_K[:,0],data_comp_K[:,1],color='b')
plt.plot(data_comp_K[:,0],N_1,color='b')

plt.scatter(data_comp_H[:,0],data_comp_H[:,1],color='r')
plt.plot(data_comp_H[:,0],N_1,color='r')

plt.scatter(data_comp_J[:,0],data_comp_J[:,1],color='g')
plt.plot(data_comp_J[:,0],N_1,color='g')

plt.scatter(data_comp_g[:,0],data_comp_g[:,1],color='black')
plt.plot(data_comp_g[:,0],N_1,color='black')

plt.scatter(data_comp_r[:,0],data_comp_r[:,1],color='yellow')
plt.plot(data_comp_r[:,0],N_1,color='yellow')

plt.scatter(data_comp_i[:,0],data_comp_i[:,1],color='orange')
plt.plot(data_comp_i[:,0],N_1,color='orange')
plt.ylim([1e-15,5e-14])
plt.show()
'''


start_time1 = time.time()

#print np.real(ifft(cont[:,1]).real), PSD(freq,np.real((ifft(cont[:,1]).real - np.min(ifft(cont[:,1]).real)))*10**(-scale))

for i_run in range(runs):
    '''The MCMC'''
    try1 = 3
    start_time = time.time()
    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
    b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
    sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
    time_model = np.arange(time1[:1],time1[-1:],5)
    model = np.zeros((1,len(time_model)))
    flux_test = flux1[0]

    #lag_thermal += (-1)**np.random.randint(2)*0.2*random.random()
    #width_thermal += (-1)**np.random.randint(2)*0.02*random.random()

    #mu_power_K += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_H += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_J += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_g += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_r += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_i += (-1)**np.random.randint(2)*0.1*random.random()
    #mu_power_z += (-1)**np.random.randint(2)*0.1*random.random()

    #sigma_power_K += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_H += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_J += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_g += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_r += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_i += (-1)**np.random.randint(2)*0.01*random.random()
    #sigma_power_z += (-1)**np.random.randint(2)*0.01*random.random()

    '''Changing the parameters'''
    change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
    weight = 0.005
    cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))
    
    #popt_K, pcov_K = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_K,sigma_power_K,lag_thermal,width_thermal,N_s_power,BB_K,A_T,T,data_comp_K),cont_days,data_comp_K[:,1])
    #popt_H, pcov_H = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_H,sigma_power_H,lag_thermal,width_thermal,N_s_power,BB_H,A_T,T,data_comp_H),cont_days,data_comp_H[:,1])
    #popt_J, pcov_J = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_J,sigma_power_J,lag_thermal,width_thermal,N_s_power,BB_J,A_T,T,data_comp_J),cont_days,data_comp_J[:,1])
    #popt_g, pcov_g = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_g,sigma_power_g,lag_thermal,width_thermal,N_s_power,BB_g,A_T,T,data_comp_g),cont_days,data_comp_g[:,1])
    #popt_r, pcov_r = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_r,sigma_power_r,lag_thermal,width_thermal,N_s_power,BB_r,A_T,T,data_comp_r),cont_days,data_comp_r[:,1])
    #popt_i, pcov_i = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_i,sigma_power_i,lag_thermal,width_thermal,N_s_power,BB_i,A_T,T,data_comp_i),cont_days,data_comp_i[:,1])
    #popt_z, pcov_z = curve_fit(lambda cont_days, N_s_power: determine_constants(cont_days,mu_power_z,sigma_power_z,lag_thermal,width_thermal,N_s_power,BB_z,A_T,T,data_comp_z),cont_days,data_comp_z[:,1])

    #N_s_power[0] = abs(popt_K[0])
    #N_s_power[1] = abs(popt_H[0])
    #N_s_power[2] = abs(popt_J[0])
    #N_s_power[3] = abs(popt_g[0])
    #N_s_power[4] = abs(popt_r[0])
    #N_s_power[5] = abs(popt_i[0])
    #N_s_power[6] = abs(popt_z[0])

    for i in range(4):
        print 'i =', i
        change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
        weight = 0.01
        cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))

    
    for i1 in range(runs1):
        '''Doing 2 loops to allow regular progress reports'''

        jump_change += 1
        slope_jump_change += 1
        #print jump
        if jump_change > 500: #250:
            #jump *= 0.7
            jump = 10**np.random.normal(-4,1,1)[0]
            jump_change = 0
            change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
            weight = 0.001
            cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))
            end_time1 = time.time()
            print float(i1)/float(runs1)*100, '%', 'jump =', jump, '0_val =', zero_val, '1_val =', \
                  one_val, '2_val =', two_val, 'slope =', slopeolder1, 'Time =', end_time1 - start_time1
            print 'lag thermal =', lag_thermal, 'Temperature =', T[0]
            start_time1 = time.time()
            zero_val = 0
            one_val = 0
            two_val = 0
            res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
            tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
            b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
            sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
            time_model = np.arange(time1[:1],time1[-1:],5)
            model = np.zeros((1,len(time_model)))
            flux_test = flux1[0]
            #param_A = [1.]
            
            #res_A_T = mini(prob_A,param_A,method='Nelder-Mead',tol=1e-18) #
            #A_T = abs(res_A_T.x[0])
            
            #print 'Constants =', A_T
        if slope_jump_change > 300000000:
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
        #data_comp_z_save = np.copy(data_comp_z)

        T_save = np.copy(T) # K
        lag_thermal_save = np.copy(lag_thermal)
        width_thermal_save = np.copy(width_thermal)
        lag_slope_power_save = np.copy(lag_slope_power)
        lag_intercept_power_save = np.copy(lag_intercept_power)
        width_slope_power_save = np.copy(width_slope_power)
        width_intercept_power_save = np.copy(width_intercept_power)
        A_T_save = np.copy(A_T)
        #N_s_power_save = np.copy(N_s_power)
        scale_save = np.copy(scale)

        cont_fft_save = np.copy(cont_fft)

        mu_power_K_save = np.copy(mu_power_K)
        mu_power_H_save = np.copy(mu_power_H)
        mu_power_J_save = np.copy(mu_power_J)
        mu_power_g_save = np.copy(mu_power_g)
        mu_power_r_save = np.copy(mu_power_r)
        mu_power_i_save = np.copy(mu_power_i)
        #mu_power_z_save = np.copy(mu_power_z)

        sigma_power_K_save = np.copy(sigma_power_K)
        sigma_power_H_save = np.copy(sigma_power_H)
        sigma_power_J_save = np.copy(sigma_power_J)
        sigma_power_g_save = np.copy(sigma_power_g)
        sigma_power_r_save = np.copy(sigma_power_r)
        sigma_power_i_save = np.copy(sigma_power_i)
        #sigma_power_z_save = np.copy(sigma_power_z)

        BB_K_save = np.copy(BB_K)
        BB_H_save = np.copy(BB_H)
        BB_J_save = np.copy(BB_J)
        BB_g_save = np.copy(BB_g)
        BB_r_save = np.copy(BB_r)
        BB_i_save = np.copy(BB_i)
        #BB_z_save = np.copy(BB_z)

        '''Power transfer functions'''
        transfer_power_K_save = np.copy(transfer_power_K)
        transfer_power_H_save = np.copy(transfer_power_H)
        transfer_power_J_save = np.copy(transfer_power_J)
        transfer_power_g_save = np.copy(transfer_power_g)
        transfer_power_r_save = np.copy(transfer_power_r)
        transfer_power_i_save = np.copy(transfer_power_i)
        #transfer_power_z_save = np.copy(transfer_power_z)
        
        '''Black Body Transfer function'''
        transfer_thermal_save = np.copy(transfer_thermal)
        
        transfer_K_save = np.copy(transfer_K)
        transfer_H_save = np.copy(transfer_H)
        transfer_J_save = np.copy(transfer_J)
        transfer_g_save = np.copy(transfer_g)
        transfer_r_save = np.copy(transfer_r)
        transfer_i_save = np.copy(transfer_i)
        #transfer_z_save = np.copy(transfer_z)
        
        '''The creation of the different transfer arrays'''
        transfer_array_K_save = np.copy(transfer_array_K)
        transfer_array_H_save = np.copy(transfer_array_H)
        transfer_array_J_save = np.copy(transfer_array_J)
        transfer_array_g_save = np.copy(transfer_array_g)
        transfer_array_r_save = np.copy(transfer_array_r)
        transfer_array_i_save = np.copy(transfer_array_i)
        #transfer_array_z_save = np.copy(transfer_array_z)

        if try1 == 0 or try1 == 1.1 or try1 == 2.1:
            change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
            weight = random.random()/float(random.randint(100,5000)) #0.05
            cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))
        if try1 == 0: # and random.random() < 0.05:
            #    #print i, i1/float(runs1), slopeolder1
            '''Changing the parameters'''
            flux1 = irfft(cont_fft)
            model, time_model, model2 = flux_model(param)
            model1, time_model1, model12 = flux_model1(param)
            
            time_model = np.add(time_model,time_model1)/2.
            model2 = np.add(model2,model12)/2.

            change = rfft(model2)
            
            #change = rfft(colour(np.random.normal(2.8,0.25,1)[0],len(cont[:,1]),np.nanmean(irfft(cont_fft)))) #change_size) #len(cont[:,1]))) rfft
            
            weight = random.random()/float(random.randint(2,2500)) #0.05
            cont_fft = np.exp((1 - weight)*np.log(cont_fft) + weight*np.log(change)) #origin2))
            #print try1
        if try1 == 1.1: # or try1 == 0:
            '''Thermal component'''
            T_direction = (-1)**np.random.randint(2,size=1) #np.array([(-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) for i in range(len(cont_days))]) # K
            lag_thermal_direction = (-1)**np.random.randint(2)
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal_direction = (-1)**np.random.randint(2)#*np.random.randint(5,size=1)

            T += T_direction*np.random.random()*5
            lag_thermal += lag_thermal_direction*random.random()*jump*5
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal += width_thermal_direction*random.random()*jump*5#*np.random.randint(5,size=1)
            
        #elif try1 == 1.2:
            '''Power component'''
            mu_power_K_direction = (-1)**np.random.randint(2)
            mu_power_H_direction = (-1)**np.random.randint(2)
            mu_power_J_direction = (-1)**np.random.randint(2)
            mu_power_g_direction = (-1)**np.random.randint(2)
            mu_power_r_direction = (-1)**np.random.randint(2)
            mu_power_i_direction = (-1)**np.random.randint(2)
            #mu_power_z_direction = (-1)**np.random.randint(2)
            
            sigma_power_K_direction = (-1)**np.random.randint(2)
            sigma_power_H_direction = (-1)**np.random.randint(2)
            sigma_power_J_direction = (-1)**np.random.randint(2)
            sigma_power_g_direction = (-1)**np.random.randint(2)
            sigma_power_r_direction = (-1)**np.random.randint(2)
            sigma_power_i_direction = (-1)**np.random.randint(2)
            #sigma_power_z_direction = (-1)**np.random.randint(2)

            mu_power_K += mu_power_K_direction*random.random()*jump
            mu_power_H += mu_power_H_direction*random.random()*jump
            mu_power_J += mu_power_J_direction*random.random()*jump
            mu_power_g += mu_power_g_direction*random.random()*jump
            mu_power_r += mu_power_r_direction*random.random()*jump
            mu_power_i += mu_power_i_direction*random.random()*jump
            #mu_power_z += mu_power_z_direction*random.random()*jump

            sigma_power_K += sigma_power_K_direction*random.random()*jump
            sigma_power_H += sigma_power_H_direction*random.random()*jump
            sigma_power_J += sigma_power_J_direction*random.random()*jump
            sigma_power_g += sigma_power_g_direction*random.random()*jump
            sigma_power_r += sigma_power_r_direction*random.random()*jump
            sigma_power_i += sigma_power_i_direction*random.random()*jump
            #sigma_power_z += sigma_power_z_direction*random.random()*jump

        #elif try1 == 1.3:
            #N_s_power_direction = (-1)**np.random.randint(2,size=7)
            A_T_direction = (-1)**np.random.randint(2)
            scale_direction = (-1)**np.random.randint(2)

            #N_s_power += N_s_power_direction*np.random.rand(7)*0.5*jump #.transpose()
            A_T += A_T_direction*random.random()*jump*0.01
            scale += scale_direction*random.random()/1000.

        if try1 == 2.1:
            T += (-1)**np.random.randint(2,size=1)*np.random.random()*5 #np.array([(-1)**np.random.randint(2,size=1)*random.random()*np.random.randint(5,size=1) for i in range(len(cont_days))]) # K
            #print np.random.randint(2)
            lag_thermal += lag_thermal_direction*random.random()*jump*5
            '''LOOK AT SOMETHING LIKE EMCEE IMPLEMENTATION IN lmfit(?)'''
            width_thermal += width_thermal_direction*random.random()*jump*5#*np.random.randint(5,size=1)

        #elif try1 == 2.2:
            mu_power_K += mu_power_K_direction*random.random()*jump
            mu_power_H += mu_power_H_direction*random.random()*jump
            mu_power_J += mu_power_J_direction*random.random()*jump
            mu_power_g += mu_power_g_direction*random.random()*jump
            mu_power_r += mu_power_r_direction*random.random()*jump
            mu_power_i += mu_power_i_direction*random.random()*jump
            #mu_power_z += mu_power_z_direction*random.random()*jump

            sigma_power_K += sigma_power_K_direction*random.random()*jump
            sigma_power_H += sigma_power_H_direction*random.random()*jump
            sigma_power_J += sigma_power_J_direction*random.random()*jump
            sigma_power_g += sigma_power_g_direction*random.random()*jump
            sigma_power_r += sigma_power_r_direction*random.random()*jump
            sigma_power_i += sigma_power_i_direction*random.random()*jump
            #sigma_power_z += sigma_power_z_direction*random.random()*jump
            
        #elif try1 == 2.3:
            A_T += A_T_direction*random.random()*jump*0.01
            #N_s_power += N_s_power_direction*np.random.rand(7)*jump #.transpose()
            scale += scale_direction*random.random()/1000.

        lag_intercept_power = abs(lag_intercept_power)
        width_intercept_power = abs(width_intercept_power)
        width_thermal = abs(width_thermal)
        lag_thermal = abs(lag_thermal)
        A_T = abs(A_T)

        #N_s_power = abs(N_s_power)
        '''Fitting the constants'''
        

        '''Looping over the course of the CLC'''
        F_N_v = np.zeros((len(freq))) #Creating the base array for the PSD
        cont_real[:,1] = irfft(cont_fft)#*10**(-scale) #Implementing change

        slope, P_v = PSD(freq,cont_real) #Finding PSD slope

        if try1 == 1.1 or try1 == 2.1: # or try1 == 0:

            '''Black Boby radiation'''
            BB_K = planck(K_cen,T)*K_width
            BB_H = planck(H_cen,T)*H_width
            BB_J = planck(J_cen,T)*J_width
            BB_g = planck(g_cen,T)*g_width
            BB_r = planck(r_cen,T)*r_width
            BB_i = planck(i_cen,T)*i_width
            #BB_z = planck(z_cen,T)*z_width

            #popt_A_T, pcov_A_T = curve_fit(lambda cont_days, A_T: determine_constants(cont_days,mu_power_K,sigma_power_K,lag_thermal,width_thermal,N_s_power[0],BB_K,A_T,T,data_comp_K),cont_days,data_comp_K[:,1])

            #A_T = abs(popt_A_T[0])

            '''Power transfer functions'''
            transfer_power_K = (1-A_T)*N_s_power[0]*create_lognorm(x_list,mu_power_K,sigma_power_K)
            transfer_power_H = (1-A_T)*N_s_power[1]*create_lognorm(x_list,mu_power_H,sigma_power_H)
            transfer_power_J = (1-A_T)*N_s_power[2]*create_lognorm(x_list,mu_power_J,sigma_power_J)
            transfer_power_g = (1-A_T)*N_s_power[3]*create_lognorm(x_list,mu_power_g,sigma_power_g)
            transfer_power_r = (1-A_T)*N_s_power[4]*create_lognorm(x_list,mu_power_r,sigma_power_r)
            transfer_power_i = (1-A_T)*N_s_power[5]*create_lognorm(x_list,mu_power_i,sigma_power_i)
            #transfer_power_z = (1-A_T)*N_s_power[6]*create_lognorm(x_list,mu_power_z,sigma_power_z)

            '''Black Body Transfer function'''
            transfer_thermal = A_T*create_lognorm(x_list,lag_thermal,width_thermal)

            transfer_K = transfer_power_K + transfer_thermal
            transfer_H = transfer_power_H + transfer_thermal
            transfer_J = transfer_power_J + transfer_thermal
            transfer_g = transfer_power_g + transfer_thermal
            transfer_r = transfer_power_r + transfer_thermal
            transfer_i = transfer_power_i + transfer_thermal
            #transfer_z = transfer_power_z + transfer_thermal
        

            '''The creation of the different transfer arrays'''
            transfer_array_K = transfer(cont_days,data_comp_K,transfer_K)
            transfer_array_H = transfer(cont_days,data_comp_H,transfer_H)
            transfer_array_J = transfer(cont_days,data_comp_J,transfer_J)
            transfer_array_g = transfer(cont_days,data_comp_g,transfer_g)
            transfer_array_r = transfer(cont_days,data_comp_r,transfer_r)
            transfer_array_i = transfer(cont_days,data_comp_i,transfer_i)
            #transfer_array_z = transfer(cont_days,data_comp_z,transfer_z)
        
        '''The fitting of the various models'''
        data_comp_K[:,2] = model_data(cont_real,data_comp_K,transfer_array_K) #The data after running through the transfer function
        data_comp_H[:,2] = model_data(cont_real,data_comp_H,transfer_array_H)
        data_comp_J[:,2] = model_data(cont_real,data_comp_J,transfer_array_J)
        data_comp_g[:,2] = model_data(cont_real,data_comp_g,transfer_array_g)
        data_comp_r[:,2] = model_data(cont_real,data_comp_r,transfer_array_r)
        data_comp_i[:,2] = model_data(cont_real,data_comp_i,transfer_array_i)
        #data_comp_z[:,2] = model_data(cont_real,data_comp_z,transfer_array_z)

        '''Determining the quality of the various fits'''
        chi2_K = np.nansum((data_comp_K[:,1] - data_comp_K[:,2])**2/data_comp_K[:,3]**2) #Finding residuals squared summed
        chi2_H = np.nansum((data_comp_H[:,1] - data_comp_H[:,2])**2/data_comp_H[:,3]**2)
        chi2_J = np.nansum((data_comp_J[:,1] - data_comp_J[:,2])**2/data_comp_J[:,3]**2)
        chi2_g = np.nansum((data_comp_g[:,1] - data_comp_g[:,2])**2/data_comp_g[:,3]**2)
        chi2_r = np.nansum((data_comp_r[:,1] - data_comp_r[:,2])**2/data_comp_r[:,3]**2)
        chi2_i = np.nansum((data_comp_i[:,1] - data_comp_i[:,2])**2/data_comp_i[:,3]**2)
        #chi2_z = np.nansum((data_comp_z[:,1] - data_comp_z[:,2])**2)

        #print data_comp_K[:,2]

        chi2 = chi2_K + chi2_H + chi2_J + chi2_g + chi2_r + chi2_i #+ chi2_z

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

        if 1 < chi2 < chi1 \
           and (data_comp_K[:,2] >= 0.).all() == True \
           and (data_comp_H[:,2] >= 0.).all() == True \
           and (data_comp_J[:,2] >= 0.).all() == True \
           and (data_comp_g[:,2] >= 0.).all() == True \
           and (data_comp_r[:,2] >= 0.).all() == True \
           and (data_comp_i[:,2] >= 0.).all() == True \
           and (data_comp_z[:,2] >= 0.).all() == True: #abs(slopeolder1) <
            #and (abs(slopeaim - slopeallow) < abs(slope) < abs(slopeaim + slopeallow) or slopechange <= 0)

            '''The first instance of acceptance'''
            #print float(i1)/float(runs1)*100, '%', try1, chi2 - chi1, slope

            chi1 = chi2
            slopeolder1 = slope
            if try1 == 3:
                print float(i1)/float(runs1)*100, '%', try1, chi2 - chi1, slope
            elif try1 == 0:
                try1 = 0
                zero_val += 1
            elif try1 == 1.1:
                try1 = 2.1
                one_val += 1
            elif try1 == 2.1:
                two_val += 1
            #elif try1 == 1.2:
            #    try1 = 2.2
            #elif try1 == 2.1 and random.random() < 0.01:
            #    try1 = 1.2
                #jump_change = 0
                #jump_change = 0
            #elif try1 == 1.3:
            #    try1 = 2.3
            #elif try1 == 2.2 and random.random() < 0.005:
            #    try1 = 1.3
            #elif try1 == 2.3 and random.random() < 0.01:
            #    try1 = 0
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
            #data_comp_z = data_comp_z_save

            T = T_save # K
            lag_thermal = lag_thermal_save
            width_thermal = width_thermal_save
            lag_slope_power = lag_slope_power_save
            lag_intercept_power = lag_intercept_power_save
            width_slope_power = width_slope_power_save
            width_intercept_power = width_intercept_power_save
            A_T = A_T_save
            #N_s_power = N_s_power_save
            scale = scale_save

            mu_power_K = np.copy(mu_power_K_save)
            mu_power_H = np.copy(mu_power_H_save)
            mu_power_J = np.copy(mu_power_J_save)
            mu_power_g = np.copy(mu_power_g_save)
            mu_power_r = np.copy(mu_power_r_save)
            mu_power_i = np.copy(mu_power_i_save)
            #mu_power_z = np.copy(mu_power_z_save)

            cont_fft = np.copy(cont_fft_save)
            
            sigma_power_K = np.copy(sigma_power_K_save)
            sigma_power_H = np.copy(sigma_power_H_save)
            sigma_power_J = np.copy(sigma_power_J_save)
            sigma_power_g = np.copy(sigma_power_g_save)
            sigma_power_r = np.copy(sigma_power_r_save)
            sigma_power_i = np.copy(sigma_power_i_save)
            #sigma_power_z = np.copy(sigma_power_z_save)

            BB_K = np.copy(BB_K_save)
            BB_H = np.copy(BB_H_save)
            BB_J = np.copy(BB_J_save)
            BB_g = np.copy(BB_g_save)
            BB_r = np.copy(BB_r_save)
            BB_i = np.copy(BB_i_save)
            #BB_z = np.copy(BB_z_save)

            '''Power transfer functions'''
            transfer_power_K = np.copy(transfer_power_K_save)
            transfer_power_H = np.copy(transfer_power_H_save)
            transfer_power_J = np.copy(transfer_power_J_save)
            transfer_power_g = np.copy(transfer_power_g_save)
            transfer_power_r = np.copy(transfer_power_r_save)
            transfer_power_i = np.copy(transfer_power_i_save)
            #transfer_power_z = np.copy(transfer_power_z_save)
            
            '''Black Body Transfer function'''
            transfer_thermal = np.copy(transfer_thermal_save)
            
            transfer_K = np.copy(transfer_K_save)
            transfer_H = np.copy(transfer_H_save)
            transfer_J = np.copy(transfer_J_save)
            transfer_g = np.copy(transfer_g_save)
            transfer_r = np.copy(transfer_r_save)
            transfer_i = np.copy(transfer_i_save)
            #transfer_z = np.copy(transfer_z_save)
            
            '''The creation of the different transfer arrays'''
            transfer_array_K = np.copy(transfer_array_K_save)
            transfer_array_H = np.copy(transfer_array_H_save)
            transfer_array_J = np.copy(transfer_array_J_save)
            transfer_array_g = np.copy(transfer_array_g_save)
            transfer_array_r = np.copy(transfer_array_r_save)
            transfer_array_i = np.copy(transfer_array_i_save)
            #transfer_array_z = np.copy(transfer_array_z_save)


            #print 'FAIL', try1, chi2 - chi1, slope

            if try1 == 3 and i1 > 1:
                try1 = 0
            elif try1 == 0: # and random.random() < 0.1:
                try1 = 1.1
            elif (try1 == 1.1 or try1 == 2.1) and random.random() < 0.01: # and random.random() < 0.1: # and random.random() < 0.2: # and random.random() < 0.3:
                try1 = 0
            #elif (try1 == 1.2 or try1 == 2.2) and random.random() < 0.2: # and random.random() < 0.5:
            #    try1 = 1.3
            #elif (try1 == 1.3 or try1 == 2.3): # and random.random() < 0.2: # and random.random() < 0.5:
            #    try1 = 0



        #if chi2 < chi1:
        #    print chi2, chi1

        #model3 = model2
        gc.collect()
    
    param_A = [N_s_a,N_s_b] #,A_T]
    
    res_A_T = mini(prob_A,param_A,method='Nelder-Mead',tol=1e-18)
    
    N_s_a = res_A_T.x[0]
    N_s_b = res_A_T.x[1]
    #A_T = abs(res_A_T.x[2])
    
    N_s_power[0] = abs(power(K_cen[0],N_s_a,N_s_b))
    N_s_power[1] = abs(power(H_cen[0],N_s_a,N_s_b))
    N_s_power[2] = abs(power(J_cen[0],N_s_a,N_s_b))
    N_s_power[3] = abs(power(g_cen[0],N_s_a,N_s_b))
    N_s_power[4] = abs(power(r_cen[0],N_s_a,N_s_b))
    N_s_power[5] = abs(power(i_cen[0],N_s_a,N_s_b))
    
    #N_s_power[6] = abs(popt_z[0])
    
    print 'Constants =', N_s_power, A_T
    
    
    #-----------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------

    '''Black Boby radiation'''
    BB_K = planck(K_cen,T)*K_width
    BB_H = planck(H_cen,T)*H_width
    BB_J = planck(J_cen,T)*J_width
    BB_g = planck(g_cen,T)*g_width
    BB_r = planck(r_cen,T)*r_width
    BB_i = planck(i_cen,T)*i_width
    #BB_z = planck(z_cen,T)*z_width
    
    #popt_A_T, pcov_A_T = curve_fit(lambda cont_days, A_T: determine_constants(cont_days,mu_power_K,sigma_power_K,lag_thermal,width_thermal,N_s_power[0],BB_K,A_T,T,data_comp_K),cont_days,data_comp_K[:,1])
    
    #A_T = abs(popt_A_T[0])
    
    '''Power transfer functions'''
    transfer_power_K = (1-A_T)*N_s_power[0]*create_lognorm(x_list,mu_power_K,sigma_power_K)
    transfer_power_H = (1-A_T)*N_s_power[1]*create_lognorm(x_list,mu_power_H,sigma_power_H)
    transfer_power_J = (1-A_T)*N_s_power[2]*create_lognorm(x_list,mu_power_J,sigma_power_J)
    transfer_power_g = (1-A_T)*N_s_power[3]*create_lognorm(x_list,mu_power_g,sigma_power_g)
    transfer_power_r = (1-A_T)*N_s_power[4]*create_lognorm(x_list,mu_power_r,sigma_power_r)
    transfer_power_i = (1-A_T)*N_s_power[5]*create_lognorm(x_list,mu_power_i,sigma_power_i)
    #transfer_power_z = (1-A_T)*N_s_power[6]*create_lognorm(x_list,mu_power_z,sigma_power_z)
    
    '''Black Body Transfer function'''
    transfer_thermal = A_T*create_lognorm(x_list,lag_thermal,width_thermal)
    
    transfer_K = transfer_power_K + transfer_thermal
    transfer_H = transfer_power_H + transfer_thermal
    transfer_J = transfer_power_J + transfer_thermal
    transfer_g = transfer_power_g + transfer_thermal
    transfer_r = transfer_power_r + transfer_thermal
    transfer_i = transfer_power_i + transfer_thermal
    #transfer_z = transfer_power_z + transfer_thermal
    
    
    '''The creation of the different transfer arrays'''
    transfer_array_K = transfer(cont_days,data_comp_K,transfer_K)
    transfer_array_H = transfer(cont_days,data_comp_H,transfer_H)
    transfer_array_J = transfer(cont_days,data_comp_J,transfer_J)
    transfer_array_g = transfer(cont_days,data_comp_g,transfer_g)
    transfer_array_r = transfer(cont_days,data_comp_r,transfer_r)
    transfer_array_i = transfer(cont_days,data_comp_i,transfer_i)
    #transfer_array_z = transfer(cont_days,data_comp_z,transfer_z)
    
    '''The fitting of the various models'''
    data_comp_K[:,2] = model_data(cont_real,data_comp_K,transfer_array_K) #The data after running through the transfer function
    data_comp_H[:,2] = model_data(cont_real,data_comp_H,transfer_array_H)
    data_comp_J[:,2] = model_data(cont_real,data_comp_J,transfer_array_J)
    data_comp_g[:,2] = model_data(cont_real,data_comp_g,transfer_array_g)
    data_comp_r[:,2] = model_data(cont_real,data_comp_r,transfer_array_r)
    data_comp_i[:,2] = model_data(cont_real,data_comp_i,transfer_array_i)
    #data_comp_z[:,2] = model_data(cont_real,data_comp_z,transfer_array_z)
    #-----------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------
        

    thermal_list.append([lag_thermal,width_thermal,T[0],A_T])
    K_list[i_run,0],K_list[i_run,1],K_list[i_run,2],K_list[i_run,3:] = mu_power_K,sigma_power_K,N_s_power[0],data_comp_K[:,2]
    H_list[i_run,0],H_list[i_run,1],H_list[i_run,2],H_list[i_run,3:] = mu_power_H,sigma_power_H,N_s_power[1],data_comp_H[:,2]
    J_list[i_run,0],J_list[i_run,1],J_list[i_run,2],J_list[i_run,3:] = mu_power_J,sigma_power_J,N_s_power[2],data_comp_J[:,2]
    g_list[i_run,0],g_list[i_run,1],g_list[i_run,2],g_list[i_run,3:] = mu_power_g,sigma_power_g,N_s_power[3],data_comp_g[:,2]
    r_list[i_run,0],r_list[i_run,1],r_list[i_run,2],r_list[i_run,3:] = mu_power_r,sigma_power_r,N_s_power[4],data_comp_r[:,2]
    i_list[i_run,0],i_list[i_run,1],i_list[i_run,2],i_list[i_run,3:] = mu_power_i,sigma_power_i,N_s_power[5],data_comp_i[:,2]
    transfer_list.append([np.sum(transfer_power_K),np.sum(transfer_power_H),np.sum(transfer_power_J),\
                          np.sum(transfer_power_g),np.sum(transfer_power_r),np.sum(transfer_power_i),\
                          np.sum(transfer_thermal)])
    continuous_list[i_run,:len(np.real(cont[:,0]))], continuous_list[i_run,len(np.real(cont[:,0])):] = np.real(cont[:,0]),cont_real[:,1]
    np.savetxt('long_run_2/NGC3783/thermal_list_1',np.array(thermal_list))
    np.savetxt('long_run_2/NGC3783/K_list_1',K_list)
    np.savetxt('long_run_2/NGC3783/H_list_1',H_list)
    np.savetxt('long_run_2/NGC3783/J_list_1',J_list)
    np.savetxt('long_run_2/NGC3783/g_list_1',g_list)
    np.savetxt('long_run_2/NGC3783/r_list_1',r_list)
    np.savetxt('long_run_2/NGC3783/i_list_1',i_list)
    np.savetxt('long_run_2/NGC3783/transfer_list_1',np.array(transfer_list))
    np.savetxt('long_run_2/NGC3783/cont_list_1',continuous_list)

    print 'K', np.mean(data_comp_K[:,1]), np.mean(data_comp_K[:,2])
    print 'H', np.mean(data_comp_H[:,1]), np.mean(data_comp_H[:,2])
    print 'J', np.mean(data_comp_J[:,1]), np.mean(data_comp_J[:,2])
    print 'g', np.mean(data_comp_g[:,1]), np.mean(data_comp_g[:,2])
    print 'r', np.mean(data_comp_r[:,1]), np.mean(data_comp_r[:,2])
    print 'i', np.mean(data_comp_i[:,1]), np.mean(data_comp_i[:,2])

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
    
    print 'cont =', cont[:,1][0]

    print 'cont_real =', cont_real[:,1][0]

    
    print 'H', data_comp_H[:,1][10],data_comp_H[:,2][10]
    print 'J', data_comp_J[:,1][10],data_comp_J[:,2][10]
    print 'K', data_comp_K[:,1][10],data_comp_K[:,2][10]
    print 'g', data_comp_g[:,1][10],data_comp_g[:,2][10]
    print 'r', data_comp_r[:,1][10],data_comp_r[:,2][10]
    print 'i', data_comp_i[:,1][10],data_comp_i[:,2][10]
    #print 'z', data_comp_z[:,1][10],data_comp_z[:,2][10]
    #plt.ylim([4e-15,1.2e-14])
    #plt.yscale('log')

    print 'cont =', cont[:,1][0]

    print 'cont_real =', cont_real[:,1][0]

    #plt.ylim([1.2e-18,0.4e-09])
    print 'H', data_comp_H[:,1][10],data_comp_H[:,2][10]
    print 'J', data_comp_J[:,1][10],data_comp_J[:,2][10]
    print 'K', data_comp_K[:,1][10],data_comp_K[:,2][10]
    print 'g', data_comp_g[:,1][10],data_comp_g[:,2][10]
    print 'r', data_comp_r[:,1][10],data_comp_r[:,2][10]
    print 'i', data_comp_i[:,1][10],data_comp_i[:,2][10]
    #print 'z', data_comp_z[:,1][10],data_comp_z[:,2][10]
    #plt.ylim([4e-15,1.2e-14])
    #plt.yscale('log')
    

    print 'mu_power_K =', mu_power_K, 'sigma_power_K =', sigma_power_K
    print 'mu_power_H =', mu_power_H, 'sigma_power_H =', sigma_power_H
    print 'mu_power_J =', mu_power_J, 'sigma_power_J =', sigma_power_J
    print 'mu_power_g =', mu_power_g, 'sigma_power_g =', sigma_power_g
    print 'mu_power_r =', mu_power_r, 'sigma_power_r =', sigma_power_r
    print 'mu_power_i =', mu_power_i, 'sigma_power_i =', sigma_power_i
    #print 'mu_power_z =', mu_power_z, 'sigma_power_z =', sigma_power_z

    print 'lag_thermal =', lag_thermal, 'width_thermal =', width_thermal
    
    print slopeolder1
    end_time = time.time()
    print 'Time =', end_time - start_time
    print 'transfer_thermal', np.sum(transfer_thermal)
    print 'transfer_power', np.sum(transfer_power_K),np.sum(transfer_power_H),np.sum(transfer_power_J), \
          np.sum(transfer_power_g),np.sum(transfer_power_r),np.sum(transfer_power_i)
    TCK.append(np.max(transfer_power_K))
    MK.append(mu_power_K)
    SK.append(sigma_power_K)
    print TCK
    print MK
    print SK
    print np.sum(transfer_power_K/N_s_power[0]), np.sum(transfer_thermal/A_T)

#print cont[:,1]
