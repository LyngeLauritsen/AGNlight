import numpy as np
from numpy import *
import scipy as sc
import scipy.optimize as opt
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from numpy.random import randint
import math
#from scipy.optimize import minimize as mini
from scipy.optimize import least_squares as ls
from scipy.optimize import leastsq
import gc
from multiprocessing import Process
from lmfit import minimize, Parameters



runs = 3000 #Number of runs (arbitrary due to the gradual updates, would probably take close to a month to run)
runs1 = 10
days_before = 200. #How far to extend the lightcurve back (days)
days_after = 100. #How far to extend the lightcurve forward (days)
step = 10. #Timestep (days)
slopeaim = -2.5 #Slope that is targeted
slope = 0
mean = 140. #Transfer function
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma
kelly_repeat = 20

# There are signal processing package in scipy that defines filters, prob. can find lognormal in there.
def lognorm(x,mu,sigma):
    '''Defines the transfer function'''
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)


x_list = np.linspace(0.01,1500,1500) #X-axis of transfer function



log_norm_transfer = []
for i in range(len(x_list)):
    '''Creates the list used to identify the degree of transfered light'''
    log_norm_transfer.append(lognorm(x_list[i],mu,sigma))
#Suggest:
# log = lognorm(np.arange(len(x_list)))

print np.nansum(log_norm_transfer)*1500/1500.

data = np.loadtxt('NOVEMBER/NOV-NGC3783-K') #'NOVEMBER/Kelly-NGC3783K')
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt') #Load data

#plt.figure()
#plt.scatter(data[:,0],data[:,1])
#plt.ylim([6e-15,1.1e-14])

plt.show()

dd_ddt_data = []
for j in range(len(data[:,0])-2):
    '''Purely a test to see the largest double derivative of the observed light curve'''
    dd_ddt_data.append(abs(((data[j+1,1] - data[j,1])/(data[j+1,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+1,0] - data[j-1,0])/2.))

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


dd_ddt_allowed = np.nanmax(dd_ddt_data)*1.4 #Identifies the maximal double derivative allowed in the continuum light curve


dd_ddt_allowed = np.nanmax(dd_ddt_data)*1.4 #Identifies the maximal double derivative allowed in the continuum light curve

'''cont denotes the model for the continuum light curve, whereas data_comp denotes the array giving indication of the
observed light curve and how the continuum light curve fitted through the transfer function compares'''
cont_days = np.arange(min(data[:,0])-days_before,max(data[:,0]+days_after),step) #Defines the spacing of the light curve
#cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]+100),step) #Defines the spacing of the light curve
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = np.nanmean(data[:,1])
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-2-5-3')
day = cont_days[0]

data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0
#data_comp = np.loadtxt('CONTINUUM/NGC3783-data_comp-slope-2-5-3')
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

'''The set of equations that forms the basis for the Kelly function (see second Kelly 2009 paper in the PDF
for the full set of equations explained.'''

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

def parameters(b,tau,sigma_tot):
    x_hat = [0]
    omega = [0.5*tau*sigma_tot**2]
    x_star = [-b*tau]
    a = []
    for i in range(len(flux1)):
        a.append(a_i(time[i+1],time[i],tau))
        x_hat.append(x_hat_i(a[i],x_hat[i],omega[i],sigma[i],x_star[i]))
        omega.append(omega_i(omega_0(tau,sigma_tot),a[i],omega[i],sigma[i]))
        x_star.append(x_star_i(flux[i+1],b,tau))
    x_hat[0] = np.mean(x_hat)
    x_hat[1] = np.mean(x_hat)
    return x_hat,omega,x_star,a

def prob(params):
    b = params['b']
    tau = params['tau']
    sigma_tot = params['sigma_tot']
    x_hat,omega,x_star,a = parameters(b,tau,sigma_tot)
    probability = 1
    for i in range(len(a)):
        part = np.log(1/np.sqrt(2*np.pi*(omega[i+1] + sigma[i+1]**2)))\
        -(1/2.)*(((x_hat[i+1] - x_star[i+1])**2.)/(omega[i+1] + sigma[i+1]**2.))
        probability += part #just use np.sum()
    return -probability

def dX(tau,sigma,dt,b,epsilon,X):
    dt = dt
    sigma = sigma
    term1 = -(1/tau)*X*dt
    term2 = sigma*np.sqrt(dt)*epsilon
    term3 =  b*dt
    return term1 + term2 + term3

def flux_model(b,tau,sigma_tot):
    '''Creating the kelly function'''
    #res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    params = Parameters()
    #tau = param[1] #res.x[1]
    #b = param[0] #res.x[0]
    #sigma_tot = param[2] #res.x[2]
    params.add('tau',value=tau)
    params.add('sigma_tot',value=sigma_tot)
    params.add('b',value=b)
    print b, tau, sigma_tot
    res = minimize(prob,params,args=())
    print b, tau, sigma_tot
    model = np.zeros((kelly_repeat,len(flux1)))
    flux_test = flux1[0]
    dx_day = np.zeros((kelly_repeat,len(flux1)))
    flux_array = np.array([flux,]*(kelly_repeat))
    dt = np.array([abs(np.diff(time)),]*kelly_repeat)
    epsilon = np.random.normal(0,1,size=np.shape(dt))
    X1 = np.array([flux1[:],]*kelly_repeat) #flux1[1:]
    change = dX(tau,sigma_tot,dt,b,epsilon,X1)
    model = flux_array[:,:-1] + change
    model2 = np.mean(model,axis=0)
    std = np.std(model,axis=0)
    time_model = time1
    dx_days = np.mean(dx_day,axis=0)
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model,time_model, model2, dx_days

def flux_model1(param):
    '''Forming the Kelly function from in reverse'''
    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1]
    b = res.x[0]
    timef = time[::-1]
    fluxf = flux[::-1]
    time1f = time1[::-1]
    flux1f = flux1[::-1]
    sigma_tot = res.x[2]
    model = np.zeros((kelly_repeat,len(flux1f)))
    flux_test = flux1f[0]
    dx_day = np.zeros((kelly_repeat,len(flux1f)))
    flux_array = np.array([fluxf,]*(kelly_repeat))
    dt = np.array([abs(np.diff(timef)),]*kelly_repeat)
    epsilon = np.random.normal(0,1,size=np.shape(dt))
    X1 = np.array([flux1f[:],]*kelly_repeat) #flux1[1:]
    change = dX(tau,sigma_tot,dt,b,epsilon,X1)
    model = flux_array[:,:-1] + change
    model2 = np.mean(model,axis=0)
    std = np.std(model,axis=0)
    time_model = time1f
    dx_days = np.mean(dx_day,axis=0)
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model[::-1],time_model[::-1], model2[::-1], dx_days[::-1]


'''
def flux_model1(param):

    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1]
    b = res.x[0]
    sigma_tot = res.x[2]
    model = np.zeros((kelly_repeat,len(flux1)))
    flux_test = flux1[0]
    dx_day = np.zeros((kelly_repeat,len(flux1)))
    timef = time1[::-1]
    fluxf = flux1[::-1]
    flux_array = np.array([fluxf[:],]*(kelly_repeat))
    dt = np.array([abs(np.diff(timeff)),]*kelly_repeat)
    print shape(flux_array), shape(dt)
    epsilon = np.random.normal(0,1,size=np.shape(dt))
    X1 = np.array([fluxf[:],]*kelly_repeat) #fluxf[1:]
    change = dX(tau,sigma_tot,dt,b,epsilon,X1)
    model = flux_array[:,:-1] + change
    model2 = np.mean(model,axis=0)
    std = np.std(model,axis=0)
    time_model = timef
    dx_days = np.mean(dx_day,axis=0)
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model,time_model[::-1], model2[::-1], dx_days[::-1]
'''

freq = []

for i in range(int(len(cont_days)*0.5)):
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
    /(np.nanmean(cont[:,1])**2*len(cont[:,0])**2)*F_N_v**2 #Finding PSD

    slope = np.polyfit(np.log10(freq),np.log10(P_v),1)[0] #Finding PSD slope

    return slope

def model_data(cont,data_comp,transfer_array):
    '''The simulated observed light curve as a result of the continuum light curve and the transfer function'''

    cont_flux = array([cont[:,1],]*len(data_comp[:,0]))

    return sum(step*cont_flux*transfer_array,axis=1)

model3 = cont[:,1]
tau = 800. #The Kelly relaxation time
sigma_tot = np.nanmean(error[:,1]) #Kelly sigma
chi1 = 1e100 #To ensure that the first change is always accepted to get it started
max_jump = 0.01 #Maximal allowed jump relative to the value at point

slopeolder1 = 0 #Temporary constant

F_N_v = np.zeros((len(freq)))
P_v = 0
print F_N_v

for i in range(runs):
    '''The MCMC'''
    for i1 in range(runs1):
        '''Doing 2 loops to allow regular progress reports'''
        print i, i1/float(runs1), slopeolder1
        for j in range(len(cont_days)-3):
            '''Looping over the course of the CLC'''
            F_N_v = np.zeros((len(freq))) #Creating the base array for the PSD
            data_comp[:,2] = 0 #Resetting the simulated OLC
            dd1_ddt = abs(((cont[j+1,1] - cont[j,1])/step - (cont[j,1] - cont[j-1,1])/step)/step/2.) #Creating the double derivatives of the original CLC
            dd2_ddt = abs(((cont[j,1] - cont[j-1,1])/step - (cont[j-1,1] - cont[j-2,1])/step)/step/2.) #Doing so for the three points surrounding the change
            dd3_ddt = abs(((cont[j+2,1] - cont[j+1,1])/step - (cont[j+1,1] - cont[j,1])/step)/step/2.)

            startpoint1 = abs(cont[j-1,1]-model3[j-1])/cont[j-1,1]/step #Finding the relative difference between the original CLC and original Kelly 2009
            startpoint2 = abs(cont[j,1]-model3[j])/cont[j,1]/step
            startpoint3 = abs(cont[j+1,1]-model3[j+1])/cont[j+1,1]/step

            change = (-1)**randint(0,2)*cont[j,1]*random.random()*max_jump #Finding the magnitude and direction of change
            cont[j,1] += change #Implementing change

            h = 0
            b = np.nanmean(cont[:,1]) #Finding mean value of new Kelly 2009
            flux1 = cont[:,1] #Updating global variable for Kelly 2009
            param = [b,tau,sigma_tot] #See above

            '''p1 = Process(target=flux_model(param))
            p2 = Process(target=flux_model1(param))
            p1.start()
            p2.start()
            model, time_model, model2, dx_days = p1.join()
            model1, time_model1, model12, dx_days1 = p2.join()''' #Unsuccessfull attempt at parrallisation

            model, time_model, model2, dx_days = flux_model(b,tau,sigma_tot) #Finding forward and backward Kelly 2009
            model1, time_model1, model12, dx_days1 = flux_model1(param)

            #print shape(model2), shape(time_model)

            time_model = np.add(time_model,time_model1)/2. #Finding final Kelly 2009
            model2 = np.add(model2,model12)/2. #See above


            model21 = model2[:] - np.nanmean(model2[:])

            slope = PSD(freq,transfer_array,cont) #Finding PSD slope

            data_comp[:,2] = model_data(cont,data_comp,transfer_array) #The data after running through the transfer function

            chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2) #Finding residuals squared summed

            endpoint1 = abs(cont[j-1,1]-model2[j-1])/cont[j-1,1]/step #See startpoint
            endpoint2 = abs(cont[j,1]-model2[j])/cont[j,1]/step
            endpoint3 = abs(cont[j+1,1]-model2[j+1])/cont[j+1,1]/step

            dd1_ddt1 = abs(((cont[j+1,1] - cont[j,1])/step - (cont[j,1] - cont[j-1,1])/step)/step/2.) #See dd1_ddt
            dd2_ddt1 = abs(((cont[j,1] - cont[j-1,1])/step - (cont[j-1,1] - cont[j-2,1])/step)/step/2.)
            dd3_ddt1 = abs(((cont[j+2,1] - cont[j+1,1])/step - (cont[j+1,1] - cont[j,1])/step)/step/2.)

            slopechange = (slope - slopeaim)**2 - (slopeolder1 - slopeaim)**2 #Is the new PSD slope a better fit.

            if chi2 <= chi1 and slopechange < 0. and dd1_ddt1 < dd_ddt_allowed and dd2_ddt1 < dd_ddt_allowed \
               and dd3_ddt1 < dd_ddt_allowed: #\ #and endpoint1 < 0.0055/2. and endpoint2 < 0.0055/2. and endpoint3 < 0.0055/2.
               # and slopechange < 0.:
                '''The first instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 1, slope, dd1_ddt1

            elif chi2 <= chi1 and slopeaim - 0.05 < slope < slopeaim + 0.05 and dd1_ddt1 < dd_ddt_allowed \
                 and dd2_ddt1 < dd_ddt_allowed and dd3_ddt1 < dd_ddt_allowed: #and endpoint1 < 0.0055/2. and endpoint2 < 0.0055/2. and endpoint3 < 0.0055/2. \
               #and dd_ddt1 < 3.e-17 and slopeaim - 0.05 < slope < slopeaim + 0.05:
                '''The second instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 2, slope

            elif chi2 <= chi1 and (startpoint1 - endpoint1) + (startpoint2 - endpoint2) \
                 + (startpoint3 - endpoint3) > 0. and (dd1_ddt1 - dd1_ddt) + (dd2_ddt1 - dd2_ddt) \
                 + (dd3_ddt1 - dd3_ddt) < 0 and random.random() < 0.005: #  and slopechange < 0:
                '''The third instance of acceptance'''
                chi1 = chi2
                slopeolder1 = slope
                print 3, slope

            else:
                '''Rejection'''
                cont[j,1] -= change

            model3 = model2
            gc.collect()
    print P_v
    np.savetxt('CONTINUUM/NGC3783-continuum-slope-2-5-6',cont) #Updating the files
    np.savetxt('CONTINUUM/NGC3783-data_comp-slope-2-5-6',data_comp)
    #P_show = np.log10(P_v)
    plt.figure()
    plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    plt.scatter(cont[:,0],cont[:,1],color='g')
    plt.plot(time_model,model2)
    plt.ylim([4e-15,1.2e-14])
    plt.show(block=False)
    '''plt.figure()
    plt.scatter(freq,P_show,color='b')
    #plt.ylim([1e-21,1e-31])
    plt.xlim([1e-8,3e-7])
    #plt.yscale('log')
    plt.xscale('log')
    plt.show(block=False)'''
    gc.collect()

print cont[:,1]
