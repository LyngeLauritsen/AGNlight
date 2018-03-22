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


runs = 3000
runs1 = 10
step = 10. #Timestep (days)
slopeaim = -3. #Slope that is targeted
slope = 0
mean = 140. #Transfer function
sigma = 0.8
mu = np.log(mean) - 1/2.*sigma

def lognorm(x,mu,sigma):
    '''Defines the transfer function'''
    sigma = float(sigma)
    mu = float(mu)
    x = float(x)
    exp = -((np.log(x)-mu)**2/(2*sigma**2))
    front = 1/(x*sigma*np.sqrt(2*np.pi))
    return front*np.exp(exp)


x_list = np.linspace(0.01,1500,1500)



log = []
for i in range(len(x_list)):
    '''Creates the list used to identify the degree of transfered light'''
    log.append(lognorm(x_list[i],mu,sigma))
    
print np.nansum(log)*1500/5000.

data = np.loadtxt('NOVEMBER/NOV-NGC3783-K') #'NOVEMBER/Kelly-NGC3783K') 
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')

dd_ddt_data = []
for j in range(len(data[:,0])-2):
    '''Purely a test to see the largest double derivative of the observed light curve'''
    dd_ddt_data.append(abs(((data[j+1,1] - data[j,1])/(data[j+1,0] - data[j,0]) - (data[j,1] - data[j-1,1])/(data[j,0] - data[j-1,0]))/(data[j+1,0] - data[j-1,0])/2.))

print np.max(dd_ddt_data)
dd_ddt_allowed = np.nanmax(dd_ddt_data)*1.4 #Identifies the maximal double derivative allowed in the continuum light curve
print data[:,0]
print data[:,1]

print error[:,1]

cont_days = np.arange(min(data[:,0])-400.,max(data[:,0]+100),step) #Defines the spacing of the light curve
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = np.nanmean(data[:,1])
#cont = np.loadtxt('CONTINUUM/NGC3783-continuum-slope-2-5-3point')
day = cont_days[0]

data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0
#data_comp = np.loadtxt('CONTINUUM/NGC3783-data_comp-slope-2-5-2')
#print data[:,2]

flux1 = cont[:,1]
time1 = cont[:,0]
sigma1 = cont[:,2] + np.nanmean(error[:,1])
time = np.insert(time1,0,0)
sigma = np.insert(sigma1,0,np.mean(sigma1))
flux = np.insert(flux1,0,np.mean(flux1))

'''The set of equations that forms the basis for the Kelly function'''

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
        a.append(a_i(time[i+1],time[i],param[1]))
        x_hat.append(x_hat_i(a[i],x_hat[i],omega[i],sigma[i],x_star[i]))
        omega.append(omega_i(omega_0(param[1],param[2]),a[i],omega[i],sigma[i]))
        x_star.append(x_star_i(flux[i+1],param[0],param[1]))
    x_hat[0] = np.mean(x_hat)
    x_hat[1] = np.mean(x_hat)
    return x_hat,omega,x_star,a

def prob(param):
    x_hat,omega,x_star,a = parameters(param)
    probability = 1
    for i in range(len(a)):
        part = np.log(1/np.sqrt(2*np.pi*(omega[i+1] + sigma[i+1]**2)))\
        -(1/2.)*(((x_hat[i+1] - x_star[i+1])**2.)/(omega[i+1] + sigma[i+1]**2.))
        probability += part
    return -probability

def dX(tau,sigma,dt,b,epsilon,X):
    dt = dt
    sigma = sigma
    return (-(1/tau)*X*dt + sigma*np.sqrt(dt)*epsilon + b*dt)

def flux_model(param):
    '''Creating the kelly function'''
    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
    b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
    sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
    model = np.zeros((20,len(flux1)))
    flux_test = flux1[0]
    dx_day = np.zeros((20,len(flux1)))
    for j in range(len(model[:,0])):
        model[j,0] = flux1[0]
        for i in range(len(flux1)-1):
            dt = (time1[i+1]-time1[i])
            epsilon = np.random.normal(0,1,1)
            dX1 = (flux1[i+1])#*(-1) + model[j,0]
            change = dX(tau,sigma_tot,dt,b,epsilon,dX1) #flux1[i+1] - flux1[i]
            model[j,i+1] = change + flux1[i] #model[j,i] # - flux_test)
            dx_day[j,i+1] = change/dt
            #print 1
    model2 = []
    time_model = []
    std = []
    dx_days = []
    for i in range(len(model[0,:])):
        model2.append(np.mean(model[:,i]))
        std.append(np.std(model[:,i]))
        time_model.append(time1[i])
        dx_days.append(np.mean(dx_day[:,i]))
        #print 1
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model,time_model, model2, dx_days

def flux_model1(param):
    '''Forming the Kelly function from right to left'''
    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
    b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
    sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
    model = np.zeros((20,len(flux1)))
    flux_test = flux1[0]
    dx_day = np.zeros((20,len(flux1)))
    timef = time1[::-1]
    fluxf = flux1[::-1]
    #print timef, time1
    for j in range(len(model[:,0])):
        model[j,0] = fluxf[0]
        for i in range(len(fluxf)-1):
            dt = abs(timef[i+1]-timef[i])
            epsilon = np.random.normal(0,1,1)
            dX1 = (fluxf[i+1])#*(-1) + model[j,0]
            change = dX(tau,sigma_tot,dt,b,epsilon,dX1) #flux1[i+1] - flux1[i]
            model[j,i+1] = change + fluxf[i] #model[j,i] # - flux_test)
            dx_day[j,i+1] = change/dt
            #print 2
    model2 = []
    time_model = []
    std = []
    dx_days = []
    for i in range(len(model[0,:])):
        model2.append(np.mean(model[:,i]))
        std.append(np.std(model[:,i]))
        time_model.append(timef[i])
        dx_days.append(np.mean(dx_day[:,i]))
        #print 2
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model,time_model[::-1], model2[::-1], dx_days[::-1]

#freq = np.empty((np.shape(time1)[0],), dtype=np.object_)
#freq.fill([])
#freq = np.frompyfunc(list,1,1)(freq)'

freq = []

for i in range(int(len(cont_days)*0.5)):
    freq.append((i+1)/(60.*60.*24*float(cont_days[len(cont_days)-1] - cont_days[0])))

print freq
print np.shape(freq)
print len(cont_days)/(60.*60.*24*float(2)*float(cont_days[len(cont_days)-1] - cont_days[0]))

model3 = cont[:,1]
tau = 800.
sigma_tot = np.nanmean(error[:,1])
chi1 = 1e10000
max_jump = 0.01

slopeolder1 = -1e1000

F_N_v = np.zeros((len(freq)))
P_v = 0
print F_N_v

for i in range(runs):
    #if i == int(0.4*runs):
    #    max_jump = 0.01
    #if i == int(0.8*runs):
    #    max_jump = 0.005
    for i1 in range(runs1):
        print i, i1/float(runs1), slopeolder1
        for j in range(len(cont_days)-5):
            #print i, j/float(len(co nt_days))
            F_N_v = np.zeros((len(freq)))
            data_comp[:,2] = 0
            dd1_ddt = abs(((cont[j+1,1] - cont[j,1])/step - (cont[j,1] - cont[j-1,1])/step)/step/2.) 
            dd2_ddt = abs(((cont[j,1] - cont[j-1,1])/step - (cont[j-1,1] - cont[j-2,1])/step)/step/2.) 
            dd3_ddt = abs(((cont[j+2,1] - cont[j+1,1])/step - (cont[j+1,1] - cont[j,1])/step)/step/2.)
            dd4_ddt = abs(((cont[j+3,1] - cont[j+2,1])/step - (cont[j+2,1] - cont[j+1,1])/step)/step/2.)
            dd5_ddt = abs(((cont[j+4,1] - cont[j+3,1])/step - (cont[j+3,1] - cont[j+2,1])/step)/step/2.)
            startpoint1 = abs(cont[j-1,1]-model3[j-1])/cont[j-1,1]/step
            startpoint2 = abs(cont[j,1]-model3[j])/cont[j,1]/step
            startpoint3 = abs(cont[j+1,1]-model3[j+1])/cont[j+1,1]/step
            startpoint4 = abs(cont[j+2,1]-model3[j+2])/cont[j+2,1]/step
            startpoint5 = abs(cont[j+3,1]-model3[j+3])/cont[j+3,1]/step
            direction = randint(0,2)
            change1 = (-1)**direction*cont[j,1]*random.random()*max_jump
            change2 = (-1)**direction*cont[j+1,1]*random.random()*max_jump
            change3 = (-1)**direction*cont[j+2,1]*random.random()*max_jump
            cont[j,1] += change1
            cont[j+1,1] += change2
            cont[j+2,1] += change3
            h = 0
            b = np.nanmean(cont[:,1])
            flux1 = cont[:,1]
            param = [b,tau,sigma_tot]
            #p1 = Process(target=flux_model(param))
            #p2 = Process(target=flux_model1(param))
            #p1.start()
            #p2.start()
            #model, time_model, model2, dx_days = p1.join()
            #model1, time_model1, model12, dx_days1 = p2.join()
            model, time_model, model2, dx_days = flux_model(param) #p1.join()
            model1, time_model1, model12, dx_days1 = flux_model1(param) #p2.join()
            #model = np.add(model,model1)/2.
            time_model = np.add(time_model,time_model1)/2.
            model2 = np.add(model2,model12)/2.
            cont1 = cont[:,1] - np.nanmean(cont[:,1])
            model21 = model2[:] - np.nanmean(model2[:])
            for i1 in range(len(freq)):
                F_N_v[i1] = np.sum(cont1[:]*np.cos(2*np.pi*freq[i1]*60.*60.*24*cont[:,0]))**2 \
                             + np.sum(cont1[:]*np.sin(2*np.pi*freq[i1]*60.*60.*24*cont[:,0]))**2
            P_v = 2*60.*60.*24*float(cont_days[len(cont_days)-1] - cont_days[0]) \
                  /(np.nanmean(cont[:,1])**2*len(cont_days)**2)*F_N_v**2
            #print P_v
            #print 1/np.nanmean(cont[:,1])**2
            slope = np.polyfit(np.log10(freq),np.log10(P_v),1)[0]
            for h in range(len(data_comp[:,0])):
                for k in range(len(cont_days)-1):
                    if cont_days[k] < data_comp[h,0]:
                        data_comp[h,2] += step*np.mean((cont[k,1],cont[k+1,1]))*log[abs(int((cont[k,0]-data_comp[h,0])))]
            #print data_comp[:,2]
            chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2)
            chi_test = random.random()
            endpoint1 = abs(cont[j-1,1]-model2[j-1])/cont[j-1,1]/step
            endpoint2 = abs(cont[j,1]-model2[j])/cont[j,1]/step
            endpoint3 = abs(cont[j+1,1]-model2[j+1])/cont[j+1,1]/step
            endpoint4 = abs(cont[j+2,1]-model2[j+2])/cont[j+2,1]/step
            endpoint5 = abs(cont[j+3,1]-model2[j+3])/cont[j+3,1]/step
            dd1_ddt1 = abs(((cont[j+1,1] - cont[j,1])/step - (cont[j,1] - cont[j-1,1])/step)/step/2.) 
            dd2_ddt1 = abs(((cont[j,1] - cont[j-1,1])/step - (cont[j-1,1] - cont[j-2,1])/step)/step/2.) 
            dd3_ddt1 = abs(((cont[j+2,1] - cont[j+1,1])/step - (cont[j+1,1] - cont[j,1])/step)/step/2.)
            dd4_ddt1 = abs(((cont[j+3,1] - cont[j+2,1])/step - (cont[j+2,1] - cont[j+1,1])/step)/step/2.)
            dd5_ddt1 = abs(((cont[j+4,1] - cont[j+3,1])/step - (cont[j+3,1] - cont[j+2,1])/step)/step/2.)
            #print dd_ddt1
            #jump = abs((cont[j-2,1]+cont[j-1,1]+cont[j+1,1]+cont[j+2,1])/4. - cont[j,1])
            #if chi2 <= chi1 and abs((cont[j,1]-cont[j-1,1])/cont[j,1])/(float(cont[j,0]-cont[j-1,0])) < 0.006 \
            #   and abs((cont[j+1,1]-cont[j,1])/cont[j,1])/float(cont[j,0]-cont[j-1,0]) < 0.006 or chi_test < 0.05:
            slopechange = (slope - slopeaim)**2 - (slopeolder1 - slopeaim)**2
            #print slope
            if chi2 <= chi1 and slopechange < 0. and dd1_ddt1 < dd_ddt_allowed and dd2_ddt1 < dd_ddt_allowed \
               and dd3_ddt1 < dd_ddt_allowed and dd4_ddt1 < dd_ddt_allowed and dd5_ddt1 < dd_ddt_allowed: #\ #and endpoint1 < 0.0055/2. and endpoint2 < 0.0055/2. and endpoint3 < 0.0055/2.
               # and slopechange < 0.:
                chi1 = chi2
                slopeolder1 = slope
                print 1, slope, dd1_ddt1
            elif chi2 <= chi1 and slopeaim - 0.05 < slope < slopeaim + 0.05 and dd1_ddt1 < dd_ddt_allowed \
                 and dd2_ddt1 < dd_ddt_allowed and dd3_ddt1 < dd_ddt_allowed and dd4_ddt1 < dd_ddt_allowed \
                 and dd5_ddt1 < dd_ddt_allowed: #and endpoint1 < 0.0055/2. and endpoint2 < 0.0055/2. and endpoint3 < 0.0055/2. \
               #and dd_ddt1 < 3.e-17 and slopeaim - 0.05 < slope < slopeaim + 0.05:
                chi1 = chi2
                slopeolder1 = slope
                print 2, slope
            elif chi2 <= chi1 and (startpoint1 - endpoint1) + (startpoint2 - endpoint2) \
                 + (startpoint3 - endpoint3) + (startpoint4 - endpoint4) + (startpoint5 - endpoint5) > 0. \
                 and (dd1_ddt1 - dd1_ddt) + (dd2_ddt1 - dd2_ddt) \
                 + (dd3_ddt1 - dd3_ddt) + (dd4_ddt1 - dd4_ddt) + (dd5_ddt1 - dd5_ddt) < 0 \
                 and random.random() < 0.05 and slopeaim - 0.3 < slope < slopeaim + 0.3: #  and slopechange < 0:
                chi1 = chi2
                slopeolder1 = slope
                print 3, slope
            #elif chi2 <= chi1 and (startpoint1 - endpoint1) + (startpoint2 - endpoint2) \
            #     + (startpoint3 - endpoint3) > 0.: # and dd_ddt1 < dd_ddt and slopeaim - 0.05 < slope < slopeaim + 0.05:
            #    chi1 = chi2
            #    slopeolder1 = slope
            #    print 4, slope
            else:
                cont[j,1] -= change1
                cont[j+1,1] -= change2
                cont[j+2,1] -= change3
            model3 = model2
            gc.collect()
    print P_v
    np.savetxt('CONTINUUM/NGC3783-continuum-slope-3-3point',cont)
    np.savetxt('CONTINUUM/NGC3783-data_comp-slope-3-3point',data_comp)
    P_show = np.log10(P_v)
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


