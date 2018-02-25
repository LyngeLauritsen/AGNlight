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


runs = 15

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

data = np.loadtxt('NOVEMBER/Kelly-NGC3783K') #'NOVEMBER/NOV-NGC3783-K')
error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')

print data[:,0]
print data[:,1]

print error[:,1]
step = 20.
cont_days = np.arange(min(data[:,0])-200.,max(data[:,0]),step)
cont = np.zeros((len(cont_days),3))
cont[:,0] = cont_days
cont[:,1] = np.nanmean(data[:,1])
day = cont_days[0]

data_comp = np.zeros((len(data[:,1]),3))
data_comp[:,0] = data[:,0]
data_comp[:,1] = data[:,1]
data_comp[:,2] = 0
#print data[:,2]

flux1 = cont[:,1]
time1 = cont[:,0]
sigma1 = cont[:,2] + np.nanmean(error[:,1])
time = np.insert(time1,0,0)
sigma = np.insert(sigma1,0,np.mean(sigma1))
flux = np.insert(flux1,0,np.mean(flux1))

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
    res = mini(prob,param,method='Nelder-Mead',tol=1e-18)
    tau = res.x[1] #1400. #res.x[1] #1400. #res.x[1]
    b = res.x[0] #9.2/tau #res.x[0] #7.5*10**(-15)/tau #res.x[0]
    sigma_tot = res.x[2] #0.01 #res.x[2] #1.3*10**(-16) #res.x[2] 
    model = np.zeros((50,len(flux1)))
    flux_test = flux1[0]
    dx_day = np.zeros((50,len(flux1)))
    for j in range(len(model[:,0])):
        model[j,0] = flux1[0]
        for i in range(len(flux1)-1):
            dt = (time1[i+1]-time1[i])
            epsilon = np.random.normal(0,1,1)
            dX1 = (flux1[i+1])#*(-1) + model[j,0]
            change = dX(tau,sigma_tot,dt,b,epsilon,dX1) #flux1[i+1] - flux1[i]
            model[j,i+1] = change + flux1[i] #model[j,i] # - flux_test)
            dx_day[j,i+1] = change/dt
    model2 = []
    time_model = []
    std = []
    dx_days = []
    for i in range(len(model[0,:])):
        model2.append(np.mean(model[:,i]))
        std.append(np.std(model[:,i]))
        time_model.append(time1[i])
        dx_days.append(np.mean(dx_day[:,i]))
    flux_model = np.array((time_model,model2,std,res.x))
    return flux_model,time_model, model2, dx_days

model3 = cont[:,1]
tau = 400.
sigma_tot = np.nanmean(error[:,1])
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
            startpoint = abs(cont[j,1]-model3[j])/cont[j,1]
            change = (-1)**randint(0,2)*cont[j,1]*random.random()*max_jump
            cont[j,1] += change
            h = 0
            b = np.nanmean(cont[:,1])
            flux1 = cont[:,1]
            param = [b,tau,sigma_tot]
            model, time_model, model2, dx_days = flux_model(param)
            for h in range(len(data_comp[:,0])):
                for k in range(len(cont_days)-1):
                    if cont_days[k] < data_comp[h,0]:
                        data_comp[h,2] += step*np.mean((cont[k,1],cont[k+1,1]))*log[abs(int((cont[k,0]-data_comp[h,0])))]
            #print data_comp[:,2]
            chi2 = np.nansum((data_comp[:,1] - data_comp[:,2])**2)
            chi_test = random.random()
            #jump = abs((cont[j-2,1]+cont[j-1,1]+cont[j+1,1]+cont[j+2,1])/4. - cont[j,1])
            #if chi2 <= chi1 and abs((cont[j,1]-cont[j-1,1])/cont[j,1])/(float(cont[j,0]-cont[j-1,0])) < 0.006 \
            #   and abs((cont[j+1,1]-cont[j,1])/cont[j,1])/float(cont[j,0]-cont[j-1,0]) < 0.006 or chi_test < 0.05:
            if chi2 <= chi1 and abs((cont[j,1]-model2[j])/cont[j,1]) < 0.02:
                chi1 = chi2
            elif chi2 <= chi1 and (startpoint - abs((cont[j,1]-model2[j])/cont[j,1])) < 0.:
                chi1 = chi2
            else:
                cont[j,1] -= change
            model3 = model2
            gc.collect()
    plt.figure()
    plt.scatter(data_comp[:,0],data_comp[:,1],color='b')
    plt.scatter(data_comp[:,0],data_comp[:,2],color='r')
    plt.scatter(cont[:,0],cont[:,1],color='g')
    plt.plot(time_model,model2)
    plt.ylim([4e-15,1.2e-14])
    plt.show(block=False)
    gc.collect()

print cont[:,1]

np.savetxt('continuum3',cont)
np.savetxt('data_comp3',data_comp)
