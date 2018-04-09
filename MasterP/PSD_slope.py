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
OLC = np.loadtxt('NOVEMBER/Kelly-NGC3783K') # ''NOVEMBER/NOV-NGC3783-K')
OLC_error = np.loadtxt('NOVEMBER/NGC3783_NOISE_K.txt')

freqOLC = []

freqCLC = []

for i in range(int(len(OLC[:,0])/2.)):
    '''The frequency determination for the PSD'''
    freqOLC.append((i+1)/(60.*60.*24*float(OLC[:,0][len(OLC[:,0])-1] - OLC[:,0][0])))

for i in range(int(len(CLC[:,0])/2.)):
    '''The frequency determination for the PSD'''
    freqCLC.append((i+1)/(60.*60.*24*float(CLC[:,0][len(CLC[:,0])-1] - CLC[:,0][0])))

def PSD(freq,cont):
    '''Code designed to determine the Power Spectral Density Slope'''
    cont1 = cont[:,1] - np.nanmean(cont[:,1]) #Removing zero frequency

    cont_zero_freq = np.transpose(np.array([cont1,]*len(freq)))
    days_cont = np.transpose(np.array([cont[:,0],]*len(freq)))
    freq_array = np.array([freq,]*len(cont[:,0]))

    F_N_v = np.nansum(cont_zero_freq*np.cos(2*np.pi*freq_array*60.*60.*24*days_cont),axis=0)**2 \
    + np.nansum(cont_zero_freq*np.sin(2*np.pi*freq_array*60.*60.*24*days_cont),axis=0)**2

    #print len(cont[:,0]), F_N_v
    P_v = 2*60.*60.*24*float(cont[:,0][len(cont[:,0])-1] - cont[0,0]) \
    /(np.nanmean(cont[:,1])**2*len(cont[:,0])**2)*F_N_v #Finding PSD

    slope = np.polyfit(np.log10(freq),np.log10(P_v),1)[0] #Finding PSD slope

    return slope, P_v

slopeOLC, P_vOLC = PSD(freqOLC,OLC)

print slopeOLC

slopeCLC, P_vCLC = PSD(freqCLC,CLC)

print slopeCLC
