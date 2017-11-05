# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 14:17:11 2017

@author: Lynge R. B. Lauritsen
"""

import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
#import aplpy as apl #The Astronomy Plotting Library for python
import astropy
import astropy.units as U
from astropy.coordinates import ICRS, Galactic, FK4, FK5, Angle, Latitude, Longitude
import astropy.constants as C
from astropy import wcs
import astropy.io.fits as fits
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
#import spectral_cube as SC
#from spectral_cube import SpectralCube as sc
from astropy.wcs import WCS
#matplotlib.use('Agg')
import matplotlib.cm as cm
#import astrometry as ast
import pyfits
import requests
import json
import os
#from sklearn import datasets, linear_model
from scipy.optimize import curve_fit
#import alipy

quasar = pyfits.getdata('NGC3783_1_10_J_new.fits')
header = pyfits.getheader('NGC3783_1_10_J_new.fits')
w = WCS('NGC3783_1_10_J_new.fits')

def dist(mid,end):
    return np.sqrt((mid[0] - end[0])**2 + (mid[1] - end[1])**2)


def Basic(quasar,header):
    '''Provides a basic list of neccesities from the fits header file for use in the later flux determination, \
    it does not serve an individual purpose beyound that'''
    AGNdata = quasar #np.rot90(np.rot90(np.fliplr(quasar)))
    #header = data[0].header
    #print header
    RAstd = header['CRVAL1']
    DECstd = header['CRVAL2']
    #print RAstd, DECstd
    pixspa = header['PIXSCALE']/(3600.)
    #print pixspa
    RAstdpix = header['CRPIX1']
    DECstdpix = header['CRPIX2']
    exptime = header['EXPTIME']
    return header,AGNdata,RAstd,DECstd,pixspa,RAstdpix,DECstdpix,exptime


def curve(quasar,header,w,center,apparature):
    '''Determines the FLUX of the stellar object given the the numpy array of the image, the header of the fits file, the astropy coordinate representation, \
    the object position and the apparature of interest'''
    header,AGNdata,RAstd,DECstd,pixspa,RAstdpix,DECstdpix,exptime = Basic(quasar,header)
    #AGNdata = np.swapaxes(AGNdata,0,1)
    y, x = np.ogrid[0:header['NAXIS1'],0:header['NAXIS1']]
    #w = np.swapaxes(w,0,1)
    #F_AGN = np.zeros((512,512,3)) #((int(2*apparature/pixspa)+1,int(2*apparature/pixspa)+1,3))
    #print center
    x1 = w.all_world2pix(center[0],center[1],0)
    x1,y1 = (x1[0]),(x1[1])
    mask = ((y-y1)**2 + (x-x1)**2) > (apparature/float(pixspa))**2
    AGNdata[mask] = float(0)
    mask = AGNdata == 0
    rows = np.flatnonzero((~mask).sum(axis=1))
    cols = np.flatnonzero((~mask).sum(axis=0))
    AGNdata = AGNdata[rows.min():rows.max()+1, cols.min():cols.max()+1]
    return AGNdata

apparature = 0.006
a = curve(pyfits.getdata('NGC3783_1_10_J_new.fits'),pyfits.getheader('NGC3783_1_10_J_new.fits'),WCS('NGC3783_1_10_J_new.fits'),(header['RA'],header['DEC']),apparature)

print a


plt.figure(figsize=(15,15))
plt.imshow(a,cmap='hot')
plt.colorbar()
plt.show()


def func(x,a,b):
    return a*x + b


app = []
flux = []


for k in range(15):
    apparature = 0.0005*(k+1)
    a = curve(pyfits.getdata('NGC3783_1_10_J_new.fits'),pyfits.getheader('NGC3783_1_10_J_new.fits'),WCS('NGC3783_1_10_J_new.fits'),(header['RA'],header['DEC']),apparature)
    #print np.shape(a),np.max(a[:,:,2])
    app.append(apparature)
    flux.append(np.sum(a))
    print apparature, np.sum(a)
    #plt.figure(figsize=(10,10))
    #plt.imshow(a)
    #plt.colorbar()
    #plt.show()

x = np.linspace(0,0.01,10)
y = func(x,float(1000000),float(0))

popt, pcov = curve_fit(func, x, y)

print popt, pcov

plt.figure()
plt.scatter(app,flux)
plt.plot(x,y)
plt.xlabel('Apparature')
plt.ylabel('Flux')
plt.xscale('log')
#plt.yscale('log')
plt.xlim([0.0005,0.01])
plt.legend()
plt.show()
