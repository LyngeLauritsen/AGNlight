# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:24:57 2017

@author: lynge
"""

import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
import aplpy as apl #The Astronomy Plotting Library for python
import astropy
import astropy.units as U
from astropy.coordinates import ICRS, Galactic, FK4, FK5, Angle, Latitude, Longitude
import astropy.constants as C
from astropy import wcs
import astropy.io.fits as fits
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import spectral_cube as SC
from spectral_cube import SpectralCube as sc
from astropy.wcs import WCS
#matplotlib.use('Agg')
import matplotlib.cm as cm
import astrometry as ast
import pyfits
import requests
import json
import os
import alipy
R = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": "sbfzysroepnryuht"})})
print R.text

quasar = pyfits.getdata('NGC3783_1_10_J_new.fits')
header = pyfits.getheader('NGC3783_1_10_J_new.fits')
w = WCS('NGC3783_1_10_J_new.fits')

#print header
#lon, lat = w.all_pix2world(309, 300, 0)
#print lon, lat
#print os.listdir('/../home/lynge/MasterP/')
#print quasar

#%%

def dist(mid,end):
    return np.sqrt((mid[0] - end[0])**2 + (mid[1] - end[1])**2)


#%%
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


#%%

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
    return AGNdata

apparature = 0.01
a = curve(pyfits.getdata('NGC3783_1_10_J_new.fits'),pyfits.getheader('NGC3783_1_10_J_new.fits'),WCS('NGC3783_1_10_J_new.fits'),(header['RA'],header['DEC']),apparature)

print a


plt.figure()
plt.imshow(a)
plt.colorbar()


#%%


app = []
flux = []


for k in range(20):
    apparature = 0.0005*(k+1)
    a = curve(pyfits.getdata('NGC3783_1_10_J_new.fits'),pyfits.getheader('NGC3783_1_10_J_new.fits'),WCS('NGC3783_1_10_J_new.fits'),(header['RA'],header['DEC']),apparature)
    #print np.shape(a),np.max(a[:,:,2])
    app.append(apparature)
    flux.append(np.sum(a))
    print apparature, a[1]
    

plt.figure()
plt.scatter(app,flux)
plt.xlabel('Apparature')
plt.ylabel('Flux')
plt.legend()
plt.show()

#pos = []
#val = []

#plt.figure()
#plt.imshow(a[:,:,2])
#plt.colorbar()

#%%

plt.figure()
plt.imshow(quasar)
#plt.axis([300,320,290,310])
plt.colorbar()
plt.show()

'''
plt.figure()
plt.imshow(a[:,:,2])
plt.axis([300,320,290,310])
plt.colorbar()
plt.show()
'''
'''
F_AGN[i+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[0]),j+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[1]),0] = w.all_pix2world(i+1,j+1,1)[0]
                F_AGN[i+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[0]),j+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[1]),1] = w.all_pix2world(i+1,j+1,1)[1]
                #F_AGN[i,j,2] = dis
                F_AGN[i+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[0]),j+int(apparature/pixspa-w.all_world2pix(header['RA'],header['DEC'],1)[1]),2] = quasar[i+1,j+1]
''' 
#%%
#print AGNdata.index(max(AGNdata))
print np.argmax(quasar[30:490,30:450])
d = 512 - 295
f = 512 - 305
g = 512 - 255
h = 512 - 265
#print d,f,g,h
#print w.all_pix2world(40,100)
#print w.all_pix2world(300,300)
#print w.all_pix2world(450,500)
#print w.all_pix2world(302,310)
#print AGNdata[295:310,305:315]
#for i in range(500):
#    print i, max(AGNdata[:,i])
#    print w.all_pix2world(i,300)
#print AGNdata[290:300,160:170]

#c = SkyCoord('header['RA'],header['DEC']') #SkyCoord('11h39m01.68s', '-37d44m18.2s')
c1 = SkyCoord('11h39m50s', '-37d44m18.2s')
#print c.ra.deg, c.dec.deg, -37.739, 174.757

#%%
#sc = SkyCoord(header['GEOLAT'],header['GEOLON'],header['GEOELEV'],frame='barycentrictrueecliptic',unit='deg',representation = 'unitspherical')
#print sc.galactic
#print sc.obsgeoloc

#astropy.coordinates.ICRS().representation_info

data1 = fits.open('NGC3783_1_10_J_new.fits')[0]
dataW = WCS(data1.header)
#sky = fits.open('NGC3783_3_10_K_sky.fits')
#correct = data - sky

'''
print data.header['CRVAL2']
print data.header['CDELT2']
print data.header['CRVAL1']
print data.header['CDELT1']
'''

ap = apl.FITSFigure(data1,figsize=(5,5))
#ap.recenter(c.ra.deg,c.dec.deg,0.0053391667)
ap.show_colorscale(aspect='auto')
ap.add_grid()
ap.show_circles(header['RA'],header['DEC'],0.0053391667)
ap.show_circles(float(w.all_pix2world(267,346,1)[0]),float(w.all_pix2world(267,346,1)[1]),0.01)

print header['RA'],header['DEC']
#ap.set_xaxis_coord_type('longitude')
#%%
print header['RA'],header['DEC']
ap = apl.FITSFigure(data1,figsize=(5,5))
ap.recenter(174.775,-37.755,0.0053391667)
ap.show_colorscale(aspect='auto')
ap.add_grid()
ap.show_circles(174.774,-37.754,0.004)

#%%
ap = apl.FITSFigure(data1,figsize=(5,5))
ap.recenter(header['RA'],header['DEC'],0.008)#174.774,-37.7540888889,0.004)
ap.show_colorscale(aspect='auto')
ap.add_grid()
ap.show_circles(header['RA'],header['DEC'],0.006)#174.774,-37.7540888889,0.004)
'''
ap = apl.FITSFigure(data,figsize=(5,5))
ap.recenter(c1.ra.deg, c1.dec.deg,0.004)
ap.show_colorscale(aspect='auto')
ap.add_grid()
ap.show_circles(c1.ra.deg, c1.dec.deg,0.004)
'''
