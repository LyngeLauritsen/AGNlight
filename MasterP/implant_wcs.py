"""implant wcs information into image header

adjust values below accordingly and the script will implant this
information into all FITS files in the same directory

michael.mommert@nau.edu, Mar 21, 2017
"""

import os
from astropy.io import fits

for filename in os.listdir('.'):
    if not '.new' in filename:
        continue

    print(filename)
    
    hdu = fits.open(filename, mode='update')
    header = hdu[0].header

    header['CRTYPE1'] = 'RA---TAN'
    header['CRTYPE2'] = 'DEC--TAN'    
    header['CRPIX1'] = header['_RPIX1'] # use center of field x axis
    header['CRPIX2'] = header['_RPIX2'] # use center of field y axis
    header['CRVAL1'] = header['_RVAL1'] # use RA of center of field
    header['CRVAL2'] = header['_RVAL2'] # use Dec of center of field
    header['CD1_1'] = header['PIXSCALE']/3600 # pixel scale (arcsec/px) / 3600
    header['CD1_2'] = 0 # 0 if no field rotation
    header['CD2_1'] = 0 # 0 if no field rotation   
    header['CD2_2'] = header['PIXSCALE']/3600. # pixel scale (arcsec/px) / 3600

    hdu.flush()

