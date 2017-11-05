import os
from astropy.io import fits

file_list = ['mscience0250.fits', 'mscience0270.fits', 'mscience0290.fits',
             'mscience0320.fits']
target_ra = [120.30651012, 120.31289098, 120.31954771, 120.32980469]
target_dec = [56.54986214, 56.55066694, 56.55146889, 56.55265498]

reference_frame_ra = target_ra[0]
reference_frame_dec = target_dec[0]

for i, filename in enumerate(file_list):

    hdu = fits.open(filename)
    jd = hdu[0].header['MIDTIMJD']

    print(hdu[0].header['CRVAL1'],hdu[0].header['CRVAL2'])
    
    ra_offset = target_ra[i] - reference_frame_ra
    dec_offset = target_dec[i] - reference_frame_dec

    print('shifting %s by %.1f/%.1f arcsec in RA/Dec' %
          (filename, ra_offset*3600, dec_offset*3600))
    
    hdu[0].header['CRVAL1'] = float(hdu[0].header['CRVAL1'])-ra_offset
    hdu[0].header['CRVAL2'] = float(hdu[0].header['CRVAL2'])-dec_offset
                               
    hdu.writeto(filename[:filename.find('.fits')]+'_shifted.fits',
                overwrite=True)
     
     
