"""merge SCAMP output .head files into .fits header

put this script in the same directory as the data and run with 'python
merge_headers.py'

michael.mommert@nau.edu, Mar 21, 2017
"""

import os
import warnings
from astropy.io import fits

warnings.filterwarnings('ignore', category=fits.card.VerifyWarning)

for filename in os.listdir('.'):
    hdu = fits.open(filename, mode='update', verify='silentfix',
                    ignore_missing_end=True)

    headfilename = filename[:filename.find('.fit')]+'.head'
    try:
        newhead = open(headfilename, 'r').readlines()
    except:
        continue

    print('merging %s into %s' % (headfilename, filename))
    
    for line in newhead:
        key = line[:8].strip()
        try:
            value = float(line[10:30].replace('\'', ' ').strip())
        except ValueError:
            value = line[10:30].replace('\'', ' ').strip()
        comment = line[30:].strip()
        if key.find('END') > -1:
            break
        hdu[0].header[key] = (str(value), comment)

    hdu.flush()
