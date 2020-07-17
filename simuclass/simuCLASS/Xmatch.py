import numpy as np
import pyfits as pf
import sys
from astropy.io import fits


""" get Xmatched coords from COSMOS "stars", 
write out index in Laigle 2015 cat
"""

a = np.loadtxt('COSMOSstars.txt')
ras,decs = a[:,1:3].T

inputcat = 'COSMOS2015_Laigle+_v1.1.fits'
hdulist = fits.open(inputcat)
tbdata = hdulist[1].data
cols = hdulist[1].columns
cols.info('name')
rac = tbdata['ALPHA_J2000']
decc = tbdata['DELTA_J2000']


XmatchName = 'Xmatch-COSMOSstars.npy'

mindist = 1. # Xmatch distance in ARCSEC (min to be considered)


inds = []

for target in range(4): #len(ras)):
    print target
    sep = 3600 * np.sqrt(((ras[target] - rac)*(np.cos(decs[target])))**2 + (decs[target] - decc)**2)

    print sep

    n = np.where(sep == sep.min())[0][0]


    if (sep[n] < mindist):
        inds.append(n)
        print 'cross match found, sep = %.2f arcsec'%sep[n]

    elif (sep[n] > mindist):
        print 'nocounterpart'

inds = np.array(inds)

print('total cross-matches = %i'%len(inds))

np.save(XmatchName,inds)



