"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module performs image stacking of radio data

"""
#--------------------------------------------------------------------------------------------



import os
import numpy as np
import astropy.io.fits as fits
import json


#Loading the parameters json file:

with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"] 

#Image stacking with PASTA 
def fixcoords(coords):
    "Making coordinates suitable to be read by PASTA code"
    s = np.loadtxt(path+coords, usecols = (1,2,3))
    ra_ = s[:,0]
    dec_ = s[:,1]
    flux_ = s[:,2]

    outfile = open(path+'PASTAcoords.txt', 'w')
    for j in range(len(ra_)):
        srcline = '%.10f %.10f %.10f'%(ra_[j],dec_[j],flux_[j])
        print(srcline, file = outfile)
    outfile.close()


fixcoords(d["Simulator_params"]["skymodel_name"]+".txt")

os.system('genstack -l "0 1 2" /scratch/users/thatomanamela/Results/PASTAcoords.txt /scratch/users/thatomanamela/Results/2Dimage.fits 31')


### Editing the stacked image to be suitable for wsclean predict
data_ref, header_ref = fits.getdata(path+"2Dimage.fits", header = True)
data_GT_stacked, header_GT_stacked = fits.getdata(path+"PASTAcoords.txt2Dimage.fits_median.fits", header = True)

header_ref["SIMPLE"]  = header_GT_stacked["SIMPLE"]                     
header_ref["BITPIX"]  = header_GT_stacked["BITPIX"]                               
header_ref["NAXIS"]   = header_GT_stacked["NAXIS"]                   
header_ref["NAXIS1"]  = header_GT_stacked["NAXIS1"]                                                  
header_ref["NAXIS2"]  = header_GT_stacked["NAXIS2"]                                                                                                                 
header_ref["BUNIT"]   = header_GT_stacked["BUNIT"] 
header_ref["CTYPE1"]  = 'RA---SIN' 
header_ref["CTYPE2"]  = 'DEC--SIN'
header_ref["CDELT1"]  = header_GT_stacked["CDELT1"]
header_ref["CRVAL1"]  = header_GT_stacked["CRVAL1"]
header_ref["CRVAL2"]  = header_GT_stacked["CRVAL2"]



fits.writeto(path+'GroundT_stacked_header_ref-model.fits', data_GT_stacked, header_ref)
