"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# Module cumputes FOV size cut

"""
#--------------------------------------------------------------------------------------------



import numpy as np
import json
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import configparser
import sys


if __name__=='__main__':
    config = configparser.ConfigParser()
    config.read(sys.argv[-1])

    path = config.get('pipeline', 'data_path') 


    FOV = config.getfloat('stacking_params', 'FOV_size_sqdeg') #[deg]
    diameter = config.getfloat('stacking_params', 'FOV_size_cut_value')*FOV #[deg]

    ra0_deg = config.getfloat('simulator_params', 'ra_deg0')
    dec0_deg = config.getfloat('simulator_params', 'dec_deg0')

    # Loading the RA and DEC from the T-RECS catalogue
    data_file = Table.read(path+ config.get('stacking_params', 'FOV_size_skyodel'),format = 'ascii')
    RA = data_file["ra_abs"]
    DEC = data_file["dec_abs"]


    # Calculating the angular separtion between the pointing center
    # and the sources:
    c1 = SkyCoord(ra0_deg*u.deg, dec0_deg*u.deg, frame='galactic')
    c2 = SkyCoord(RA*u.deg, DEC*u.deg,  frame='galactic')
    sep = c1.separation(c2)

    # Mask for a certian diameter of observation 
    mask = np.argwhere((sep.degree < (diameter)))

    #New RA and DEC textfile:
    outfile = open(path + "fov_cut_coords.csv", 'w')
    for i in range(len(RA[mask])):
        srcline = '%.10f,%.10f'%(RA[mask][i],DEC[mask][i])

        print>> outfile, srcline

    outfile.close()
