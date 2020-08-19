"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module decides on the stacking depth, on the given telescopes resolution element and 
the number of sources being observed. 

"""
#--------------------------------------------------------------------------------------------

import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import configparser
import sys




def stacking_depth(cat,res_element, full_imagenoise):
    """
    This function decides on the stacking depth, on the given resolution element and 
    the number of sources being observed. 
    
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Parameters:
    ----------
    cat: #str, catalogue of sources being stacked
    res_element: #int, resolution element per source
    full_imagenoise: #float, the total image noise of the full image
    """
    cat_ = Table.read(cat, format = "ascii")
    Number_of_sources = len(cat_["peak_flux"])
    stacking_depth = Number_of_sources/res_element
    stacking_depth = int(stacking_depth)
    
    ## Making sure Not dominated by faint sources only, but including bright soirces too:
    above_stacked_noise_level = full_imagenoise/np.sqrt(stacking_depth) # select sources below the noise BUT: above the stacked noise level
    
    mask_faint = np.argwhere((5*above_stacked_noise_level<=cat_["integrated_flux"])&(cat_["integrated_flux"]<=full_imagenoise*5))
    #mask_faint = np.argwhere((cat_["integrated_flux"]<=10e-6))
    mask_faint = np.ravel(mask_faint) # flattening from 2D to 1D

    mask_bright = cat_["integrated_flux"].argsort()[-100:][::-1] # selecting the indexes of the brightest 100 srcs
    mask_bright = np.ravel(mask_bright) # flattening from 2D to 1D

    
    indices = []
    faint_sources_only  = []
    np.random.seed(1)
    for i in range(stacking_depth):
        value = np.random.choice(mask_faint)
        indices.append(value)
        faint_sources_only.append(value)
    
    
    for j in range(mask_bright.size):
        value2 = np.random.choice(mask_bright)
        indices.append(value2)
    
    
    ra = cat_["ra_abs"]
    dec = cat_["dec_abs"]
    

    outfile_csv = open(path+'coords.csv', 'w')
    for k in faint_sources_only:
        srcline_csv = '%.10f,%.10f'%(ra[k],dec[k])
        print>>outfile_csv,srcline_csv
    outfile_csv.close()
    
    
    
    # Making new model image consists of stacking depth sources plus all bright sources
    trecs_cat_ = Table.read(cat, format = "ascii")
    trecs_ra = trecs_cat_["ra_abs"]
    cat_ra = cat_["ra_abs"][indices]
    
    indices2 = []
    for l in range(len(cat_ra)):
        value3 = np.argwhere(cat_ra[l] == trecs_ra)
        indices2.append(value3)
        
    indices_true = []
    for e in range(len(indices2)):
        value4 = indices2[e][0][0]
        indices_true.append(value4)
        
    
    outfile = open('simuclass/simuCLASS/T_recs_catalogues/catalogue_SFGs_stacking_depth.txt', 'w')

    print>> outfile, '#lon            lat            size           flux            e1             e2             gamma1         gamma2'

    for m in indices_true:
        srcline ='%.10f           %.10f            %.10f           %.10f           %.10f       %.10f             %.10f             %.10f'%(trecs_cat_['ra_offset'][m],trecs_cat_['dec_offset'][m],trecs_cat_['size'][m],(trecs_cat_['integrated_flux'][m])*1e3,trecs_cat_['e1'][m],trecs_cat_['e2'][m],trecs_cat_['g1'][m],trecs_cat_['g2'][m])
        print>> outfile, srcline

    outfile.close()
    
    
    # catalogue of sizes and flux density for faint sources
    outfile_ = open(path+'faint_sourcesOnly.txt', 'w')

    print>> outfile_, '#size           flux'

    for t in faint_sources_only:
        srcline_ ='%.10f           %.10f'%(trecs_cat_['size'][t],trecs_cat_['integrated_flux'][t])
        print>> outfile_, srcline_

    outfile_.close()

    
#==================
# Running tasks
#==================



if __name__=='__main__':
    config = configparser.ConfigParser()
    config.read(sys.argv[-1])
    path = config.get('pipeline', 'data_path') 

    if config.getboolean('stacking_params', 'FOV_size_cut') == True: #cumputes FOV size cut for stacking across a certain FOV of interest

        FOV = config.getfloat('stacking_params', 'FOV_size_sqdeg') #[deg]
        diameter = config.getfloat('stacking_params', 'FOV_size_cut_value')*FOV #[deg]

        ra0_deg = config.getfloat('simulator_params', 'ra_deg0')
        dec0_deg = config.getfloat('simulator_params', 'dec_deg0')

        # Loading the RA and DEC from the T-RECS catalogue
        data_file = Table.read(config.get('stacking_params', 'stacking_depth_skymodel_name'), format = "ascii")
        RA,DEC = data_file["ra_abs"], data_file["dec_abs"]


        # Calculating the angular separtion between the pointing center and the sources:
        c1 = SkyCoord(ra0_deg*u.deg, dec0_deg*u.deg, frame='galactic')
        c2 = SkyCoord(RA*u.deg, DEC*u.deg,  frame='galactic')
        sep = c1.separation(c2)

        # Mask for a certian diameter of observation 
        mask = np.argwhere((sep.degree <= (diameter)))
        mask = np.ravel(mask)

        #New RA and DEC textfile:

        data = Table({'ra_offset': data_file["ra_offset"][mask],
                      'dec_offset': data_file["dec_offset"][mask],
                      'dec_abs': DEC[mask],
                       'ra_abs': RA[mask],
                      'integrated_flux': data_file["integrated_flux"][mask],
                      'size': data_file["size"][mask],
                      'peak_flux': data_file["peak_flux"][mask],
                     'e1': data_file["e1"][mask],
                     'e2': data_file["e2"][mask],
                     'g1': data_file["g1"][mask],
                     'g2': data_file["g2"][mask]},
                     names=['ra_offset','dec_offset','dec_abs','ra_abs','integrated_flux','size','peak_flux',\
                           'e1','e2','g1','g2'])

        ascii.write(data, path +'fov_cut_coords.txt', format='csv', fast_writer=False, overwrite=True) 

        # Run stacking depth function
        stacking_depth(path+'fov_cut_coords.txt',
                       config.getfloat('stacking_params', 'res_element_per_source'),
                       config.getfloat('stacking_params', 'im_noise_Jy'))



    else:
        stacking_depth(config.get('stacking_params', 'stacking_depth_skymodel_name'),
                       config.getfloat('stacking_params', 'res_element_per_source'),
                       config.getfloat('stacking_params', 'im_noise_Jy'))
