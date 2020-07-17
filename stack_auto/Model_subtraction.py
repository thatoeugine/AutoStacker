"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module performs visibility stacking of radio data

"""
#--------------------------------------------------------------------------------------------

import stacker
import stacker.image
import stacker.uv
import stacker.modsub
import numpy as np
import json
import casa


#Loading the parameters json file:

with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"].encode("ascii","ignore") 



def stack_clean(vis_, imagename_, niter_, threshold_, imsize_, cell_, weightings, rbst):
    """Function removes bright sources for stacking Clean depth. 

    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Parameters:
    ----------
    vis_: str # measurement set file name
    imagename_: str # Full vis image name
    niter_: int # number of cleaning iterations
    threshold_: str # Flux level to stop cleaning, in units Jy
    imsize_: int # the size of the image 
    cell_: str # image pixel size in units of arcsecs
    weightings: str # clean weighting
    rbst: float # clean robustness
    """
    # Creating full image of the observation
    
    casa.clean(path+vis_, path+imagename_, imagermode='csclean',\
               cell = cell_, imsize=imsize_, pbcor = True, niter=niter_, minpb=0.05,\
               threshold = threshold_, psfmode='hogbom', weighting = weightings, robust = rbst)
    
    
    # First produce component list of model
    stacker.modsub.cl_from_im(path+imagename_+'.model', path+'full.cl')
    # and subtract the component list from the data.
    stacker.modsub.modsub(path+'full.cl',path+vis_, path+'residual.ms',\
                          primarybeam='constant') # No pb since .model is not pbcorrected.



#==================
# Running tasks
#==================

print('CLEANing bright sources ...........')

    
stack_clean(d["Simulator_params"]["msfile_name"].encode("ascii","ignore"),
            d["Imaging_params"]["imagename"].encode("ascii","ignore"),
            d["Imaging_params"]["niter"],
            d["Imaging_params"]["threshold"].encode("ascii","ignore"),
            d["Imaging_params"]["imsize"],
            d["Imaging_params"]["cell"].encode("ascii","ignore"),
            d["Imaging_params"]["weighting"].encode("ascii","ignore"),
            d["Imaging_params"]["robust"])
