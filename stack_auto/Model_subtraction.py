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
import ConfigParser




if __name__=='__main__':
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[-1])

    print('CLEANing out bright sources ...........')

    path = config.get('pipeline', 'data_path') 
    vis_ = config.get('simulator_params', 'msfile_name')
    imagename_ = config.get('casa_imager', 'imagename')
    niter_ = config.getint('casa_imager', 'niter')
    threshold_ = config.get('casa_imager', 'threshold')
    imsize_ =  config.getint('casa_imager', 'imsize')
    cell_ = config.get('casa_imager', 'cell')
    weightings = config.get('casa_imager', 'weighting')
    rbst = config.getfloat('casa_imager', 'robust')

    # Creating full image of the observation
    
    casa.clean(path+vis_, path+imagename_, imagermode='csclean',\
               cell = cell_, imsize=imsize_, pbcor = True, niter=niter_, minpb=0.05,\
               threshold = threshold_, psfmode='hogbom', weighting = weightings, robust = rbst)
    
    
    # First produce component list of model
    stacker.modsub.cl_from_im(path+imagename_+'.model', path+'full.cl')
    # and subtract the component list from the data.
    stacker.modsub.modsub(path+'full.cl',path+vis_, path+'residual.ms',\
                          primarybeam='constant') # No pb since .model is not pbcorrected.
