"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module images the stacked visibilities/ generates the stacked image
"""
#--------------------------------------------------------------------------------------------




import casa
import numpy as np
import json


#Loading the parameters json file:

with open('params.json') as json_data:
    d = json.load(json_data)
    
path = d["Simulator_params"]["working_directory"].encode("ascii","ignore") 

stampsize = d["Imaging_params"]["stampsize"]
cell_ = d["Imaging_params"]["cell"].encode("ascii","ignore")
vis_ = d["Imaging_params"]["stack_msfile_name"].encode("ascii","ignore")

# Image the uv-stacked data to produce an image.
casa.clean(vis=path+vis_, imagename=path+'uvstacked',cell = cell_,\
           imsize = stampsize, mask = [int(stampsize/2)-2, int(stampsize/2)-2,
              int(stampsize/2)+2, int(stampsize/2)+2])

casa.exportfits(path+'uvstacked.image', path+'uvstacked.fits') #exporting casa image to fts