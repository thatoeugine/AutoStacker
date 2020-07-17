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
import stacker.uv
import numpy as np
import json



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                  Stacking by removing bright sources first.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Loading the parameters json file:

with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"].encode("ascii","ignore") 

def stack_clean(coords,FOV_size_cut,outputms):
    """Function stacks cleaned visibility file (i.e bright sources have been removed). 

    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Parameters:
    ----------
    vis: str # measurement set file name
    coords: str # stacking positions from previous observation, and its in a csv file
    """
    
    if(FOV_size_cut == True):
        Coords = stacker.readCoords(path+coords, unit='deg') # Reading stacking positiions
        stacker.uv.stack(Coords, path +'residual.ms', path+outputms,datacolumn='corrected',\
                         primarybeam= 'constant', use_cuda = True) # Visibility stacking
        
    else:
    # Stacking the Residual visibility file
        Coords = stacker.readCoords(path+coords, unit='deg') # Reading stacking positiions
        stacker.uv.stack(Coords, path+'residual.ms', path+outputms, datacolumn='corrected',\
                         primarybeam= 'constant', use_cuda = True) # Visibility stacking



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Normal stacking procedure without removing bright sources that act as noise in the observation.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def stack_Noclean(vis_,coords,FOV_size_cut,outputms):
    """Function stacking without removing bright sources. 

    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Parameters:
    ----------
    vis: str # measurement set file name    niter: int # number of cleaning iterations
    coords: str # stacking positions from previous observation, and its in a csv file
    """

    if(FOV_size_cut == True):
        Coords = stacker.readCoords(path+'fov_cut_coords.csv', unit='deg') # Reading stacking positiions
        stacker.uv.stack(Coords, path + vis_, path+ outputms,datacolumn='data',\
                         primarybeam= 'constant', use_cuda = True) # Visibility stacking

        
    else:
        Coords = stacker.readCoords(path+coords, unit='deg') # Reading stacking positiions
        stacker.uv.stack(Coords, path + vis_, path+ outputms, datacolumn='data',\
                         primarybeam= 'constant', use_cuda = True) # Visibility stacking


#==================
# Running tasks
#==================

print('visibility stacking ...........')


if d["Imaging_params"]["stack_clean?"] ==True:
    
    stack_clean('coords.csv', d["Stacking_params"]["FOV_size_cut?"],
                d["Imaging_params"]["stack_msfile_name"].encode("ascii","ignore"))
    
else:
    
    stack_Noclean(d["Simulator_params"]["msfile_name"].encode("ascii","ignore"),
                  'coords.csv',
                  d["Stacking_params"]["FOV_size_cut?"],
                  d["Imaging_params"]["stack_msfile_name"].encode("ascii","ignore"))