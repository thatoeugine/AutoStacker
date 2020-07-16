import os
import json 
import sys
import subprocess
import configparser


#Loading the parameters json file:

with open('params.json') as json_data:
    d = json.load(json_data)
    
# working directory
PATH = d["Simulator_params"]["working_directory"]
os.system('mkdir %s'%PATH)
# casa container
CASA_CONTAINER = 'singularity exec /data/exp_soft/containers/casa-stable-5.6.2-2-2020-04-21.simg'
# kern3 container
KERN3_CONTAINER = 'singularity exec /data/exp_soft/containers/kern3.img'
# source finding container
SOURCE_FINDING_CONTAINER = 'singularity exec /data/exp_soft/containers/sourcefinding.img'
#stacker container
STACKER_CONTAINER_cuda = 'singularity exec /idia/software/containers/stacker-casa-divorce-cuda.simg'
STACKER_CONTAINER = 'singularity exec /data/exp_soft/containers/stacker-dev.simg'
# File names used here
uvstacked_ms = d["Imaging_params"]["stack_msfile_name"]
original_ms = d["Simulator_params"]["msfile_name"]

config = configparser.ConfigParser()
config.read(d["Simulator_params"]["trecs_params_file"])
skymodel= config.get('pipeline', 'project_name')+'_truthcat.txt'
skymodel_stacking_depth = d["Stacking_params"]["stacking_depth_skymodel_name"]
modelimage_ = config.get('field', 'fitsname')
modelimage = modelimage_.replace("-model.fits"," ")
trecs_params_file = d["Simulator_params"]["trecs_params_file"]

#======================================================================================================
#                                     Simulating with simms+meqtrees+Simuclass
#======================================================================================================

os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/simple_simulate/simulate.py')
os.system(CASA_CONTAINER+' '+'casa --nologger -c stack_auto/simple_simulate/noise_sefd.py')

if d["Stacking_params"]["No._of_srcs_of_choice?"] ==True:
    os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/stacking_depth_proto.py')
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+trecs_params_file)
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+PATH+skymodel)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)

elif d["Stacking_params"]["stacking_depth?"] ==True:
    os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/stacking_depth.py')
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+trecs_params_file)
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+PATH+skymodel)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)
    

else:
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+trecs_params_file)
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+skymodel_stacking_depth)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)


os.system(CASA_CONTAINER+' '+'casa --nologger -c stack_auto/Moving_vis_from_model.py')
          

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                       Stacking
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

os.system(STACKER_CONTAINER+' '+'casa --nologger -c stack_auto/Model_subtraction.py')


        #image stacking:
os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/im_stacking2.py')

        #vis stacking:
#os.system('CUDA_VISIBLE_DEVICES=0')
os.system(STACKER_CONTAINER_cuda+' '+'python stack_auto/vis_stacking.py')
os.system(CASA_CONTAINER+' '+'casa --nologger -c stack_auto/Image_stacked_data.py')

####
    #Predicting the visibilities to uv-stacked for comparison with Im-stacked
####
os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+'Imstacked'+' '+PATH+uvstacked_ms)

#Plotting results
#os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/plots_stacking_images.py')
#os.system(SOURCE_FINDING_CONTAINER+' '+'python stack_auto/plot_stacking_results.py')
