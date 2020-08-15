import os
import sys
import configparser


#Loading the parameters file:
config = configparser.ConfigParser()
config.read(sys.argv[1])
    
# working directory
PATH = config.get('pipeline', 'data_path')
os.system('mkdir %s'%PATH)
# source finding container
SOURCE_FINDING_CONTAINER = 'singularity exec /data/exp_soft/containers/sourcefinding.img'
# kern3 container
KERN3_CONTAINER = 'singularity exec /data/exp_soft/containers/kern3.img'
# File names used here
uvstacked_ms = config.get('casa_imager', 'stack_msfile_name')  
original_ms = config.get('simulator_params', 'msfile_name')

skymodel= config.get('pipeline', 'project_name')+'_truthcat.txt'
skymodel_stacking_depth = config.get('stacking_params', 'stacking_depth_skymodel_name')
modelimage_ = config.get('field', 'fitsname')
modelimage = modelimage_.replace("-model.fits"," ")
#======================================================================================================
#                                     Simulating with simms+meqtrees+Simuclass
#======================================================================================================

if config.getboolean('pipeline', 'dosimulate'):
    from stack_auto.simple_simulate.Simulator import runSimulate, runAddNoise
    runSimulate(config, sys.argv[1])
    runAddNoise(config, sys.argv[1])


from stack_auto.Stacker import runFoVcut, runImstack, runImageStackData, runStack_clean, runMoveVis, runPlotStackingResults, runPlotStackingImages, runStackingDepth, runStackingDepthProto, runUVStack

if config.getboolean('stacking_params', 'No._of_srcs_of_choice') == True:
    runStackingDepthProto(config, sys.argv[1])
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+sys.argv[1])
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+PATH+skymodel)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)

if config.getboolean('stacking_params', 'stacking_depth') == True:
    runStackingDepth(config, sys.argv[1])
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+sys.argv[1])
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+PATH+skymodel)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)
    

else:
    os.system(KERN3_CONTAINER+' '+'python simuclass/simuCLASS.py'+' '+sys.argv[1])
    os.system(SOURCE_FINDING_CONTAINER+' '+'python simuclass/Edit_to_Tigger_frendly.py'+' '+skymodel_stacking_depth)
    
    ####
    #Predicting the visibilities
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+modelimage+' '+PATH+original_ms)


runMoveVis(config, sys.argv[1])
          

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                       Stacking
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if config.getboolean('pipeline', 'dostacking') == True:
    
    if config.getboolean('stacking_params', 'stack_clean') == True:
        runStack_clean(config, sys.argv[1])

    else:
        pass


            #image stacking:
    runImstack(config, sys.argv[1])

            #vis stacking:
    #os.system('CUDA_VISIBLE_DEVICES=0')

    runUVStack(config, sys.argv[1])

    runImageStackData(config, sys.argv[1])

    ####
        #Predicting the visibilities to uv-stacked for comparison with Im-stacked
    ####
    os.system(SOURCE_FINDING_CONTAINER+' '+'wsclean -mem 100 -predict -name'+' '+PATH+'Imstacked'+' '+PATH+uvstacked_ms)

    #Plotting results
    runPlotStackingImages(config, sys.argv[1])
    runPlotStackingResults(config, sys.argv[1])

else:
    'Simulation Complete'
