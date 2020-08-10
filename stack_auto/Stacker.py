import os

# casa container
CASA_CONTAINER = 'singularity exec /data/exp_soft/containers/casa-stable-5.6.2-2-2020-04-21.simg'
# kern3 container
KERN3_CONTAINER = 'singularity exec /data/exp_soft/containers/kern3.img'
# source finding container
SOURCE_FINDING_CONTAINER = 'singularity exec /data/exp_soft/containers/sourcefinding.img'
#stacker container
STACKER_CONTAINER_cuda = 'singularity exec /idia/software/containers/stacker-casa-divorce-cuda.simg'
STACKER_CONTAINER = 'singularity exec /data/exp_soft/containers/stacker-dev.simg'

def runFoVcut(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/FOV_cut_coords.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
    
def runImstack(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/im_stacking2.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,parameter_filename)
    print(cmd)
    os.system(cmd)
    
    
def runImageStackData(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} casa --nologger --log2term -c {1}/stack_auto/Image_stacked_data.py {1}/{2}'.format(CASA_CONTAINER,run_dir, parameter_filename)
    print(cmd)
    os.system(cmd)
    

def runStack_clean(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} casa --nologger --log2term -c {1}/stack_auto/Model_subtraction.py {1}/{2}'.format(STACKER_CONTAINER,run_dir, parameter_filename)
    print(cmd)
    os.system(cmd)
    

def runMoveVis(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} casa --nologger --log2term -c {1}/stack_auto/Moving_vis_from_model.py {1}/{2}'.format(CASA_CONTAINER,run_dir,\
                                                                                                     parameter_filename)
    print(cmd)
    os.system(cmd)
    
    
def runPlotStackingResults(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/plot_stacking_results.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
    
def runPlotStackingImages(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/plots_stacking_images.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
    
def runStackingDepth(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/stacking_depth.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,parameter_filename)
    print(cmd)
    os.system(cmd)
    
def runStackingDepthProto(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/stacking_depth_proto.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,parameter_filename)
    print(cmd)
    os.system(cmd)
    
    
def runUVStack(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/vis_stacking.py {1}/{2}'.format(STACKER_CONTAINER_cuda,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
