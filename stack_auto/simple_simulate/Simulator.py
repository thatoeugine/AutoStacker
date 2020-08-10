import os

# source finding container
SOURCE_FINDING_CONTAINER = 'singularity exec /data/exp_soft/containers/sourcefinding.img'
# casa container
CASA_CONTAINER = 'singularity exec /data/exp_soft/containers/casa-stable-5.6.2-2-2020-04-21.simg'

def runSimulate(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} python {1}/stack_auto/simple_simulate/simulate.py {1}/{2}'.format(SOURCE_FINDING_CONTAINER,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
    
def runAddNoise(config, parameter_filename):
    
    run_dir = os.getcwd()

    cmd = '{0} casa --nologger --log2term -c {1}/stack_auto/simple_simulate/noise_sefd.py {1}/{2}'.format(CASA_CONTAINER,run_dir,
                                                           parameter_filename)
    print(cmd)
    os.system(cmd)
    
