#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
This module simulate TRECS sources with MeerKAT.

"""
#--------------------------------------------------------------------------------------------



import os
import configparser
import sys

if __name__=='__main__':
    config = configparser.ConfigParser()
    config.read(sys.argv[-1])
        
    path = config.get('pipeline', 'data_path') 
    
    print "BEGINNING SIMULATION............................"
    
    # creating empty MS table:
    msoutname = path+config.get('simulator_params', 'msfile_name')
    ra0_deg = config.getfloat('simulator_params', 'ra_deg0')
    dec0_deg = config.getfloat('simulator_params', 'dec_deg0')
    obs_length = config.getint('simulator_params', 'obs_length_hr')
    dt = config.getint('simulator_params', 'dtime_s')
    freq0 = config.get('simulator_params', 'freq0_MHz')
    dfreq = config.get('simulator_params', 'dfreq_MHz')
    nchan = config.getint('simulator_params', 'nchan')
    
    simms_cmd = 'simms -T MeerKAT -t ascii -cs itrf \
                       -n {0} \
                       -ra {1}deg \
                       -dec {2}deg \
                       -pl "XX XY YX YY" \
                       -st {3} \
                       -dt {4} \
                       -f0 {5} \
                       -df {6} \
                       -nc {7} stack_auto/simple_simulate/meerkat.itrf.txt'.format(msoutname,
                       ra0_deg, dec0_deg,obs_length,dt,freq0,dfreq,nchan)
    os.system(simms_cmd) # this creates an EMPTY Measurement Set
    
    
    # now simulate an observation of a fake sky model you've created.
    skymodelname = "stack_auto/simple_simulate/point_src.txt"
    # The following will write the visiblities into the MODEL column of the empty MS you created. 
    meqtrees_cmd = 'meqtree-pipeliner.py --mt 32 \-c stack_auto/simple_simulate/tdlconf.profiles [turbo-sim] tiggerlsm.filename=%s ms_sel.msname=%s ms_sel.tile_size=1000000 ms_sel.output_column=DATA stack_auto/simple_simulate/turbo-sim.py =_simulate_MS'%(skymodelname,msoutname)
    print(meqtrees_cmd)
    os.system(meqtrees_cmd)
    
print "SIMULATION COMPLETE,"


