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
import json


######### loading json parameters######
with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"].encode("ascii","ignore") 

print "BEGINNING SIMULATION............................"

# creating empty MS table:
msoutname = path+d["Simulator_params"]["msfile_name"].encode("ascii","ignore")
ra0_deg = d["Simulator_params"]["ra_deg0"]
dec0_deg = d["Simulator_params"]["dec_deg0"]
obs_length = d["Simulator_params"]["obs_length_hr"]
dt = d["Simulator_params"]["dtime_s"]
freq0 = d["Simulator_params"]["freq0_MHz"].encode("ascii","ignore")
dfreq = d["Simulator_params"]["dfreq_MHz"].encode("ascii","ignore")
nchan = d["Simulator_params"]["nchan"]

simms_cmd = 'simms -T MeerKAT -t ascii -cs itrf -n %s -ra %sdeg -dec %sdeg -pl "XX XY YX YY" -st %i -dt %i -f0 %s -df %s -nc %i stack_auto/simple_simulate/meerkat.itrf.txt'%(msoutname, ra0_deg, dec0_deg,obs_length,dt,freq0,dfreq,nchan)
os.system(simms_cmd) # this creates an EMPTY Measurement Set


# now simulate an observation of a fake sky model you've created.
skymodelname = "stack_auto/simple_simulate/point_src.txt"
# The following will write the visiblities into the MODEL column of the empty MS you created. 
meqtrees_cmd = 'meqtree-pipeliner.py --mt 32 -c stack_auto/simple_simulate/tdlconf.profiles [turbo-sim] tiggerlsm.filename=%s ms_sel.msname=%s ms_sel.tile_size=1000000 ms_sel.output_column=DATA stack_auto/simple_simulate/turbo-sim.py =_simulate_MS'%(skymodelname,msoutname)
print(meqtrees_cmd)
os.system(meqtrees_cmd)
