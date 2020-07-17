# README

Python module for simulating and stacking (uv and image) MeerKAT data. A modified version of [simuCLASS](https://bitbucket.org/itrharrison/simuclass/src/master/) is used in this module to simualate a radio interferometer sky model of the [TRECS](https://arxiv.org/abs/1805.05222) realist catalogue.

### Requirements:
ILIFU cluster


### To run:
```python
python Mastercode.py
```
**Note:**

1) The module consists of three config files: ***full_trecs.ini*** for generating a large TRECS sky model for a certain field of view of choice, ***stackingdepth2.ini*** for generating a small TRECS sky model based on the number of sources of choice and ***params.json*** for controlling the entire module. An example of a fully documented ***.ini*** cofig file can be found [here](https://bitbucket.org/itrharrison/simuclass/src/master/example_verbose.ini). For ***.json*** config file most of the parameters are selfexplainotry, the simulator_params should match with those of the .ini file.
2) for stacking paramters, 

```json
{"Simulator_params":{
"working_directory": "/scratch/users/thatomanamela/Results/testing/",
"trecs_params_file": "stackingdepth2.ini",  # .ini being used 

"Stacking_params":{
"stacking_depth?": false,  # should be kept false when "No._of_srcs_of_choice?" is true
"res_element_per_source": 100,  # the number of resolution elements per source
"im_noise_Jy": 3e-6,      # MeerKAT's Full image noise 
"stacking_depth_skymodel_name": "stack_auto/trecs-simu_truthcat.txt",
"FOV_size_cut?": false,
"FOV_size_sqdeg": 0.84,   # MeerKAT's FOV size in sq-deg
"FOV_size_cut_value": 0.5,  # sources will be generated within 50% of the FOV if "FOV_size_cut?" is true
"No._of_srcs_of_choice?": true,  # when resolutions elements per source is not prefered
"No._of_srcs": 10,    # flux density of choices provided "No._of_srcs_of_choice?" is true
"flux_density_Jy": 1e-6, # flux density of choices provided "No._of_srcs_of_choice?" is true
"src_size_arcsec": 0 # when 0 TRECS source sizes are used
}
```
