# README

Python module for simulating and stacking (uv and image) MeerKAT data. A modified version of [simuCLASS](https://bitbucket.org/itrharrison/simuclass/src/master/) is used in this module to simualate a radio interferometer sky model of the [TRECS](https://arxiv.org/abs/1805.05222) realist catalogue.

### Requirements:
ILIFU cluster


### To run:
```python
python Mastercode.py
```
**Note:**

The module consists of three config files: ***full_trecs.ini*** for generating a large TRECS sky model for a certain field of view of choice, ***stackingdepth2.ini*** for generating a small TRECS sky model based on the number of sources of choice and ***params.json*** for controlling the entire module. 1) An example of a fully documented ***.ini*** cofig file can be found [here](https://bitbucket.org/itrharrison/simuclass/src/master/example_verbose.ini)


```json

"Stacking_params":{
"stacking_depth?": false,
"res_element_per_source": 100,
"im_noise_Jy": 3e-6,
"stacking_depth_skymodel_name": "stack_auto/trecs-simu_truthcat.txt",
"FOV_size_cut?": false,
"FOV_size_sqdeg": 0.84,
"FOV_size_cut_value": 0.5,
"No._of_srcs_of_choice?": true,
"No._of_srcs": 10,
"flux_density_Jy": 1e-6,
"src_size_arcsec": 0 
}
```
