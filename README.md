# README

Python module for simulating and stacking (uv and image) MeerKAT data. A modified version of [simuCLASS](https://bitbucket.org/itrharrison/simuclass/src/master/) is used in this module to simulate a radio interferometer sky model of the [TRECS](https://arxiv.org/abs/1805.05222) realist catalogue.

---

### Requirements:
[ILIFU](http://docs.ilifu.ac.za/#/) cluster

---

### To run:
```python
python Mastercode.py
```
---

**Note:**

The module consists of three config files: ***full_trecs.ini*** for generating a large TRECS sky model for a certain field of view of choice, ***stackingdepth2.ini*** for generating a small TRECS sky model based on the stacking depth and ***params.json*** for controlling the entire module. An example of a fully documented ***.ini*** config file can be found [here](https://bitbucket.org/itrharrison/simuclass/src/master/example_verbose.ini) and the ***params.json*** config file documentation can be found [here](https://github.com/thatoeugine/uv-stacking_MeerKAT_sim_data.wiki.git), in a case where the default parameters need to be changed.
