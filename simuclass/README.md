# README #

About:

SimuCLASS is a modular pipeline for simulating radio interferometer observations,
from sky model to deconvolved images, built in support of the SuperCLASS experiment:
http://www.e-merlin.ac.uk/legacy/projects/superclass.html

Not really intended for public consumption yet, but if you want to give it a go you
are welcome!

Prerequisites:

* Astropy
* CASA
* GalSim
* WSCLEAN (optional)
* AIPS (optional)

Usage:

python simuCLASS.py example.ini

Components:

* exportdata, importdata : code for import/export of data between CASA and AIPS
* skymodel : code to create a simulated sky model using GalSim
* simulatedata : code to create simulated visibility data from a skymodel and telescope configuration
* imager : code to run image deconvloution, currently CASA, AIPS and WSCLEAN
* thumbnailer : code to generate thumbnails of individual simulated sources (e.g. for shape measurement)

Contact:

Ian Harrison, ian 'dot' harrison 'dash' 2 'at' manchester 'dot' ac 'dot' uk

To Do:

* Add calibration errors
* Add simulated RFI
* Excised and un-excised