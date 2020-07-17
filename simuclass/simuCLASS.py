import pdb
import numpy as np
import sys
import os
#import pyfits as fits
import ConfigParser
import subprocess

config = ConfigParser.ConfigParser()
config.read(sys.argv[1])

# ToDo: a bit of arithmetic on the ra and dec to give to runskymodel (as new config parameter)


if config.getboolean('pipeline', 'doskymodel'):
  from skymodel.skymodel import runSkyModel
  runSkyModel(config)

if config.getboolean('pipeline', 'dosimdata'):
  from simulatedata.simdata import runSimulateData
  runSimulateData(config)#, sys.argv[1])

if config.getboolean('pipeline', 'doimagedata'):
  from imager.imager import runNWImager, runCASAClean, runWSClean
  if config.get('imager', 'type') == 'casaclean':
    runCASAClean(config, sys.argv[1])
  elif config.get('imager', 'type') == 'nwimager':
    runNWImager(config, sys.argv[1])
  elif config.get('imager', 'type') == 'wsclean':
    runWSClean(config)
    if config.getboolean('imager', 'dopostagestamps'):
      from thumbnailer.thumbnailer import makeThumbnails
      makeThumbnails(config)
  else:
    'You picked an unknown imager!'
