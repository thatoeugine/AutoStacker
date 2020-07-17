import numpy as np
from primarybeam.primarybeam import *
from skymodel.skymodel_tools import setup_wcs

from astropy.io import fits
from astropy import units as uns

import ConfigParser
import sys

if __name__=='__main__':

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[1])
  
  nu_min = config.getfloat('observation', 'lowest_frequency')
  nu_width = config.getfloat('observation', 'total_bandwidth')
  n_chan = config.getint('observation', 'n_channels')
  delta_nu = nu_width/n_chan
  
  pix_scale = config.getfloat('skymodel', 'pixel_scale')*uns.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')*uns.arcmin
  
  n_pix = fov.to(uns.arcsec).value / pix_scale.value
    
  ra_arr = np.linspace(-fov.value/2,fov.value/2,n_pix)
  dec_arr = np.linspace(-fov.value/2,fov.value/2,n_pix)
  
  xx, yy = np.meshgrid(ra_arr, dec_arr)
  rr = np.sqrt(xx**2 + yy**2)
  
  w_fourd = setup_wcs(config, ndim=4, nu_axis=True)
  output_header = w_fourd.to_header()
  output_cube_fname = config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'pb_cube.fits'
  
  for i in range(n_chan):
    sys.stdout.write('\rGenerating channel {0} of {1}...'.format(i+1, n_chan))
    nu = nu_min + i*delta_nu
    pb_image = pbfunc(rr, nu)
    
    fits.append(output_cube_fname, pb_image, output_header)
    sys.stdout.flush()
  print('\n...done')
