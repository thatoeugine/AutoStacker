import numpy as np
import pdb
import sys
import ConfigParser
import time

import galsim

from astropy.io import fits
from astropy.table import Table
from astropy import wcs as ast_wcs

from skymodel_tools import write4dImage

rseed = 123456


def runSkyModel(config):
  # image properties
  data_path = config.get('pipeline', 'data_path')
  pixel_scale = config.getfloat('skymodel', 'pixel_scale')*galsim.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')*galsim.arcmin
  image_size = int((fov/galsim.arcmin)/(pixel_scale/galsim.arcmin))

  ra_field = config.get('field', 'field_ra')
  ra_field_gs = galsim.HMS_Angle(ra_field)
  dec_field = config.get('field', 'field_dec')
  dec_field_gs = galsim.DMS_Angle(dec_field)
  
  cat_file_name = config.get('field', 'catalogue')
  print('Loading catalogue from {0} ...'.format(cat_file_name))
  cat =  fits.getdata(cat_file_name)
  nobj = len(cat)
  
  cat_wcs = ast_wcs.WCS(naxis=2)
  cat_wcs.wcs.crpix = [image_size/2, image_size/2]
  cat_wcs.wcs.cdelt = [pixel_scale/galsim.degrees, pixel_scale/galsim.degrees]
  cat_wcs.wcs.crval = [0.e0, 0.e0]
  cat_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
  
  
  gal_ra = cat['latitude']
  gal_dec = cat['longitude']
  gal_e1 = cat['e1']
  gal_e2 = cat['e2']
  gal_flux = cat['I1400'] #mjy
  gal_r0 = cat['size']/2.
  g1 = 0
  g2 = 0
  print('...done.')
  
  full_image = galsim.ImageF(image_size, image_size, scale=pixel_scale)
  im_center = full_image.bounds.trueCenter()
  sky_center = galsim.CelestialCoord(ra=ra_field_gs, dec=dec_field_gs)

  # - on dx's since the ra axis is flipped.
  dudx = -pixel_scale / galsim.arcsec
  dudy = 0.
  dvdx = 0.
  dvdy = pixel_scale / galsim.arcsec
  image_center = full_image.trueCenter()
  affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=full_image.trueCenter())
  wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
  full_image.wcs = wcs

  tstart=time.time()
  
  nobj=200
  
  for i in range(nobj):
  
    sys.stdout.write('\rAdding source {0} of {1} to skymodel...'.format(i+1, nobj))
    

    gal = galsim.Exponential(scale_radius=gal_r0[i], flux=gal_flux[i])

    ellipticity = galsim.Shear(e1=gal_e1[i],e2=gal_e2[i])
    shear = galsim.Shear(g1=g1[i],g2=g2[i])
    total_shear = ellipticity + shear

    gal = gal.shear(total_shear)
    
    x, y = cat_wcs.wcs_world2pix(gal_ra[i], gal_dec[i], 0)
    x = float(x)
    y = float(y)
    
    # Account for the fractional part of the position:
    ix = int(np.floor(x+0.5))
    iy = int(np.floor(y+0.5))
    offset = galsim.PositionD(x-ix, y-iy)
    
    stamp = gal.drawImage(scale=pixel_scale/galsim.arcsec, offset=offset)
    stamp.setCenter(ix, iy)
    
    bounds = stamp.bounds & full_image.bounds
    full_image[bounds] += stamp[bounds]
    sys.stdout.flush()
  
  tend = time.time()
  print('\n...done in {0} seconds.'.format(tend-tstart))
  all_gals_fname = data_path+config.get('field', 'fitsname')
  print('Writing image data to {0} ...'.format(all_gals_fname))
  image_data = full_image.array
  write4dImage(all_gals_fname, image_data,
               pixel_scale / galsim.degrees,
               obs_ra=ra_field_gs / galsim.degrees,
               obs_dec=dec_field_gs / galsim.degrees,
               obs_freq=config.getfloat('observation', 'lowest_frequency'))
  print('...done.')
  
  print('runSkyModel complete.')

if __name__=='__main__':

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[1])
  
  runSkyModel(config)
