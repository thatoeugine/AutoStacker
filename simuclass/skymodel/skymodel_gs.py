'''
Script to convert a T-RECS catalogue into a sky model FITS file.
Forms part of the larger simuCLASS pipeline for creating simulated superCLASS
observations.

Usage:
python skymodel_gs.py example.ini

Prequisities:
skymodel_tools.py
galsim
astropy

Contact:
Ian Harrison
ian.harrison-2@manchester.ac.uk

'''
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

def runSkyModel(config):
  '''Simulate a sky model from a T-RECS catalogue.

  Parameters
  ----------
  config : configparser
    ConfigParser configuration containing necessary sections.

  '''

  data_path = config.get('pipeline', 'data_path')

  # Set some image properties
  pixel_scale = config.getfloat('skymodel', 'pixel_scale')*galsim.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')*galsim.arcmin
  image_size = int((fov/galsim.arcmin)/(pixel_scale/galsim.arcmin))

  ra_field = config.get('field', 'field_ra')
  ra_field_gs = galsim.HMS_Angle(ra_field)
  dec_field = config.get('field', 'field_dec')
  dec_field_gs = galsim.DMS_Angle(dec_field)
  
  # Load the catalogue
  cat_file_name = config.get('field', 'catalogue')
  print('Loading catalogue from {0} ...'.format(cat_file_name))
  cat =  fits.getdata(cat_file_name)

  
  # Set up a WCS for the catalogue
  cat_wcs = ast_wcs.WCS(naxis=2)
  cat_wcs.wcs.crpix = [image_size/2, image_size/2]
  cat_wcs.wcs.cdelt = [pixel_scale/galsim.degrees, pixel_scale/galsim.degrees]
  cat_wcs.wcs.crval = [0.e0, 0.e0]
  cat_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
  
  cat_fov_max = float(fov/galsim.arcmin)/(2.*60.) # cat is generated on 1 degsq
  fov_cut = (abs(cat['latitude']) < cat_fov_max)*(abs(cat['longitude']) < cat_fov_max)
  
  #pdb.set_trace()
  
  cat = cat[fov_cut]
  if config.getboolean('skymodel', 'highfluxcut'):
    highflux_cut = cat['I1400']*1.e-3 < 500.e-6
    cat = cat[highflux_cut]
  if config.getboolean('skymodel', 'lowfluxcut'):
    lowflux_cut = cat['I1400']*1.e-3 > 25.e-6
    cat = cat[lowflux_cut]
  if config.getboolean('skymodel', 'highsizecut'):
    highsize_cut = cat['size']/2. < 10
    cat = cat[highsize_cut]
  if config.getboolean('skymodel', 'lowsizecut'):
    lowsize_cut = cat['size']/2. > 0.75
    cat = cat[lowsize_cut]
  
  gal_ra = cat['latitude']
  gal_dec = cat['longitude']
  gal_e1 = cat['e1']
  gal_e2 = cat['e2']
  gal_flux = cat['I1400']*1.e-3 #mjy convert to Jy
  gal_r0 = cat['size']/2.# ????? Factor of 2?
  g1 = 0
  g2 = 0
  
  if config.get('skymodel', 'fluxscale')=='constant':
    gal_flux = np.ones_like(gal_flux)*100e-6
  
  nobj = len(cat)
  if config.getint('skymodel', 'ngals') > -1:
    nobj = config.getint('skymodel', 'ngals')
  
  ix_arr = np.ones(nobj)
  iy_arr = np.ones(nobj)
  print('...done.')
  
  # Create the galsim image
  full_image = galsim.ImageF(image_size, image_size, scale=pixel_scale)
  im_center = full_image.bounds.trueCenter()
  sky_center = galsim.CelestialCoord(ra=ra_field_gs, dec=dec_field_gs)

  # Create a WCS for the galsim image
  dudx = -pixel_scale / galsim.arcsec # - on dx's since the ra axis is flipped.
  dudy = 0.
  dvdx = 0.
  dvdy = pixel_scale / galsim.arcsec
  image_center = full_image.trueCenter()
  affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=full_image.trueCenter())
  wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
  full_image.wcs = wcs

  tstart=time.time()
  
  # Draw the galaxies onto the galsim image
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
    ix_arr[i] = ix
    iy_arr[i] = iy
    offset = galsim.PositionD(x-ix, y-iy)
    
    # Create the sub-image for this galaxy
    stamp = gal.drawImage(scale=pixel_scale/galsim.arcsec, offset=offset)
    stamp.setCenter(ix, iy)
    
    # Add the sub-image to the full iamge
    bounds = stamp.bounds & full_image.bounds
    full_image[bounds] += stamp[bounds]
    sys.stdout.flush()
  
  tend = time.time()
  print('\n...done in {0} seconds.'.format(tend-tstart))
  all_gals_fname = data_path+config.get('field', 'fitsname')
  print('Writing image data to {0} ...'.format(all_gals_fname))

  # Extract the numpy array from the galsim image
  image_data = full_image.array

  # Write out the image with the 4D FITS header correct for e.g. casa simulation
  write4dImage(all_gals_fname, image_data,
               pixel_scale / galsim.degrees,
               obs_ra=ra_field_gs / galsim.degrees,
               obs_dec=dec_field_gs / galsim.degrees,
               obs_freq=config.getfloat('observation', 'lowest_frequency'))
  print('...done.')
  
  if config.getboolean('skymodel', 'im3cat'):
    np.savetxt(config.get('pipeline', 'data_path')+'im3cat.txt', np.column_stack([np.arange(nobj), ix_arr, iy_arr]))
  
  print('runSkyModel complete.')

if __name__=='__main__':

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[1])
  
  runSkyModel(config)
