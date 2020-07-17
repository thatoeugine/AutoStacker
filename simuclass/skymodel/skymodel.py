'''
Script to convert a T-RECS catalogue into a sky model FITS file.
Forms part of the larger simuCLASS pipeline for creating simulated superCLASS
observations.

Usage:
python skymodel.py example.ini

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
import cPickle as pickle

import galsim

from astropy.io import fits
from astropy.table import Table
from astropy import units as uns
from astropy import wcs as ast_wcs

from scipy import random
from matplotlib import pyplot as plt

from skymodel_tools import setup_wcs
from skymodel_tools import write4dImage
from primarybeam.primarybeam import *

big_fft_params = galsim.GSParams(maximum_fft_size=81488)
arcsectorad = (1.*uns.arcsec).to(uns.rad).value
degtoarcsec = (1.*uns.deg).to(uns.arcsec).value

# theta is in trecs catalogue
# Rs > 0.5 is FRII and is in trecs catalogue

  
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
  
  w_twod = setup_wcs(config, ndim=2)
  w_fourd = setup_wcs(config, ndim=4)
  header_twod = w_twod.to_header()
  header_fourd = w_fourd.to_header()
  
  # Load the catalogue
  cat_file_name = config.get('field', 'catalogue')
  print('Loading catalogue from {0} ...'.format(cat_file_name))
  cat = Table()
  if config.get('skymodel', 'catalogue_origin') == 'trecs':
    cat_read =  Table.read(cat_file_name, format='ascii')
    
    cat['ra_offset'] = cat_read['lon'] # deg
    cat['ra_offset'].unit = 'deg'
    
    cat['dec_offset'] = cat_read['lat'] # deg
    cat['dec_offset'].unit = 'deg'

    
    cat['dec_abs'] = dec_field_gs / galsim.degrees + cat['dec_offset']
    dec_abs_radians = cat['dec_abs']*galsim.degrees / galsim.radians
    cat['ra_abs'] = ra_field_gs / galsim.degrees + cat['ra_offset']/np.cos(np.asarray(dec_abs_radians, dtype=float))

    cat['integrated_flux'] = cat_read['flux']*1.e-3 # Jy
    cat['integrated_flux'].unit = 'Jy'

    cat['size'] = cat_read['size'] # arcsec
    cat['size'].unit = 'arcsec'

    cat['peak_flux'] = cat['integrated_flux'] / (2.*cat['size']*arcsectorad)
    cat['peak_flux'].unit = 'Jy'

    cat['e1'] = cat_read['e1']
    cat['e2'] = cat_read['e2']

    cat['g1'] = cat_read['gamma1'] # 0
    cat['g2'] = cat_read['gamma2'] # 0

  elif config.get('skymodel', 'catalogue_origin') == 'pybdsm':
    cat_read =  Table.read(cat_file_name, format='fits')

    cat['ra_abs'] = cat_read['RA'] # deg
    cat['dec_abs'] = cat_read['DEC'] # deg

    cat['dec_offset'] = cat['dec_abs'] - dec_field_gs / galsim.degrees
    dec_abs_radians = cat['dec_abs']*galsim.degrees / galsim.radians
    cat['ra_offset'] = (cat['ra_abs'] - ra_field_gs / galsim.degrees)*np.cos(np.asarray(dec_abs_radians, dtype=float))
    
    cat['integrated_flux'] = cat_read['Total_flux'] # Jy

    cat['size'] = cat_read['Maj']*degtoarcsec # deg
    cat['size'].unit = 'arcsec'

    cat['peak_flux'] = cat_read['Peak_flux'] # Jy

    cat['q'] = cat_read['Min']/cat_read['Maj']
    cat['position_angle'] = np.arctan2(cat_read['Min'], cat_read['Maj'])
    cat['mod_e'] = (1. - cat['q']**2.)/(1. + cat['q']**2.)

    cat['e1'] = cat['mod_e']*np.cos(2.*cat['position_angle'])
    cat['e2'] = cat['mod_e']*np.sin(2.*cat['position_angle'])

    cat['g1'] = 0.e0
    cat['g2'] = 0.e0

  # fov cut
  ra_offset_max = 0.9*(fov/2) / galsim.degrees
  dec_offset_max = 0.9*(fov/2) / galsim.degrees
  
  fov_cut = (abs(cat['ra_offset']) < ra_offset_max)*(abs(cat['dec_offset']) < dec_offset_max)
  cat = cat[fov_cut]
  
  # flux cuts
  if config.getboolean('skymodel', 'highfluxcut'):
    highflux_cut = cat['peak_flux'] < config.getfloat('skymodel', 'highfluxcut_value')
    cat = cat[highflux_cut]
  
  if config.getboolean('skymodel', 'lowfluxcut'):
    lowflux_cut = cat['peak_flux'] > config.getfloat('skymodel', 'lowfluxcut_value')
    cat = cat[lowflux_cut]
  
  if config.getboolean('skymodel', 'highsizecut'):
    highsize_cut = cat['size'] < config.getfloat('skymodel', 'highsizecut_value')
    cat = cat[highsize_cut]
  
  if config.getboolean('skymodel', 'lowsizecut'):
    lowsize_cut = cat['size'] > config.getfloat('skymodel', 'lowsizecut_value')
    cat = cat[lowsize_cut]
    
  if config.get('skymodel', 'sizescale')=='constant':
    cat['size'] = np.ones_like(cat['size'])*config.getfloat('skymodel', 'sizescale_constant_value')

  # number of sources, on grid if requested
  if config.getboolean('skymodel', 'grid'):
    nobj = int(np.sqrt(config.getint('skymodel', 'ngals')))**2.
    cat['ra_offset'] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
    cat['dec_offset'] = np.linspace(-ra_offset_max, ra_offset_max, nobj)
  else:
    nobj = len(cat)
    if config.getint('skymodel', 'ngals') > -1:
      nobj = config.getint('skymodel', 'ngals')
      cat = cat[:nobj]
  
  # flux range
  if config.get('skymodel', 'fluxscale')=='constant':
    cat['integrated_flux'] = np.ones_like(cat['integrated_flux'])*config.getfloat('skymodel', 'fluxscale_constant_value')
    cat['peak_flux'] = cat['integrated_flux'] / (2.*cat['size']*arcsectorad)
    
  # scale flux
  cat['integrated_flux'] = cat['integrated_flux']*config.getfloat('skymodel', 'flux_factor')
  cat['peak_flux'] = cat['peak_flux']*config.getfloat('skymodel', 'flux_factor')
  
  # write out catalogue
  cat.write(config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'_truthcat.txt', format='ascii')
  
  ix_arr = np.ones(nobj)
  iy_arr = np.ones(nobj)
  print('...done.')
  
  # Create the galsim image
  full_image = galsim.ImageF(image_size, image_size, scale=pixel_scale/galsim.arcsec)
  im_center = full_image.bounds.trueCenter()
  sky_center = galsim.CelestialCoord(ra=ra_field_gs, dec=dec_field_gs)

  # Create a WCS for the galsim image
  full_image.wcs, origin = galsim.wcs.readFromFitsHeader(header_twod)

  tstart=time.time()
  
  # Draw the galaxies onto the galsim image
  for i,cat_gal in enumerate(cat):
    
    sys.stdout.write('\rAdding source {0} of {1} to skymodel...'.format(i+1, nobj))
    
    # choose the profile
    if config.get('skymodel', 'galaxy_profile')=='exponential':
      gal = galsim.Exponential(scale_radius=cat_gal['size']/2., flux=cat_gal['integrated_flux'], gsparams=big_fft_params)
    
    elif config.get('skymodel', 'galaxy_profile')=='gaussian':
      gal = galsim.Gaussian(fwhm=cat_gal['size'], flux=cat_gal['integrated_flux'], gsparams=big_fft_params)
    
    elif config.get('skymodel', 'galaxy_profile')=='matched-exponential':
      gauss_gal = galsim.Gaussian(fwhm=cat_gal['size'], flux=cat_gal['integrated_flux'])
      gal = galsim.Exponential(half_light_radius=gauss_gal.getHalfLightRadius(), flux=cat_gal['integrated_flux'], gsparams=big_fft_params)
      del gauss_gal

    # calculate the total ellipticity
    ellipticity = galsim.Shear(e1=cat_gal['e1'],e2=cat_gal['e2'])
    shear = galsim.Shear(g1=cat_gal['g1'],g2=cat_gal['g2'])
    if config.getboolean('skymodel', 'doshear'):
      total_shear = ellipticity + shear
    else:
      total_shear = ellipticity

    gal = gal.shear(total_shear)
    
    x, y = w_twod.wcs_world2pix(cat_gal['ra_abs'], cat_gal['dec_abs'], 0,)
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
  
  if config.getboolean('skymodel', 'doagn'):
    # Load the catalogue
    cat_file_name = config.get('field', 'agncatalogue')
    print('Loading catalogue from {0} ...'.format(cat_file_name))
    cat =  Table.read(cat_file_name, format='ascii')

    ra_offset = cat['lon(deg)']
    dec_offset = cat['lat(deg)']
    
    ra_offset_max = 0.9*(fov/2) / galsim.degrees
    dec_offset_max = 0.9*(fov/2) / galsim.degrees
    
    fov_cut = (abs(ra_offset) < ra_offset_max)*(abs(dec_offset) < dec_offset_max)
    #pdb.set_trace()
    cat = cat[fov_cut]
    
    if config.getboolean('skymodel', 'highfluxcut'):
      highflux_cut = cat['flux(mJy)']*1.e-3 < 100.e-6
      cat = cat[highflux_cut]
    if config.getboolean('skymodel', 'lowfluxcut'):
      lowflux_cut = cat['flux(mJy)']*1.e-3 > 25.e-6
      cat = cat[lowflux_cut]
    if config.getboolean('skymodel', 'highsizecut'):
      highsize_cut = cat['size(arcsec)']/2. < 10
      cat = cat[highsize_cut]
    if config.getboolean('skymodel', 'lowsizecut'):
      lowsize_cut = cat['size(arcsec)']/2. > 0.75
      cat = cat[lowsize_cut]
    
    if config.getboolean('skymodel', 'grid'):
      nobj = int(np.sqrt(config.getint('skymodel', 'ngals')))**2.
      
      gal_ra_offset = np.linspace(-ra_offset_max, ra_offset_max, nobj)
      gal_dec_offset = np.linspace(-ra_offset_max, ra_offset_max, nobj)
    else:
      gal_ra_offset = cat['lon(deg)']
      gal_dec_offset = cat['lat(deg)']
      nobj = len(cat)
      if config.getint('skymodel', 'ngals') > -1:
        nobj = config.getint('skymodel', 'ngals')
    
    gal_dec = dec_field_gs / galsim.degrees + gal_dec_offset
    gal_dec_radians = (gal_dec*galsim.degrees) / galsim.radians
    gal_ra = ra_field_gs / galsim.degrees + gal_ra_offset/np.cos(np.asarray(gal_dec_radians, dtype=float))
    
    gal_e1 = cat['e1']
    gal_e2 = cat['e2']
    gal_flux = cat['flux(mJy)']*1.e-3 #mjy convert to Jy
    gal_size = cat['size(arcsec)']/(2.)
    g1 = cat['gamma1']
    g2 = cat['gamma2']
    rs = cat['Rs']
    
    if config.get('skymodel', 'fluxscale')=='constant':
      gal_flux = np.ones_like(gal_flux)*100e-6
    
    ra_abs, dec_abs = w_twod.wcs_world2pix(gal_ra, gal_dec, 0,)
    truthcat = np.column_stack([gal_ra, gal_dec, ra_abs, dec_abs, gal_e1, gal_e2, gal_flux, gal_size, g1, g2])
    np.savetxt(config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'-agn_truthcat.txt', truthcat[:nobj])
    
    ix_arr = np.ones(nobj)
    iy_arr = np.ones(nobj)
    posang = random.uniform(0,2.*np.pi, nobj)
    print('...done.')
    
    # Draw the galaxies onto the galsim image
    for i in range(nobj):
    
      sys.stdout.write('\rAdding agn source {0} of {1} to skymodel...'.format(i+1, nobj))
      
      if (rs[i] < 0.01) or (cat_gal['size'] < config.getfloat('skymodel', 'pixel_scale')/2):
        x, y = w_twod.wcs_world2pix(cat_gal['ra_abs'], cat_gal['dec_abs'], 0,)
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
        
        cen = stamp.array.shape
        # Add the hotspots as single pixel point sources
        stamp.array[cen[0]/2, cen[1]/2] += cat_gal['integrated_flux']
        
        # Add the sub-image to the full iamge
        bounds = stamp.bounds & full_image.bounds
        full_image[bounds] += stamp[bounds]
      
      else:
        lobe_flux = cat_gal['integrated_flux']*0.99
        hs_flux = cat_gal['integrated_flux'] - lobe_flux
        hs1_flux = hs_flux/3.
        hs2_flux = hs_flux/3.
        hs3_flux = hs_flux/3.
              
        hs_offset = rs[i]*cat_gal['size']
        lobe_offset = cat_gal['size']*0.6
        
        lobe1 = galsim.Gaussian(sigma=cat_gal['size']*0.25, flux=lobe_flux/2., gsparams=big_fft_params)
        lobe2 = galsim.Gaussian(sigma=cat_gal['size']*0.25, flux=lobe_flux/2., gsparams=big_fft_params)

        lobe1 = lobe1.shear(e1=0.3,e2=0)
        lobe2 = lobe2.shear(e1=0.3,e2=0)

        lobe1 = lobe1.shift(-lobe_offset,0)
        lobe2 = lobe2.shift(lobe_offset,0)

        gal = lobe1 + lobe2

        gal = gal.rotate(posang[i]*galsim.radians)      

        total_shear = galsim.Shear(g1=cat_gal['g1'], g2=cat_gal['g2'])      
              
        gal = gal.shear(total_shear)
        
        x, y = w_twod.wcs_world2pix(cat_gal['ra_abs'], cat_gal['dec_abs'], 0,)
        x = float(x)
        y = float(y)
        
        # Account for the fractional part of the position:
        ix = int(np.floor(x+0.5))
        iy = int(np.floor(y+0.5))
        ix_arr[i] = ix
        iy_arr[i] = iy
        offset = galsim.PositionD(x-ix, y-iy)
        hs_offset_pixels = hs_offset*pixel_scale/galsim.arcsec
        hs_ix_offset = hs_offset*np.sin(posang[i]) / (pixel_scale/galsim.arcsec)
        hs_iy_offset = hs_offset*np.cos(posang[i]) / (pixel_scale/galsim.arcsec)
        
        # Create the sub-image for this galaxy
        stamp = gal.drawImage(scale=pixel_scale/galsim.arcsec, offset=offset)
        stamp.setCenter(ix, iy)
        
        cen = stamp.array.shape
        # Add the hotspots as single pixel point sources
        stamp.array[cen[0]/2, cen[1]/2] += hs1_flux
        stamp.array[cen[0]/2+hs_ix_offset, cen[1]/2+hs_iy_offset] += hs2_flux
        stamp.array[cen[0]/2-hs_ix_offset, cen[1]/2-hs_iy_offset] += hs3_flux
        
        if config.getboolean('skymodel', 'pickleagn'):
          pickle.dump(stamp.array, open('agn_{0}.p'.format(i), 'wb'))
        
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
  #write4dImage(all_gals_fname, image_data,
  #             pixel_scale / galsim.degrees,
  #             obs_ra=ra_field_gs / galsim.degrees,
  #             obs_dec=dec_field_gs / galsim.degrees,
  #             obs_freq=config.getfloat('observation', 'lowest_frequency'))
  
  if config.getboolean('primarybeam', 'dopb'):
    
    nstokes = config.getint('primarybeam', 'nstokes')
    nfreq = config.getint('primarybeam', 'nfreq')
    bw = config.getfloat('observation', 'total_bandwidth')
    base_freq = config.getfloat('observation', 'lowest_frequency')
    freq_width = bw / nfreq
  
    image_cube = np.empty((nstokes,nfreq)+image_data.shape)    
        
    stokes_list = ['I', 'Q', 'U'][:nstokes]
    freq_list = base_freq + np.arange(nfreq)*freq_width
    
    for i_stokes, stokes in enumerate(stokes_list):
      for i_freq, freq in enumerate(freq_list):
        image_cube[i_stokes,i_freq] = image_data*primary_beam(config, freq)

    hdu = fits.PrimaryHDU(image_cube, header=header_fourd)

  else:
    hdu = fits.PrimaryHDU(np.expand_dims(np.expand_dims(image_data, axis=0), axis=0), header=header_fourd)
  
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(all_gals_fname, clobber=True)
  
  print('...done.')
  
  if config.getboolean('skymodel', 'im3cat'):
    np.savetxt(config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'_im3cat.txt', np.column_stack([np.arange(nobj), ix_arr, iy_arr]))
    
  # Write out the image with the 4D FITS header correct for e.g. casa simulation
  write4dImage(all_gals_fname, image_data,
              pixel_scale / galsim.degrees,
               obs_ra=ra_field_gs / galsim.degrees,
               obs_dec=dec_field_gs / galsim.degrees,
               obs_freq=config.getfloat('observation', 'lowest_frequency'))
  
  print('runSkyModel complete.')

if __name__=='__main__':

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[1])
  
  runSkyModel(config)
    