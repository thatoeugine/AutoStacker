import numpy as np
import ConfigParser
import sys
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import wcs as ast_wcs
import pdb
import galsim

from astropy import units as uns

def makeThumbnails(config):

  # read in catalogue
  # extract thumbnails at all points in catalogue

  data_path = config.get('pipeline', 'data_path')

  # Set some image properties
  pixel_scale = config.getfloat('skymodel', 'pixel_scale')*galsim.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')*galsim.arcmin
  image_size = int((fov/galsim.arcmin)/(pixel_scale/galsim.arcmin))

  ra_field = config.get('field', 'field_ra')
  ra_field_gs = galsim.HMS_Angle(ra_field)
  dec_field = config.get('field', 'field_dec')
  dec_field_gs = galsim.DMS_Angle(dec_field)
  
  # Load the wcs
  header = fits.getheader(config.get('pipeline', 'data_path')+config.get('field', 'fitsname'))
  w_twod = ast_wcs.WCS(header)
  
  # Load the catalogue
  cat_file_name = config.get('field', 'catalogue')
  print('Loading catalogue from {0} ...'.format(cat_file_name))
  cat =  fits.getdata(cat_file_name)
  ra_offset = cat['longitude']
  dec_offset = cat['latitude']
  
  ra_offset_max = 0.9*(fov/2) / galsim.degrees
  dec_offset_max = 0.9*(fov/2) / galsim.degrees
  
  fov_cut = (abs(ra_offset) < ra_offset_max)*(abs(dec_offset) < dec_offset_max)
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
  
  gal_ra_offset = cat[' longitude']
  gal_dec_offset = cat['latitude']
  
  gal_dec = dec_field_gs / galsim.degrees + gal_dec_offset
  gal_dec_radians = (gal_dec*galsim.degrees) / galsim.radians
  gal_ra = ra_field_gs / galsim.degrees + gal_ra_offset/np.cos(np.asarray(gal_dec_radians, dtype=float))
  
  gal_e1 = cat['e1']
  gal_e2 = cat['e2']
  gal_flux = cat['I1400']*1.e-3 #mjy convert to Jy
  gal_r0 = cat['size']/2.# ????? Factor of 2?
  g1 = 0
  g2 = 0
  
  
  if config.get('skymodel', 'fluxscale')=='constant':
    gal_flux = np.ones_like(gal_flux)*np.median(gal_flux)*1.e4
  
  nobj = len(cat)
  if config.getint('skymodel', 'ngals') > 0:
    nobj = config.getint('skymodel', 'ngals')
  
  print('...done.')

  #pdb.set_trace()

  for i in range(nobj):
  
    plt.close('all')
    
    sys.stdout.write('\rCreating thumbnail for source {0} of {1}...'.format(i+1, nobj))

    npix_stamp = 10*int((gal_r0[i])/(pixel_scale/galsim.arcsec))
    x, y, _, _ = w_twod.wcs_world2pix(gal_ra[i], gal_dec[i], 0, 0, 0)
    x = float(x)
    y = float(y)
    
    # Account for the fractional part of the position:
    ix = int(np.floor(x+0.5))
    iy = int(np.floor(y+0.5))

    plt.figure(i, figsize=(9, 7.5))

    plt.subplot(221)
    im = fits.getdata(config.get('pipeline', 'data_path')+config.get('field', 'fitsname'))
    im = im[0,0]

    stamp = im[iy-npix_stamp/2:iy+npix_stamp/2,ix-npix_stamp/2:ix+npix_stamp/2]
    plt.imshow(stamp, cmap='afmhot', origin='lower', interpolation='nearest')

    plt.subplot(222)
    im = fits.getdata(config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'.wsclean-dirty.fits')

    im = im[0,0,:,::-1]
    stamp = im[iy-npix_stamp/2:iy+npix_stamp/2,ix-npix_stamp/2:ix+npix_stamp/2]
    plt.imshow(stamp, cmap='afmhot', origin='lower', interpolation='nearest')
    #pdb.set_trace()
    plt.subplot(223)
    if config.get('imager','type')=='wsclean':
      im = fits.getdata(config.get('pipeline', 'data_path')+\
                config.get('pipeline', 'project_name')+'.wsclean-image.fits')
      im = im[0,0,:,::-1]
    elif config.get('imager','type')=='casa':
      im = fits.getdata(config.get('pipeline', 'data_path')+\
                config.get('pipeline', 'project_name')+'.casa-cs.image.fits')
      im = im[0,0,:,::-1]
      size_diff = im.shape[0] - image_size
      im = im[size_diff/2:-size_diff/2,size_diff/2:-size_diff/2]
                
    stamp = im[iy-npix_stamp/2:iy+npix_stamp/2,ix-npix_stamp/2:ix+npix_stamp/2]
    plt.imshow(stamp, cmap='afmhot', origin='lower', interpolation='nearest')
    plt.subplot(224)
    if config.get('imager','type')=='wsclean':
      im = fits.getdata(config.get('pipeline', 'data_path')+\
                config.get('pipeline', 'project_name')+'.wsclean-model.fits')
      im = im[0,0,:,::-1]
    elif config.get('imager','type')=='casa':
      im = fits.getdata(config.get('pipeline', 'data_path')+\
                config.get('pipeline', 'project_name')+'.casa-cs.model.fits')
      im = im[0,0,:,::-1]
      size_diff = im.shape[0] - image_size
      im = im[size_diff/2:-size_diff/2,size_diff/2:-size_diff/2]

    stamp = im[iy-npix_stamp/2:iy+npix_stamp/2,ix-npix_stamp/2:ix+npix_stamp/2]
    plt.imshow(stamp, cmap='afmhot', origin='lower', interpolation='nearest')
    plt.savefig(config.get('pipeline', 'data_path')+config.get('pipeline', 'project_name')+'source_{0}.png'.format(i), bbox_inches='tight', dpi=160)
    
    sys.stdout.flush()
    
  print('\n...done.')


if __name__=='__main__':

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[1])

  makeThumbnails(config)


