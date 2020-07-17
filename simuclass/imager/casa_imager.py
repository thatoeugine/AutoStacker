import ConfigParser
import pdb
import os
if __name__=='__main__':
  config = ConfigParser.ConfigParser()
  config.read(sys.argv[-1])
  msname = config.get('pipeline', 'project_name')+'.ms'
  msname = config.get('pipeline', 'data_path')+msname
  imagename = msname+'.image'
  fov = config.getfloat('imager', 'field_of_view')
  pixel_scale = config.getfloat('imager', 'pixel_scale')
  niter_value  = config.getint('imager', 'number_of_iterations')
  cell_size =  [str(pixel_scale)+'arcsec']
  imsizes  =  (fov * 60.)/pixel_scale#config.getint('configuration', 'image_size')#(fov * 60.)/pixel_scale
  
  imsizes = 20050

  clean(vis = msname, imagename = msname+'.clean', niter = niter_value, imsize = [int(imsizes)+200, int(imsizes)+200], cell = cell_size, wprojplanes = -1)

  viewer(infile=msname+'.image',outfile=msname+'.image',outscale=-2.0,outformat="pdf",gui=False)
  
  exportfits(imagename=msname+'.image', fitsimage=msname+'.image.fits',
             overwrite=True)

  viewer(infile=msname+'.clean.image',outfile=msname+'.clean.image',outscale=-2.0,outformat="pdf",gui=False)
  
  exportfits(imagename=msname+'.clean.image', fitsimage=msname+'.clean.image.fits',
             overwrite=True)
