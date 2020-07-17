import numpy as np

import galsim

def pb_aips(r, nu, pbparm=[0,1, -4.968, 99.99, -72.78, 0.72, 7.806]):
  '''
  r is in arcmin
  nu is in hz
  '''
  x = (r*nu/1.e9)**2.
  
  retVar = 1.e0 + x*pbparm[2]/1.e3 + x**2.*pbparm[3]/1.e7 + \
           x**3.*pbparm[4]/1.e10 + x**4.*pbparm[5]/1.e13 + x**5.*pbparm[6]/1.e16
           
  return retVar
  
def primary_beam(config, freq):

  parm_dict = {
               'jackson' : [0,1, -4.968, 99.99, -72.78, 0.72, 7.806],
               'jvla' : [0,1, -1.343, 6.579, -1.186, 0., 0.]
               }
  
  pixel_scale = config.getfloat('skymodel', 'pixel_scale')*galsim.arcsec
  fov = config.getfloat('skymodel', 'field_of_view')*galsim.arcmin
  model_name = config.get('primarybeam', 'pb_model')
  image_size = int((fov/galsim.arcmin)/(pixel_scale/galsim.arcmin))
  
  fov_rad = 0.5*fov / galsim.arcmin
  
  x_arr = np.linspace(-fov_rad, fov_rad, image_size)
  
  xx, yy = np.meshgrid(x_arr, x_arr)
  
  rr = np.sqrt(xx**2 + yy**2)
  
  pb = pb_aips(rr, freq, pbparm=parm_dict[model_name])
  
  return pb
  
if __name__=='__main__':
  
  x_arr = np.linspace(-7.5,7.5,1024)
  
  xx, yy = np.meshgrid(x_arr, x_arr)
  
  rr = np.sqrt(xx**2 + yy**2)
  
  pb = pb_aips(rr, 1.4e9)
