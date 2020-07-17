'''
Code to derive AGN properties for drawing from Anna Bonaldi's
T-RECS catalogue.
'''
import numpy as np
import pandas as pd
from astropy.cosmology import Planck13 as cosm
import astropy.units as uns
import galsim
'''
data_dir = '/local/scratch/harrison/simuCLASS/catalogues/'

sfg_cat = pd.read_csv(data_dir+'catalogue_SFGs_complete_v2.fits.txt', header=0, delim_whitespace=True)

agn_cat = pd.read_csv(data_dir+'catalogue_AGNs_complete_v2.fits.txt', header=0, delim_whitespace=True)

agn_projected_physical_size = agn_cat['sizeKpc']*np.sin(agn_cat['angle(deg)'])*uns.kpc

agn_projected_angular_size = agn_projected_physical_size / cosm.angular_diameter_distance(agn_cat['redshift'])
'''

image_size = 64
pixel_scale = 0.05

gal_flux = 1
gal_r0 = 1
g1 = 0
g2 = 0
rs = 0.8
posang = np.pi/4.

full_image = galsim.ImageF(image_size, image_size, scale=pixel_scale)

lobe_flux = gal_flux*0.99
hs_flux = gal_flux - lobe_flux
hs1_flux = hs_flux/3.
hs2_flux = hs_flux/3.
hs3_flux = hs_flux/3.
      
hs_offset = rs*gal_r0
lobe_offset = gal_r0*0.6
      
lobe1 = galsim.Gaussian(sigma=gal_r0*0.25, flux=lobe_flux/2.)
lobe2 = galsim.Gaussian(sigma=gal_r0*0.25, flux=lobe_flux/2.)

lobe1 = lobe1.shear(e1=0.3,e2=0)
lobe2 = lobe2.shear(e1=0.3,e2=0)

lobe1 = lobe1.shift(-lobe_offset,0)
lobe2 = lobe2.shift(lobe_offset,0)

gal = lobe1 + lobe2

gal = gal.rotate(posang*galsim.radians)      

total_shear = galsim.Shear(g1=g1, g2=g2)      
      
gal = gal.shear(total_shear)

x = image_size/2
y = image_size/2
x = float(x)
y = float(y)

# Account for the fractional part of the position:
ix = int(np.floor(x+0.5))
iy = int(np.floor(y+0.5))
ix_arr = ix
iy_arr = iy
offset = galsim.PositionD(x-ix, y-iy)
hs_offset_pixels = hs_offset*pixel_scale
hs_ix_offset = hs_offset*np.sin(posang[i])/pixel_scale
hs_iy_offset = hs_offset*np.cos(posang[i])/pixel_scale

# Create the sub-image for this galaxy
stamp = gal.drawImage(scale=pixel_scale, offset=offset)
stamp.setCenter(ix, iy)

# Add the sub-image to the full iamge
bounds = stamp.bounds & full_image.bounds
full_image[bounds] += stamp[bounds]

cen = stamp.array.shape
# Add the hotspots as single pixel point sources
stamp.array[cen[0]/2, cen[1]/2] += hs1_flux
stamp.array[cen[0]/2+hs_ix_offset, cen[1]/2+hs_iy_offset] += hs2_flux
stamp.array[cen[0]/2-hs_ix_offset, cen[1]/2-hs_iy_offset] += hs3_flux

'''
for gal in agn_projected_angular_size:
  lobe1 = galsim.Gaussian(sigma=lober0, flux=lobeflux)
  lobe2 = galsim.Gaussian(sigma=lober0, flu=lobeflux)
  
  hs1 = Galsim.TopHat(1.e-10, flux=centralflux)
  hs2 = Galsim.TopHat(1.e-10, flux=hsflux)
  hs3 = Galsim.TopHat(1.e-10, flux=hsflux)
  
  lobe1 = lobe1.shear(e1=lobe_e1, e2=lobe_e2)
  lobe2 = lobe2.shear(e1=lobe_e1, e2=lobe_e2)
  lobe1 = lobe1.shift(lobe_offset, 0)
  lobe2 = lobe2.shift(-lobe_offset, 0)
  
  hs2.Shift(hs_offset, 0)
  hs3.Shift(-hs_offset, 0)
  
  lobes = lobe1 + lobe2
  
  lobes_image = lobes.rotate(posang*galsim.radians)
  
  lobes.drawImage(scale=pixel_scale)
  
  spots = hs1 + hs2 + hs3
  
  spots = spots.rotate(posang*galsim.radians)
  
  spots_image = spots.drawImage(scale=pixel_scale)
  
  full_image = spots_image + lobes_image
'''
