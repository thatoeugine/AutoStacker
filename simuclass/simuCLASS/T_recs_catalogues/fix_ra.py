"""
Code edits the original T-RECS catalogues by:
    1) subtrating all the RA/longitudes greater than 0.5 deg by 360 deg
    2) Adding some additiional parameters gamma1 and gamma2 for simuclass friendly
"""

import numpy as np
import sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as pl
import astropy.io.fits as pf

input_name = sys.argv[1]

output_name = input_name.replace(".txt","_flipped_ra.txt")

cat_read =  pf.getdata(input_name)


RA = np.zeros(len(cat_read["longitude"]))


for i in range(len(cat_read["size"])):
    ra = cat_read["longitude"][i]
    if ra > 0.5:
        ra-=360
        RA[i] = ra
    else:
        RA[i] = ra

#### Plots for checking if the edit worked###        
pl.hist(cat_read["longitude"])
pl.savefig("before.png")

pl.figure()
pl.hist(RA)
pl.savefig("after.png")
############################################


if input_name == 'catalogue_SFGs_complete_deep.fits':

    outfile = open('catalogue_SFGs.txt', 'w')

    print>> outfile, '#lon            lat            size           flux           redshift       e1             e2             gamma1         gamma2'


    for i in range(len(RA)):
        srcline ='%.10f           %.10f            %.10f           %.10f           %.10f       %.10f             %.10f             %.10f         %.10f'%(RA[i],cat_read['latitude'][i],cat_read['size'][i],cat_read['I1400'][i],cat_read['redshift'][i],cat_read['e1'][i],cat_read['e2'][i],0,0)
        print>> outfile, srcline

    outfile.close()


elif input_name == 'catalogue_AGNs_complete_deep.fits':
    
    outfile = open('catalogue_AGNs.txt', 'w')

    print>> outfile, '#lon(deg)            lat(deg)            size(arcsec)           flux(mJy)           redshift       e1             e2             gamma1         gamma2            Rs'


    for i in range(len(RA)):
        srcline ='%.10f           %.10f            %.10f           %.10f           %.10f       %.10f             %.10f             %.10f         %.10f            %.10f'%(RA[i],cat_read['latitude'][i],cat_read['size'][i],cat_read['I1400'][i],cat_read['redshift'][i],0,0,0,0,cat_read['Rs'][i])
        print>> outfile, srcline

    outfile.close()
    
else:
    print "Error unknown catalogue"
    sys.exit()
