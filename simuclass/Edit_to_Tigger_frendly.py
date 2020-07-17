"""
Code edits the output catalogue from simuclass to Tigger friendly.
Note code is only made for the SFGs catalogue for now....!!! 
"""


import numpy as np
import sys
import astropy.io.fits as pf
from astropy.table import Table

input_name = sys.argv[1]

output_name = input_name.replace(".txt","_fixed.txt")

cat_read =  Table.read(input_name,format='ascii')

#Major  axis equals to the size of the disc: 
Major_axis = cat_read['size']
        
#Minor axis & Position angle determined using ellipse formular from Tunbridge et al. (2016):   
e = np.sqrt(cat_read['e1']**2 + cat_read['e2']**2) # absolute ellipticity   
        
Minor_axis = 2*((Major_axis/2)*np.sqrt((1-e)/(e+1)))
    
# Position Angle of Galaxy [degrees]:
PA= np.random.uniform(0,180,len(cat_read['size']))


outfile = open(output_name, 'w')
print>> outfile, '#format: name ra_d dec_d i emaj_s emin_s pa_d'


for i in range(len(cat_read['ra_abs'])):
    srcline = 'SrcGauss%04d %s %s %.10f %.10f %.10f %.10f'\
    %(i,cat_read['ra_abs'][i],cat_read['dec_abs'][i],cat_read['integrated_flux'][i],Major_axis[i],Minor_axis[i],PA[i])

    print>> outfile, srcline

outfile.close()
