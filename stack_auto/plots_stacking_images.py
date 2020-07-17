"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# Module plots and computes the uncertainities of the uv stacked image

"""
#--------------------------------------------------------------------------------------------

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import astropy.io.fits as pf
from matplotlib import pyplot as plt
plt.style.use("_classic_test")
from scipy import stats
import aplpy
import json
from matplotlib import rc,rcParams
rc('text', usetex=True)
# activate latex text rendering
rc('axes', linewidth=2)
rc('font', weight='bold')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"].encode("ascii","ignore") 

#Plotting the resultant vis_stacked image:

data = pf.getdata(path+'uvstacked.fits')
f = aplpy.FITSFigure(path+'uvstacked.fits')
f.show_grayscale(vmin=data[0,0,:,:].min(),vmax=data[0,0,:,:].max(), stretch='sqrt')
#f.show_contour(image, levels=10)
f.show_colorscale(vmin=data[0,0,:,:].min(),vmax=data[0,0,:,:].max(),stretch='linear',cmap='cubehelix')
f.add_colorbar(axis_label_text=r'\textbf{Flux density [Jy/beam]}',location='right' , pad = 0.2)
f.set_theme('publishable')
f.set_axis_labels_font(size=15, weight=15)
f.set_tick_color('black')
f.save(path+'uvstacked.pdf')

       
#Computing the MAD of a fits file image:

data = pf.getdata(path+'uvstacked.fits')
SEM = stats.sem(data,axis= None)
MAD = (SEM/1.4826)*1e6 # +-micro J
print "MAD vis_stack=  -+ ", MAD, "uJ"




#Plotting the resultant im_stacked image:

data = pf.getdata(path+'GroundT_stacked_header_ref-model.fits')
f = aplpy.FITSFigure(path+'Imstacked-model.fits')
f.show_grayscale(vmin=data.min(),vmax=data.max(), stretch='sqrt')
#f.show_contour(image, levels=10)
f.show_colorscale(vmin=data.min(),vmax=data.max(),stretch='linear',cmap='cubehelix')
f.add_colorbar(axis_label_text=r'\textbf{Flux density [Jy/beam]}',location='right' , pad = 0.2)
f.set_theme('publishable')
f.set_axis_labels_font(size=15, weight=15)
f.set_tick_color('black')
f.save(path+'GroundTruth_stacked.pdf')

       
#Computing the MAD of a fits file image:

data = pf.getdata(path+'GroundT_stacked_header_ref-model.fits')
SEM = stats.sem(data,axis= None)
MAD = (SEM/1.4826)*1e6 # +-micro J
print "MAD im_stack =  -+ ", MAD, "uJ"
