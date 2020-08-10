"""
@author: Thato Manamela

Uinversity of Pretoria, Department of Physics, Asronomy Research Group.

"""
#____________________________________________________________________________________________

"""
# This module images the stacked visibilities/ generates the stacked image
"""
#--------------------------------------------------------------------------------------------




import casa
import numpy as np
import ConfigParser
import sys

if __name__=='__main__':
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[-1])
    
    path = config.get('pipeline', 'data_path')
    stampsize = config.getint('casa_imager', 'stampsize')
    cell_ = config.get('casa_imager', 'cell')
    vis_ = config.get('casa_imager', 'stack_msfile_name')

    # Image the uv-stacked data to produce an image.
    casa.clean(vis=path+vis_, imagename=path+'uvstacked',cell = cell_,\
               imsize = stampsize, mask = [int(stampsize/2)-2, int(stampsize/2)-2,
                  int(stampsize/2)+2, int(stampsize/2)+2])

    casa.exportfits(path+'uvstacked.image', path+'uvstacked.fits') #exporting casa image to fts
