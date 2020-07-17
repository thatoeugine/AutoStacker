#Modules modifies the simuclass simulates image from 4D to 2D for PASTA  friendly.

import sys
import astropy.io.fits as pf
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

data_path = config.get('pipeline', 'data_path')
all_gals_fname = data_path+config.get('field', 'fitsname')


data_, header_ = pf.getdata(all_gals_fname, header = True)

del header_["NAXIS3"]
del header_["NAXIS4"]

del header_["CTYPE3"]
del header_["CTYPE4"]

del header_["CRVAL3"]
del header_["CRVAL4"]

del header_["CDELT3"]
del header_["CDELT4"]

del header_["CRPIX3"]
del header_["CRPIX4"]

del header_["CUNIT3"]
del header_["CUNIT4"]

del header_["PC01_03"]

del header_["PC01_04"]

del header_["PC02_04"]

del header_["PC02_03"]

del header_["PC03_01"]
del header_["PC03_02"]
del header_["PC03_03"]
del header_["PC03_04"]

del header_["PC04_01"]
del header_["PC04_02"]
del header_["PC04_03"]
del header_["PC04_04"]



data_2 = data_[0,0,:,:]


pf.writeto(filename = data_path+"2Dimage.fits", data = data_2, header =  header_, overwrite = True)

