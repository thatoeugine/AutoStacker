import numpy as np
import cPickle as pickle
import os
import time
import pdb
import ConfigParser

#from emerlin_antenna import *
#from vla_antenna import *

def runSimulateData(config):#, parameter_filename):
  
  run_dir = os.getcwd()
  
  #cmd = 'casa --nogui --log2term -c {0}/simulatedata/simdata.py {0}/{1}'.format(run_dir, parameter_filename)
  #cmd = 'casa --nogui --log2term -c {0}/simulatedata/simdata.py {0}/{1}'.format(run_dir, 'temp_config.p')
  #print(cmd)
  os.chdir(config.get('pipeline', 'data_path'))
  temp_dir = os.getcwd()
  pickle.dump(config, open(temp_dir+'/temp_config.p', 'wb'))
  cmd = 'casa --nogui --log2term -c {0}/simulatedata/simdata.py {1}/{2}'.format(run_dir, temp_dir, 'temp_config.p')
  print(cmd)
  os.system(cmd)


if __name__=='__main__':
  
  #config = ConfigParser.ConfigParser()
  #config.read(sys.argv[-1])
  config = pickle.load(open(sys.argv[-1], 'rb'))
  #config.read('example.ini')
  n_ifs = config.getint('observation', 'n_IFs')
  bw = config.getfloat('observation', 'total_bandwidth')
  base_freq = config.getfloat('observation', 'lowest_frequency')
  channel_width = config.getfloat('observation', 'channel_width')
  n_chan = config.getint('observation', 'n_channels')
  channel_separation = config.getfloat('observation','channel_separation')
  msname = config.get('pipeline', 'project_name')+'.ms'
  #msname = config.get('pipeline', 'data_path')+msname
  imagename = msname+'.im'
  
  if_width = bw / n_ifs
  
  '''
  from simulatedata.fixImage import fixImage
  
  fixImage(datafile=config.get('field', 'fitsname'),
  				 centre_freq=base_freq,
  				 bandwidth=bw*n_chan*n_ifs,
  				 n_channels=n_chan*n_ifs,
  				 pixel_scale=config.getfloat('skymodel', 'pixel_scale'),
  				 Ra_central=config.get('field', 'field_ra'),
  				 Dec_central=config.get('field', 'field_dec'),
  				 input_casa_image=imagename)
  '''
  fitsimage = config.get('field', 'fitsname')
  if not config.getboolean('primarybeam', 'dopb'):
    importfits(fitsimage=fitsimage, imagename=imagename, overwrite=True)
  
  '''
  from vla_antenna import *
  sm.setspwindow(spwname = 'VLA-spw',
                 freq = '1.4GHz',
                 deltafreq = '3MHz',
                 freqresolution = '3MHz',
                 nchannels = 14,
                 stokes = 'LL RR')

  posvla = me.observatory('vla')
  sm.setconfig(telescopename = 'VLA-B',
               x = vla_xx,
               y = vla_yy,
               z = vla_zz,
               dishdiameter = diam.tolist(),
               mount = 'alt-az',
               coordsystem = 'local',
               referencelocation=posvla)
             
  sm.setfeed(mode='perfect R L')
  sm.setfield(sourcename = 'field1',
              sourcedirection = 'J2000 0h0m0.0s +90d0m0.0s')
  obs_date = '2011/03/21'
  ref_time = me.epoch('IAT',obs_date)
  sm.settimes(integrationtime = '5s',
              usehourangle = True,
              referencetime=ref_time)
            
  sm.observe('field1', 'VLA-spw',
             starttime = '0s', stoptime = '7200s')
  
  '''
####################################################################################################################    
#emerlin_lovell_config
  num_dishes = 7
  diameters=[25., 32., 25., 25., 25., 76., 25.]
  diameters_25=[25., 25., 25., 25., 25., 25., 25.]
  emerlin_diam = np.zeros((num_dishes,), np.float64) + diameters
  diam_25 = np.zeros((num_dishes,), np.float64) + diameters_25
  emerlin_xx=[3923069.1710, 3919982.7520, 3859711.5030, 3828714.5130, 3822473.3650, 3822252.6430, 3817176.5610]
  emerlin_yy=[146804.3680, -2651.9820, 201995.0770, 169458.9950, 153692.3180, 153995.6830, 162921.1790]
  emerlin_zz=[5009320.5280, 5013849.8260, 5056134.2510, 5080647.7490, 5085851.3030, 5086051.4430, 5089462.0570]
###################################################################################################################
    
    
    
    
    
  if (config.get('observation', 'uvcoverage')=='simulate'):
    print('making uv coverage!')
    sm.open(msname)
    for i in np.arange(1,n_ifs+1):
      fr = base_freq + if_width*(i-1)
      print(str(channel_width)+'Hz')
      sm.setspwindow(spwname = 'IF'+str(i),
                     freq = str(fr)+'Hz', # starting frequency
                     deltafreq = str(channel_separation)+'Hz', # increment per chan
                     freqresolution = str(channel_width)+'Hz', # width per chan
                     nchannels = n_chan,
                     #stokes = 'I Q U V')
                     )

    if config.get('observation', 'telescope')=='e-merlin':
      observatory = 'e-MERLIN'
      posemerlin = me.observatory(observatory)
      sm.setconfig(telescopename = 'e-MERLIN',
                   x = emerlin_xx,
                   y = emerlin_yy,
                   z = emerlin_zz,
                   dishdiameter = emerlin_diam.tolist(),
                   mount = 'alt-az',
                   coordsystem = 'global')
      sm.setfeed(mode='perfect R L')
    elif config.get('observation', 'telescope')=='jvla':
      posvla = me.observatory('vla')
      sm.setconfig(telescopename = 'VLA',
                   x = vla_xx,
                   y = vla_yy,
                   z = vla_zz,
                   dishdiameter = vla_diam.tolist(),
                   mount = 'alt-az',
                   coordsystem = 'global')
      sm.setfeed(mode='perfect R L')
    '''
    elif config.get('observation', 'telescope')=='both':
      sm.setconfig(telescopename = 'e-MERLIN+JVLA',
                   x = vla_xx+emerlin_xx,
                   y = vla_yy+emerlin_yy,
                   z = vla_zz+emerlin_zz,
                   dishdiameter = vla_diam.tolist()+emerlin_diam.tolist(),
                   mount = 'alt-az',
                   coordsystem = 'global')
      sm.setfeed(mode='perfect R L')
    '''
    source_dec_casa = config.get('field', 'field_dec').split(':')
    source_dec_casa = source_dec_casa[0]+'d'+source_dec_casa[1]+'m'+source_dec_casa[2]+'s'
    source_dirn = me.direction('J2000', config.get('field', 'field_ra'), source_dec_casa)
    sm.setfield(sourcename = config.get('field', 'name'),
                sourcedirection = source_dirn)
                #sourcedirection = 'J2000 10h30m0.0s +68d0m0.0s')
    
    obs_date = time.strftime('%Y/%m/%d', time.gmtime())
    ref_time = me.epoch('IAT', obs_date)
    sm.settimes(integrationtime = config.get('observation', 't_int')+'s',
                usehourangle = True,
                referencetime = ref_time)

    for i in np.arange(1,n_ifs+1):
      sm.observe(config.get('field', 'name'), 'IF'+str(i),
                 starttime = '0s',
                 stoptime = config.get('observation', 'observation_time')+'s')
  else:
    print('loading uv coverage!')
    msname = config.get('observation', 'uvcoverage_ms_file')
    sm.openfromms(msname)
  
  sm.predict(imagename=imagename)
  
  if config.get('observation', 'noisemode') == 'uniform':
    sm.setnoise(mode='simplenoise', simplenoise=config.get('observation', 'uniform_noise')+'Jy')
    sm.corrupt()
  elif config.get('observation', 'noisemode') == 'real':
    rms_noise_list = pickle.load(open(config.get('observation', 'noise_file'), 'rb'))
    for bl_spw in rms_noise_list:
      sm.openfromms(msname)
      #this_selection = {'baseline' : bl_spw[0], 'spw': bl_spw[1]}
      print bl_spw
      ant1 = bl_spw[0].split('&')[0]
      ant2 = bl_spw[0].split('&')[1]
      
      # CASA antenna naming schemes are horrible and inconsistent!
      # Convert from Antenna 'Name' to 'ID'
      ant1 = ant1.replace('1', '0')
      ant1 = ant1.replace('2', '1')
      ant1 = ant1.replace('5', '4')
      ant1 = ant1.replace('6', '5')
      ant1 = ant1.replace('7', '6')
      ant1 = ant1.replace('8', '7')
      ant1 = ant1.replace('9', '8')
      
      ant2 = ant2.replace('1', '0')
      ant2 = ant2.replace('2', '1')
      ant2 = ant2.replace('5', '4')
      ant2 = ant2.replace('6', '5')
      ant2 = ant2.replace('7', '6')
      ant2 = ant2.replace('8', '7')
      ant2 = ant2.replace('9', '8')
      
      this_selection = 'ANTENNA1=={0} and ANTENNA2=={1}'.format(ant1, ant2)
      this_rms = rms_noise_list[bl_spw]
      print this_selection
      print bl_spw[1]
      print str(this_rms)+'Jy'
      sm.setdata(spwid=int(bl_spw[1]), msselect=this_selection)
      sm.setnoise(mode='simplenoise', simplenoise=str(this_rms)+'Jy')
      sm.corrupt()
      sm.done()
  if config.get('imager', 'type')=='nwimager':
    exportuvfits(vis=msname,
               fitsfile=config.get('pipeline', 'project_name')+'.uvfits')
