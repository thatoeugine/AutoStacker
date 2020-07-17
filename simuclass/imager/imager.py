import ConfigParser
import pdb
import os
from astropy.io import fits

def runNWImager(config, parameter_filename):
  
  run_dir = os.getcwd()
  
  cmd = 'source /aips/LOGIN.SH && /usr/local/bin/ParselTongue.old {0}/imager/nw_imager.py {0}/{1}'.format(run_dir,
                                                           parameter_filename)
  print(cmd)
  os.system(cmd)
  
def runCASAClean(config, parameter_filename):
  
  run_dir = os.getcwd()
  
  cmd = 'casapy --nogui --log2term -c {0}/imager/casa_imager.py {0}/{1}'.format(run_dir, parameter_filename)
  print(cmd)
  os.chdir(config.get('pipeline', 'data_path'))
  os.system(cmd)
  
def runWSClean(config):
  
  run_dir = os.getcwd()
  
  #source_cmd = ['bash', '-c', 'source /usr/local/lofar-2.12/init_env_release.sh & source /usr/local/lofar-2.12/lofar/release/lofarinit.sh']
  #source_cmd_2 = ['bash', '-c', 'source /usr/local/lofar-2.12/lofar/release/lofarinit.sh']
  #proc = subprocess.Popen(source_cmd, stdout=subprocess.PIPE)
  #proc = subprocess.Popen(source_cmd_2, stdout=subprocess.PIPE)
  
  input_image_hdr = fits.getheader(config.get('pipeline', 'data_path')+config.get('field', 'fitsname'))
  naxis1 = input_image_hdr['NAXIS1']
  naxis2 =  input_image_hdr['NAXIS2']
  #source_cmd = 'source /usr/local/lofar-2.12/lofar/release/lofarinit.sh; source /usr/local/lofar-2.12/init_env_release.sh;'
  source_cmd = 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/gcc-4.9.3/lib64:/usr/local/gcc-4.9.3/lib/;'
  
  if config.get('imager', 'algo')=='cs':
    cmd = 'wsclean -j {0} \
                   -name {1} \
                   -size {2} {3} \
                   -scale {4}asec \
                   -niter {5} \
                   -threshold {6} \
                   -weight {7} \
                   -mgain 0.8 \
                   -multiscale \
                   {8} '.format(
                  config.getint('imager', 'n_cores'),
                  config.get('pipeline', 'data_path')+\
                    config.get('pipeline', 'project_name')+'.wsclean',
                  naxis1,
                  naxis2,
                  config.getfloat('skymodel', 'pixel_scale'),
                  config.getint('imager', 'number_of_iterations'),
                  config.getfloat('imager', 'threshold'),
                  config.get('imager', 'weight'),
                  config.get('pipeline', 'project_name')+'.ms'
                  )
  else:
    cmd = 'wsclean -j {0} \
                   -name {1} \
                   -size {2} {3} \
                   -scale {4}asec \
                   -niter {5} \
                   -threshold {6} \
                   -weight {7} \
                   -nonegative \
                   {8} '.format(
                  config.getint('imager', 'n_cores'),
                  config.get('pipeline', 'data_path')+\
                    config.get('pipeline', 'project_name')+'.wsclean',
                  naxis1,
                  naxis2,
                  config.getfloat('skymodel', 'pixel_scale'),
                  config.getint('imager', 'number_of_iterations'),
                  config.getfloat('imager', 'threshold'),
                  config.get('imager', 'weight'),
                  config.get('pipeline', 'project_name')+'.ms'
                  )
  os.chdir(config.get('pipeline', 'data_path'))
  print(cmd)
  os.system(source_cmd+cmd)
  os.chdir(run_dir)

if __name__=='__main__':
  
  from AIPS import AIPS, AIPSDisk
  from AIPSTask import AIPSTask, AIPSList
  from AIPSData import AIPSUVData, AIPSImage, AIPSCat
  from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData

  from eMERLIN_tasks import *

  config = ConfigParser.ConfigParser()
  config.read(sys.argv[-1])
  aips_no = config.getint('pipeline', 'aips_number')
  AIPS.userno = aips_no
  indisk = config.getint('pipeline', 'aips_indisk')
  thisdir=0

  print AIPSCat(indisk)

  # fitld
  runfitld(config.get('pipeline', 'project_name')+'.uvfits', indisk, thisdir)
  #pdb.set_trace()
  pca = AIPSCat(indisk)
  fitsfil = pca[indisk][-1]

  # fix the antenna table and run indxr
  runtbin(data, config.get('configuration', 'antenna_file'))
  runextdest(data, 'AN', 1)
  runtabed(data, 'AN', 2, data, 1)
  runindxr(data)

  # run the imager
  chesboard64wide(uvdata, fovradius)
  
  for fitsfil in pca[indisk]:
    if (target in fitsil.name) and ('UV' in fitsfil.klass):
      uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      widemaps(uvdata,indisk,fovradius,target,noiter)
  
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('ICL' in fitsfil.klass):
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      flat(imdata,fovradius,target)  

  runflatn(imdata, fovradius, target)

def runflatn(imdata, fovradius, widetarget):
    print 'Flattening facets'
    
    cellpitch = (imdata.header.bmin * 3600) / 3
    
    flatn=AIPSTask('FLATN')
    suboutname='SUB' + widetarget[-4:-1]+widetarget[-1]
    flatn.inname=suboutname
    flatn.inseq=0
    flatn.outname=widetarget
    flatn.outclass='FLATN'
    flatn.imsize[1]=256+(fovradius*60*2/cellpitch)
    flatn.imsize[2]=256+(fovradius*60*2/cellpitch)
    flatn.nfield=64
    flatn.inseq=0
    flatn.nmaps=1
    flatn.go()

def runextdest(indata, inext, invers):
  extdest = AIPSTask('EXTDEST')
  extdest.indata = indata
  extdest.inext = inext
  extdest.invers = invers
  extdest.go()
  
def runtabedCopy(indata, inext, invers, outdata, outvers):
  tabed = AIPSTask('TABED')
  tabed.indata = indata
  tabed.inext = inext
  tabed.invers = invers
  tabed.outdata = outdata
  tabed.outvers = outvers
  tabed.optype = 'COPY'
  tabed.go()
  
def runtbin(outdata, infile):
  tbin = AIPSTask('TBIN')
  tbin.outdata = outdata
  tbin.intext = infile
  tbin.go()
  
def runuvconfrompoint(antenna_file, outname,
                     start_freq, start_wavelength,
                     dec, min_ha, max_ha,
                     t_int, channel_inc, n_chan,
                     blockage,
                     w_term=0):
  uvcon = AIPSTask('UVCON')
  uvcon.infile = antenna_file
  #uvcon.in2name = fitsfil.name
  #uvcon.in2class = fitsfil.klass
  #uvcon.in2seq = fitsfil.seq
  #uvcon.in2disk = 1
  uvcon.outname = outname
  #uvcon.nmaps = 1
  #uvcon.cmethod = 'DFT'
  #uvcon.cmodel = 'COMP'
  uvcon.smodel[1:] = 1, 0, 0, 0, 0, 0, 0
  uvcon.aparm[1:] = start_freq, start_wavelength, \
                    dec, min_ha, max_ha, 0,\
                    t_int, channel_inc, n_chan, \
                    blockage
  uvcon.bparm[1:] = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  uvcon.do3dimag = w_term
  #uvcon.inp()
  uvcon.go()

def runuvconfromfits(fitsfil, antenna_file, outname,
                     start_freq, start_wavelength,
                     dec, min_ha, max_ha,
                     t_int, channel_inc, n_chan,
                     blockage,
                     w_term=0):
  uvcon = AIPSTask('UVCON')
  uvcon.infile = antenna_file
  uvcon.in2name = fitsfil.name
  uvcon.in2class = fitsfil.klass
  uvcon.in2seq = fitsfil.seq
  uvcon.in2disk = 1
  uvcon.outname = outname
  uvcon.nmaps = 1
  uvcon.cmethod = 'DFT'
  uvcon.cmodel = 'IMAG'
  uvcon.aparm[1:] = start_freq, start_wavelength, \
                    dec, min_ha, max_ha, 0, \
                    t_int, channel_inc, n_chan, \
                    blockage
  uvcon.bparm[1:] = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  uvcon.do3dimag = w_term
  uvcon.inp()
  uvcon.go()

def runaxdefine(fitsfil, naxis, axtype, axinc, axval, axref):
  
  axdefine = AIPSTask('AXDEFINE')
  axdefine.inname = fitsfil.name
  axdefine.inclass = fitsfil.klass
  axdefine.inseq = int(fitsfil.seq)
  axdefine.indisk = int(1)
  axdefine.naxis = naxis
  axdefine.axtype = axtype
  axdefine.axinc = axinc
  axdefine.axval[1:] = axval, axval
  axdefine.axref = int(axref)
  axdefine.inp()
  axdefine.go()
  
def addFreqStokesAxes(fitsfil,
                      ra, dec,
                      pix_scale,
                      min_freq,
                      bandwidth):
  '''Give correct axes definitions to input fits file image.
  NOTE: to alter, require aips wizardy
  uv = wizuvdata
  uv.header.append(stuff)
  uv.update()
  '''
  runaxdefine(fitsfil, 1, 'RA---SIN', pix_scale, ra, 1)
  runaxdefine(fitsfil, 1, 'DEC--SIN', pix_scale, dec, 1)
  runaxdefine(fitsfil, 1, 'FREQ', min_freq, bandwidth, 1)
  runaxdefine(fitsfil, 1, 'STOKES', 1, 1, 1)
  

