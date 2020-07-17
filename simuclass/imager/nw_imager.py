import ConfigParser
import pdb

from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData

from nw_imager_tasks import *

if __name__=='__main__':
  config = ConfigParser.ConfigParser()
  config.read(sys.argv[-1])

  #Base Path is where we keep the eMERLIN_tasks directory
  sys.path.append(config.get('pipeline', 'base_path'))
  from eMERLIN_tasks.eMERLIN_tasks import *

  aips_no = config.getint('pipeline', 'aips_number')
  print aips_no
  AIPS.userno = aips_no
  indisk = config.getint('pipeline', 'aips_indisk')
  thisdir=0
  data_path = config.get('pipeline', 'data_path')
  data_file = data_path+config.get('pipeline', 'project_name')+'.uvfits'
  outfilename = config.get('pipeline', 'project_name')+'.fits'
  fovradius = config.getfloat('configuration', 'fov_radius')
  noiter= config.getint('configuration', 'number_of_iterations')
  target  = config.get('field', 'name')#'obsfield'
  option = config.getboolean('pipeline', 'clear_aips')

  print AIPSCat(indisk)
  print config.get('pipeline', 'project_name')+'.uvfits'

  #remove all files previously in aips catalogue
  removeall(option, target, indisk)

  # fitld
  runfitld2(data_file, indisk, thisdir, target)
  pca = AIPSCat(indisk)
  fitsfil = pca[indisk][-1]

  data = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)

  print 'fix antenna'
  # fix the antenna table and run indxr
  data.zap_table('AN', 1)
  runtbin(data, config.get('configuration', 'antenna_file'))
  #runextdest(data, 'AN', 1)
  #runtabedCopy(data, 'AN', 2, data, 1)
  runindxr(data)

  # run the imager
  print 'run chessboard chop'
  chessboard64wide(data, fovradius)  

  #clean the facets
  # now image all the uvfiles, base the imaging size larger than each square
  pca = AIPSCat(indisk)
  print 'run widemaps'
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('UV' in fitsfil.klass): #
      uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      #widemaps(uvdata,indisk,fovradius,target,noiter)
      widebandmaps(uvdata,indisk,fovradius,target,noiter)

  #'''
  #now use APCLN  added nwrigley 20141106
  pca = AIPSCat(indisk)
  for fitsfil in pca[indisk]:          
    if (target in fitsfil.name) and ('IIM' in fitsfil.klass):
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      apclean(imdata,target)
  #'''


  print 'run trimmer'
  # trim and flatten the images
  pca = AIPSCat(indisk)
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('ICL' in fitsfil.klass): #
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      #flat(imdata,fovradius,target)  
      trimmer(imdata,fovradius,target)

  print 'run flatn'
  runflatn(imdata, fovradius, target)

  #'''
  print 'run primary beam correction'
  #beam correction
  pca = AIPSCat(indisk)
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('FLATN' in fitsfil.klass):
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      print "Initialising primary beam correction.."
      pbcorr(imdata)
  #'''

  print 'export fits image'
  pca = AIPSCat(indisk)
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('FLATN' in fitsfil.klass):
      if os.path.isfile(os.path.join(data_path, outfilename)):
        os.system('rm -f '+os.path.join(data_path, outfilename+'.flatn'))
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      runfittp(imdata, data_path, outfilename+'.flatn')
    if (target in fitsfil.name) and ('PBCOR' in fitsfil.klass):
      if os.path.isfile(os.path.join(data_path, outfilename)):
        os.system('rm -f '+os.path.join(data_path, outfilename+'.pbcor'))
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      runfittp(imdata, data_path, outfilename+'.pbcor')

      
  '''
  print 'create catalogue with SAD'
  #Catalogue
  pca = AIPSCat(indisk)
  for fitsfil in pca[indisk]:
    if (target in fitsfil.name) and ('PBCOR' in fitsfil.klass):#('IM' in fitsfil.klass):
      imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
      SAD(imdata)
  '''

  print "Your AIPS catalog now looks like this:"
  print AIPSCat(indisk)
	

