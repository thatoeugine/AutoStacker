'''
Script to export a SuperCLASS UV coverage for import to
a CASA Measurement Set.

UVCOPs and FITTPs individual IFs for later reading in to CASA.
'''
import ConfigParser
import pdb

from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData

from eMERLIN_tasks import *

def runuvcop_per_if(uvdata, bif, eif, indisk):
  outname = uvdata.name + '.B{0}E{1}'.format(bif, eif)

  uvcop = AIPSTask('UVCOP')
  uvcop.indata = uvdata
  uvcop.outname = outname
  uvcop.outdisk = indisk
  uvcop.bif = int(bif)
  uvcop.eif = int(eif)
  uvcop.go()
  return outname

def runExportSuperCLASS(config, parameter_filename):
  
  run_dir = os.getcwd()
  
  cmd = 'parseltongue {0}/exportdata/exportdata.py {0}/{1}'.format(run_dir,
                                                                   parameter_filename)
  print(cmd)
  os.system(cmd)

if __name__=='__main__':
  config = ConfigParser.ConfigParser()
  config.read(sys.argv[-1])
  
  aips_no = config.getint('export', 'aips_number')
  AIPS.userno = aips_no
  indisk = config.getint('export', 'aips_indisk')
  sourcename = config.get('export', 'source_name')

  pca = AIPSCat(indisk)
  plist = pca[indisk]
  source_files = (item for item in plist if item['name'] == sourcename)
  fitsfil = source_files.next()
  
  while fitsfil.klass != config.get('export', 'source_class'):
    fistfil = (item for item in plist if item['name'] == sourcename).next()

  uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)

  if_list = config.get('export', 'if_list')
  if_list = if_list.split(',')
  
  #pdb.set_trace()

  for intfreq in if_list:
    outname = runuvcop_per_if(uvdata, intfreq, intfreq, indisk)
    pca = AIPSCat(indisk)
    plist = pca[indisk]
    fitsfil = (item for item in plist if item['name'] == outname).next()
    
    oneif_data = AIPSUVData(fitsfil.name, fitsfil.klass, indisk, fitsfil.seq)
    runfittp(oneif_data, config.get('export', 'output_dir'), config.get('export', 'name_root')+'.IF{0}'.format(intfreq))

