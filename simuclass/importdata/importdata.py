'''
Script to import a SuperCLASS UV coverage into a CASA Measurement Set

This should be run *after* a set of UVCOP-ed individual IF fits
files have been exported from AIPS using exportdata.
'''
import ConfigParser
import pdb

from recipes.setOrder import setToCasaOrder

def runImportSuperCLASS(config, parameter_filename):
  
  run_dir = os.getcwd()
  
  cmd = 'casapy --nogui --log2term -c {0}/importdata/importdata.py {0}/{1}'.format(run_dir, parameter_filename)
  print(cmd)
  os.chdir(config.get('import', 'data_path'))
  os.system(cmd)

if __name__=='__main__':
  
  config = ConfigParser.ConfigParser()
  config.read(sys.argv[-1])
  
  if_list = config.get('import', 'if_list')
  if_list = if_list.split(',')

  concatlist = []

  for intfreq in if_list:
    if_name = config.get('import', 'name_root')+'.IF{0}'.format(intfreq)
    
    importuvfits(fitsfile=if_name,
                  #scanreindexgap_s=60,
                  vis=if_name+'_oneif.ms')

    setToCasaOrder(inputMS=if_name+'_oneif.ms',
                   outputMS=if_name+'_srt.ms')

    fixvis(vis=if_name+'_srt.ms',
           reuse=False,
           outputvis=if_name+'_srtfix.ms')

    concatlist.append(if_name+'_srtfix.ms')

  concat(vis=concatlist,
         concatvis=config.get('import', 'name_root')+'_allif.ms',
         timesort=True)

  os.system('rm -rf *_srtfix.ms *_srt.ms *_oneif.ms')
