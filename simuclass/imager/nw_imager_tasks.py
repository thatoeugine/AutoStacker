'''
Define some tasks for Nick Wrigley imager...
'''


def runfitld2(datain, indisk, thisdir, sourcename):
  fitld = AIPSTask('FITLD')
  fitld.datain = datain    # "datain" is correct here (use "indata" in most other places)
  fitld.outdisk = indisk
  fitld.digicor = -1
  fitld.douvcomp = -1
  fitld.clint = 8/60
  fitld.bif = 1
  fitld.eif = 0
  fitld.outdata = AIPSUVData(sourcename,'UVDATA',indisk,thisdir)
  fitld.go()


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


def removeall(option, target, indisk):
  if option == True:
    print 'removing all in aips pcat'
    temp_pca = AIPSCat(indisk)
    for fitsfil in temp_pca[indisk]:
      suboutname='SUB' + target[-4:-1]+target[-1]
      if (target in fitsfil.name) and ('UV' in fitsfil.type): #
        uvdata = AIPSUVData(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
        uvdata.zap()
      if ((target in fitsfil.name) or (suboutname in fitsfil.name)) and ('MA' in fitsfil.type):
        imdata = AIPSImage(fitsfil.name,fitsfil.klass,indisk,fitsfil.seq)
        imdata.zap()