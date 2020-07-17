# Define the AIPS tasks to use in the e-MERLIN pipeline.
# mkargo 2011
# Last updated 20140814

import os, sys, math, numpy, re
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData


def runfitld(datain, indisk, thisdir):
	fitld = AIPSTask('FITLD')
	fitld.datain = datain		# "datain" is correct here (use "indata" in most other places)
	fitld.outdisk = indisk
	fitld.digicor = -1
	fitld.douvcomp = -1
	fitld.clint = 8/60
	fitld.bif = 1
	fitld.eif = 0
	fitld.outdata = AIPSUVData('TMP','UVDATA',indisk,thisdir)
	fitld.go()


def runmsort(indata):
	print 'Running MSORT.'
	msort = AIPSTask('MSORT')
	msort.indata = indata
	msort.outdata = indata
	msort.sort = 'TB'
	msort.go()

def runindxr(indata):
	print 'Running INDXR.'
	indxr = AIPSTask('INDXR')
	indxr.cparm[3] = 8/60
	indxr.indata = indata
	indxr.go()

def runfring(indata, calsour, timer, docalib, gainuse, flagver, doband, bpver, refant, refantlist, solint, aparm, dparm, snver, snr, bchan, echan):
	print 'Running FRING on ', format(calsour)
	fring = AIPSTask('FRING')
	fring.indata = indata
	fring.calsour[1:] = calsour
	fring.timer[1:] = timer
	fring.docalib = docalib
	fring.gainuse = gainuse
	fring.flagver = flagver
	fring.doband = doband
	fring.bpver = bpver
	fring.refant = refant
	fring.search[1:] = refantlist
	fring.solint = solint
	fring.cmethod = 'DFT'
	fring.aparm[1:] = aparm
	fring.aparm[7] = snr
	fring.dparm[1:] = dparm
	fring.snver = snver
	fring.bchan = bchan
	fring.echan = echan
	fring.go()

def runsnplt(indata, inext, invers, sources, nplots, optype, dotv):
	snplt = AIPSTask('SNPLT')
	snplt.indata = indata
	snplt.inext= inext
	snplt.invers = invers
	snplt.sources[1:] = sources
	snplt.nplots = nplots
	snplt.optype = optype
	snplt.dotv = dotv
	snplt.go()

def runlwpla(indata, outfile, qa):
	plver=0
	outfile = re.sub(r"\s+", '-', outfile)
	for tab in indata.tables :
		if 'PL' in tab[1] :
			plver = plver + 1
	if os.path.exists(outfile) :
		os.remove(outfile)
	lwpla = AIPSTask('LWPLA')
	lwpla.indata = indata
	lwpla.plver = 1
	lwpla.inver = plver
	lwpla.outfile = outfile
	lwpla.msgkill = -10
	if qa == 1:
		lwpla.dparm[5] = 2
		lwpla.docolor = 1
		lwpla.plcolors[10][1:] = 1,1,1
	lwpla.go()
	# Remove PL tables
	indata.zap_table('PL', -1)

def runclcal(indata, sources, calsour, opcode, interpol, snver, inver, gainver, gainuse, refant):
	clcal = AIPSTask('CLCAL')
	clcal.indata = indata
	clcal.sources[1:] = sources
	clcal.calsour[1:] = calsour
	clcal.opcode = opcode
	clcal.interpol = interpol
	clcal.snver = snver
	clcal.inver = inver
	clcal.gainver = gainver
	clcal.gainuse = gainuse
	clcal.refant = refant
	clcal.go()


def runpossm(indata, sources, timer, anten, basel, aparm, bparm, bif, eif, bchan, echan, docalib, gainuse, flagver, stokes, doband, bpver, codetype, solint, nplots, dotv, freqid):
	possm = AIPSTask('POSSM')
	possm.indata = indata
	possm.sources[1:] = sources
	possm.timer[1:] = timer
	possm.aparm[1:] = aparm
	possm.bparm[1:] = bparm
	i = 1
	for ant in anten:
		possm.antennas[i] = int(ant)
		i = i + 1
	i = 1
	for ant in basel:
		possm.baseline[i] = int(ant)
		i = i + 1
	possm.bif = bif
	possm.eif = eif
	possm.bchan = bchan
	possm.echan = echan
	possm.docalib = docalib
	possm.gainuse = gainuse
	possm.flagver = flagver
	possm.stokes = stokes
	possm.doband = doband
	possm.bpver = bpver
	possm.codetype = codetype
	possm.solint = solint
	possm.nplots = nplots
	possm.dotv = dotv
#	possm.freqid = 1
	possm.go()

def runquickpossm(indata, source):
	possm = AIPSTask('POSSM')
	possm.indata = indata
	possm.source[1:] = source
	possm.aparm[9] = 1
	possm.nplots = 6
	possm.solint = 30 # plot every 30min
	possm.stokes = 'HALF'
	possm.flagver = 0
	possm.go()


def runbpass(indata, calsour, refant, bpassprm, soltyp, solint, docal, specindx):
#	bpass = AIPSTask('BPASS', version='OLD')
	bpass = AIPSTask('BPASS')
#	bpass.version='OLD'	# Only if on the JODRELL system!
#	bpass.uvrange[1:] = 0,0
	bpass.indata = indata
	bpass.calsour[1:] = calsour
	bpass.refant = refant
	bpass.bpassprm[1:] = bpassprm
#	bpass.bif = bif
#	bpass.eif = eif
	bpass.soltype = soltyp
	bpass.solint = solint
	bpass.docal = docal
	bpass.specindx = specindx
	bpass.go()

def runsetjy(indata, sources, bif, eif, zerosp, optype):
	setjy = AIPSTask('SETJY')
	setjy.indata = indata
	setjy.sources[1:] = sources
	setjy.bif = bif
	setjy.eif = eif
	setjy.zerosp[1:] = zerosp, 0, 0, 0
	setjy.optype = optype
	setjy.inp()
	setjy.go()

def runcalib(indata, calsour, timer, uvrang, docalib, gainuse, flagver, doband, bpver, cmethod, refant, solint, aparm, doflag, soltype, solmode, minamper, minphser, cparm, snver, antwt, weightit):
	calib = AIPSTask('CALIB')
	calib.indata = indata
	calib.calsour[1:] = calsour
	calib.timer[1:] = timer
	calib.uvrange[1:] = uvrang
	calib.docalib = docalib
	calib.gainuse = gainuse
	calib.flagver = flagver
	calib.doband = doband
	calib.bpver = bpver
	calib.cmethod = cmethod
	calib.refant = refant
	calib.solint = solint
	calib.cparm[1:] = cparm
	if(solmode == 'P'):
		calib.aparm[1] = 3
	else:
		calib.aparm[1] = 4
	calib.aparm[6] = 3
	calib.aparm[7] = 5
	calib.cparm[1] = 10
	calib.doflag = doflag
	calib.soltype = soltype
	calib.solmode = solmode
	calib.minamper = minamper
	calib.minphser = minphser
	calib.snver = snver
	calib.antwt[1:] = antwt
	calib.weightit = weightit
	calib.go()


def runselfcalib(indata, in2data, calsour, timer, uvrang, docalib, gainuse, flagver, doband, bpver, cmethod, refant, solint, aparm, doflag, soltype, solmode, minamper, minphser, cparm, snver, antwt, weightit):
	# CALIB procedure for self-calibration
	calib = AIPSTask('CALIB')
	calib.indata = indata
	calib.in2data = in2data
	calib.calsour[1:] = calsour
	calib.timer[1:] = timer[1:]
	calib.uvrange[1:] = uvrang[1:]
	calib.docalib = docalib
	calib.gainuse = gainuse
	calib.flagver = flagver
	calib.doband = doband
	calib.bpver = bpver
	calib.cmethod = cmethod
	calib.refant = refant
	calib.solint = solint
#	calib.aparm[1:] = aparm[1:]
	if(solmode == 'P'):
		calib.aparm[1] = 3
	else:
		calib.aparm[1] = 4
	calib.aparm[6] = 3
	calib.aparm[7] = 5
	calib.cparm[1] = 10
	calib.doflag = doflag
	calib.soltype = soltype
	calib.solmode = solmode
	calib.minamper = minamper
	calib.minphser = minphser
	calib.cparm[1:] = cparm[1:]
	calib.snver = snver
	calib.antwt[1:] = antwt[1:]
	calib.weightit = weightit
	calib.go()



def runuvplt(indata, sources, stokes, timer, anten, basel, bif, eif, docalib, gainuse, flagver, doband, bpver, aparm, bparm, doweight, refant, dotv):
	uvplt = AIPSTask('UVPLT')
	uvplt.indata = indata
	uvplt.sources = sources
	uvplt.stokes = stokes
	uvplt.timerang = timer
	uvplt.antennas[1:] = anten
	uvplt.baseline[1:] = basel
	uvplt.bif = bif
	uvplt.eif = eif
	uvplt.docalib = docalib
	uvplt.gainuse = gainuse
	uvplt.flagver = flagver
	uvplt.doband = doband
	uvplt.bpver = bpver
	uvplt.aparm[1:] = aparm[1:]
	uvplt.bparm[1:] = bparm[1:]
	uvplt.doweight = doweight
	uvplt.refant = refant
	uvplt.dotv = dotv
	uvplt.go()



def runvplot(indata, sources, stokes, timer, anten, basel, bif, eif, bchan, echan, docalib, gainuse, flagver, doband, bpver, aparm, bparm, refant, dotv, nplots):
	vplot = AIPSTask('VPLOT')
	vplot.indata = indata
	vplot.sources[1:] = sources
	vplot.stokes = stokes
	vplot.timerang[1:] = timer
	vplot.antennas[1:] = anten
	vplot.baseline[1:] = basel
	vplot.bchan = bchan
	vplot.echan = echan
	vplot.avgchan = 1
	vplot.bif = bif
	vplot.eif = eif
	vplot.docalib = docalib
	vplot.gainuse = gainuse
	vplot.flagver = flagver
	vplot.doband = doband
	vplot.bpver = bpver
	vplot.aparm[1:] = aparm
	vplot.bparm[1:] = bparm
	vplot.crowded = 0
	vplot.refant = refant
	vplot.dotv = dotv
	vplot.xinc = 1
	vplot.nplots = nplots
	vplot.go()

def runquickvplot(uvdata):
	vplot = AIPSTask('VPLOT')
	vplot.indata = uvdata
	vplot.bchan = 3*uvdata.header['naxis'][2] /8 
	vplot.bchan = 5*uvdata.header['naxis'][2] /8
	vplot.solint = 0.5
	vplot.crowded = 3
	vplot.do3col = 1
	vplot.bparm[1] = 0
	vplot.bparm[2]=-2
	vplot.nplot = 6
	vplot.go()


def rungetjy(indata, sources, calsour, bif, eif, snver):
	getjy = AIPSTask('GETJY')
	getjy.indata = indata
	getjy.sources[1:] = sources
	getjy.calsour[1:] = calsour
	getjy.bif = bif
	getjy.eif = eif
	getjy.snver = snver
	getjy.go()


def runimagr(indata, sources, docalib, gainuse, flagver, doband, bpver, bchan, echan, nchav, chinc, cellsiz, imsiz, niter, dotv, outdisk):
	imagr = AIPSTask('IMAGR')
	imagr.indata = indata
	imagr.sources[1:] = sources
	source = str(sources[0])
	if len(source)>12 :
		source = source[0:12]
	imagr.outname = source
	imagr.outdisk = outdisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1:] = cellsiz
	imagr.imsize[1:] = imsiz
	imagr.niter = niter
	imagr.dotv = dotv
	imagr.go()


def runkntr(image, factor, userno):
	imean = AIPSTask('IMEAN')
	imean.userid = userno
	imean.indata = image
	imean()
	kntr = AIPSTask('KNTR')
	kntr.dogrey = -1
	kntr.dovect = -1
	kntr.docont = 1
	kntr.indata = image
	kntr.blc[1:] = 0,0
	kntr.trc[1:] = 0,0
#	kntr.clev = 3*imean.pixstd
#	kntr.levs[1:] = -1,1,2,4,6,8,16,32,64,128
	kntr.cbpl = 16
	kntr.dotv = -1
	kntr.plev = 1
	rms = imean.pixstd
	peak = image.header.datamax
	# Set the kntr levels based on dynamic range
	firstlev = factor*abs(rms/peak)*300.
	kntr.levs[1] = firstlev * (-1.)
	kntr.levs[2] = firstlev
	i = 3
	while (kntr.levs[i-1]*2. < 100. and i<=30):
		kntr.levs[i] = 2 * kntr.levs[i - 1]
		i += 1
	kntr.go()


def aipsuvname(aipsdata):
	return aipsdata.name + '.' + aipsdata.klass + '.' + str(aipsdata.seq)


def runfittp(uvdata, fittpdir, fittpfile):
	fittp = AIPSTask('FITTP')
	fittp.indata = uvdata
	fittp.doall = -1
	fittp.intype = ''
	fittp.outtape = 1
	srcname = re.sub(r"\s+", '-', fittpfile)
	fittp.dataout = os.path.join(fittpdir, srcname)
	print "Saving to disk: ", uvdata.name, uvdata.klass, uvdata.seq
	fittp.go()


def runtasav(uvdata, fittpdir, fittpfile, indisk):
	# Backup the FG tables after flagging
	tasav = AIPSTask('TASAV')
	tasavfil = AIPSUVData(uvdata.name, 'TASAV', indisk, uvdata.seq)
	tasav.indata = uvdata
	tasav.outdata = tasavfil
	tasav.outseq = uvdata.seq
	print "Saving tables for ", uvdata.name, uvdata.klass, uvdata.seq
	tasav.go()
	runfittp(tasavfil, fittpdir, fittpfile)
	tasavfil.zap()

def runsousp(indata, sources):
	sousp = AIPSTask('SOUSP')
	sousp.indata = indata
	sousp.sources[1:] = sources
	sousp.order = 1
	sousp.dotv = -1
	sousp.go()
	return sousp.specindx

def runuvfix(datafile, inname, indisk, thisdir):
	uvfix = AIPSTask('UVFIX')
	klass = datafile.klass
	uvfix.indata = datafile
	uvfix.outdata = AIPSUVData(datafile.name,'UVFIX',indisk,thisdir)
	uvfix.go()
	datafile.zap()
	datafile = AIPSUVData(inname,'UVFIX',indisk,thisdir)
	datafile.rename(inname, klass, thisdir)

def runuvsrt(datafile, inname, indisk, thisdir):
	uvsrt = AIPSTask('UVSRT')
	klass = datafile.klass
	uvsrt.indata = datafile
	uvsrt.outdata = AIPSUVData(datafile.name,'UVSRT',indisk,datafile.seq)
	uvsrt.sort = 'TB'
	uvsrt.go()
	datafile.zap()
	datafile = AIPSUVData(inname,'UVSRT',indisk,thisdir)
	datafile.rename(inname,klass,thisdir)

def flagLOMK(datafile):
	# WRITE A FLAG FOR MK2-LT SPACING
	uvflg = AIPSTask('UVFLG')
	print "UVFLG: FLAGGING MK2-LT Baseline " + datafile.name + '.' + datafile.klass + '.' + format(datafile.seq)
	uvflg.indata = datafile
	uvflg.outfgver = 1
	uvflg.bif = 1
	uvflg.eif = 0
	uvflg.opcode = 'FLAG'
	uvflg.reason = 'Lovell - Mk2 baseline'
	uvflg.antennas[1] = 1 #Lovell
	uvflg.baseline[1] = 2 #MK2
	uvflg.go()

def runsplat(uvdata, outchan, tint, sbandl, sbandu, smootha, smoothb, smoothc, indisk):
	splat = AIPSTask('SPLAT')
	splat.indata = uvdata
	splat.outname = uvdata.name
	splat.outdisk = indisk
	splat.outclass = 'SPLAT'
	splat.outseq = uvdata.seq
	splat.solint = 0
	if outchan != -1 :
		aparm1 = 3
		splat.channel = uvdata.header['naxis'][2] / outchan
		splat.chinc = uvdata.header['naxis'][2] / outchan
	if tint != -1 :
		splat.solint = tint / 60
	splat.bif = sbandl # SBLANL -- Lower sub-band for extraction
	splat.eif = sbandu # SBANDU - Upper sub-band for extraction
	splat.douvcomp = -1
	splat.aparm[1] = 3
	splat.aparm[2] = 1
	# Smoothing options!?!?
	splat.smooth[1] = smootha # 1 = hanning default
	splat.smooth[2] = smoothb # 4 = hanning default
	splat.smooth[3] = smoothc # 1 = hanning default
	splat.flagver = 0
	splat.go()

def runprtmsg(prtask, outprint):
	prtmsg = AIPSTask('PRTMSG')
	clrmsg = AIPSTask('CLRMSG')
	prtmsg.prtask = prtask
	prtmsg.docrt = -1
	prtmsg.outprint = outprint
	# RUN VERB
	prtmsg.go()
	clrmsg()

def runuvhgm(indata, sources, anten, basel, SEFDbif, SEFDeif, SEFDbchan, SEFDechan):
	uvhgm = AIPSTask('UVHGM')
	uvhgm.indata = indata
	uvhgm.sources[1:] = sources
	uvhgm.antennas[1:] = anten
	uvhgm.baseline[1:] = basel
	uvhgm.bif = SEFDbif
	uvhgm.eif = SEFDeif
	uvhgm.flagver = 0
	uvhgm.stokes ='HALF'
	uvhgm.bchan = SEFDbchan
	uvhgm.echan = SEFDechan
	uvhgm.doall = 1
	uvhgm.axtype ='H'
	uvhgm.timerang[1:] = 0, 0, 0, 0, 0, 0, 0, 0
#	uvhgm.timerang[1:] = 0, 0, 0, 0, 3, 0, 0, 0
#	uvhgm.pixrange[1:] = -20, 20
	uvhgm.dotv = -1
	uvhgm.docal = 1
	uvhgm.gainuse = 0
	uvhgm.pixrange[1:] = -10, 10
	uvhgm.doband = -1
	uvhgm.nboxes = 1000
	uvhgm.go()

def runlistr(indata,sources,optype,docrt,outprint):
	listr = AIPSTask('LISTR')
	listr.indata = indata
	listr.sources[1:] = sources
	listr.optype = optype
	listr.docrt = docrt
	listr.outprint = outprint
	listr.go()



################################################################################
# Other (non-AIPS) tasks
################################################################################

def get_tab(uvdata, table):
	# find the number of tables of a certain type
	ver = 0
	for i in range(len(uvdata.tables)) :
		if table in uvdata.tables[i][1] :
			ver = uvdata.tables[i][0]
	print "HIGHEST TABLE OF TYPE", table, "is", ver
	return ver


def get_ant_num(uvdata, refant_name):
	# convert antenna name to number
	antab = uvdata.table('AN',1)
	for row in antab :
		if refant_name in row.anname :
			return row.nosta


def set_refant(uvdata, refant):
	# Check the reference antenna is present, set to sensible default if not.
	# Returns referance antenna name, refant number, and a list of antennas present in the data.
	searchants = ['Mk2', 'Pi', 'Da', 'Kn', 'Cm', 'Lo', 'De']
	antab = uvdata.table('AN',1)
	refantn = 0
	refantn = get_ant_num(uvdata, refant)
	if not refantn :
		print "Warning: No refant specified or requested reference antenna not present."
		for item in searchants:
			if refantn:
				break
			for row in antab :
				if item in row.anname :
					print "Warning: Using reference antenna", item
					refant = item
					refantn = get_ant_num(uvdata, refant)
					break

	refantlist = []
	for antenna in searchants:
		antnum = get_ant_num(uvdata, antenna)
		if antnum:
			refantlist.append(antnum)

	return (refant, refantn, refantlist)


def nodot(item):
	"""Filter out hidden files."""
	return item[0] != '.'



def time2hms(seconds):
	# Function to convert seconds to hh:mm:ss.ss format, returns a string
	h=int(seconds/3600)
	m=int(seconds % 3600)/60
	s=seconds-(h*3600)-(m*60)
	h=`h`
	m=`m`
	s="%4.2f" % s
	hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
	return hms




def inputinput(prompt):
	retries=1
	complaint='Yes or no, please!'
	while True :
		ok = raw_input(prompt)
		if ok == '' :
			print "You wanted to continue but were too lazy to say so."
			break
		if ok in ('y', 'ye', 'yes','Y'):
			print "Thank you, continuing."
			break
		if ok in ('n', 'N', 'no', 'nop', 'nope'):
			print "OK, quitting."
			sys.exit()
		retries = retries - 1
		if retries < 0:
			print "Don't be daft."
			sys.exit()
		print complaint
	return


def findmaxb(uvdata):
	maxbaseline = 0
	antab = uvdata.table('AN',1)
	for row in antab :
		for row2 in antab :
			xsep = row.stabxyz[0] - row2.stabxyz[0]
			ysep = row.stabxyz[1] - row2.stabxyz[1]
			zsep = row.stabxyz[2] - row2.stabxyz[2]
			hypxy = math.sqrt((xsep * xsep) + (ysep * ysep))
			hypxyz = math.sqrt((zsep * zsep) + (hypxy * hypxy))
			if hypxyz > maxbaseline :
				maxbaseline = hypxyz
	cellsize = (1.22 * (300000000 / uvdata.header.crval[2]) / maxbaseline) / 3.141592 * 180 * 3600 / 5
	print "maxbaseline = ", maxbaseline, "cellsize = ", cellsize
	return cellsize,cellsize


########################################################################################################################################
# Flux-calibration calculations
########################################################################################################################################


def dfluxpy(freq,uvdata):

	"Function to calculate 3C286 values, modified from Danielle's dfluxpy.py."
	"Updated 20140506 to include new coefficients and a calculation for projected"
	"baseline length from inspecting the u's and v's."


	data=WizAIPSUVData(uvdata.name, uvdata.klass, uvdata.disk, uvdata.seq)

	antab = data.table('AN',1)
	numLO = 10
	numM2 = 10
	for row in antab:			# Ignore the Lovell baselines
		if 'Lo' in row.anname:
			numLO = row.nosta
		if 'Mk2' in row.anname:
			numM2 = row.nosta

	u=[]
	v=[]
	proj=[]
	basel = []
	baseline = 10000000
	for visibility in data:
		if ( visibility.baseline != [numLO, numM2] ) and ( visibility.baseline != [numM2, numLO] ):
			u.append(visibility.uvw[0])
			v.append(visibility.uvw[1])
			newbasel = (299792458.0/data.header.crval[2]) * math.sqrt((visibility.uvw[0] ** 2) + (visibility.uvw[1] ** 2))
			proj.append( newbasel )
			if newbasel < baseline:
				baseline = newbasel
				basel = visibility.baseline


	antab = data.table('AN',1)
	for row in antab:
		if row.nosta == basel[0]:
			ant1=row.anname

	for row in antab:
		if row.nosta == basel[1]:
			ant2=row.anname

	print "For projected baseline", baseline/1000, "km, between", ant1, "and", ant2


	# Perley & Butler 2012 values
	A = 1.2515
	B = -0.4605
	C = -0.1715
	D = 0.0336


	log10f = (math.log(freq)/2.3025851) - 3.0; # Why the -3? Because freq has to be GHz for the formula to work.
	log_flux = A + B*log10f + C*log10f*log10f + D*log10f*log10f*log10f
	vlaflux = math.pow(10.0,log_flux)

	ref_bl_length = 11236.79 # MK-TA separation in metres.
	ref_freq = 5000.0
	ref_rho = 0.04

	bl_length = baseline

	frac = (freq / ref_freq) * (bl_length / ref_bl_length)
	rho = frac * frac * ref_rho
	merlinflux = vlaflux / (1.0 + rho)

	print "\tfor IF with freq =", freq, ", e-MERLIN flux =", merlinflux

	return merlinflux



########################################################################################################################################
# Flagging tasks
########################################################################################################################################


def flagmask(uvdata):

	"Function to apply flags from a mask of known bad frequencies, derived by the e-MERLIN team."

	afile = 'flagmask512.fg'
	if os.path.isfile(afile) :
		flagfile = open(afile, "r")
	else :
		print "Error:" + afile + "does not exist, don't know what to flag."
		print "Error: If you do not have this file, ask the e-MERLIN team."
		sys.exit()

	start = 0

	space = re.compile(r'\s')
	newline = re.compile(r'\n')

	uvflg = AIPSTask('UVFLG')

	fqtab=uvdata.table('FQ',1)
	band_strt = uvdata.header.crval[2]
	nchan = uvdata.header.naxis[2]
	if uvdata.header.naxis[3] > 1:
		ifstep = fqtab[0].total_bandwidth[0]
		band_stop = uvdata.header.crval[2] + (len(fqtab[0].total_bandwidth) * ifstep)
		chan_wdth = fqtab[0].total_bandwidth[0] / nchan
	else:
		ifstep = fqtab[0].total_bandwidth
		band_stop = uvdata.header.crval[2] + ifstep
		chan_wdth = fqtab[0].total_bandwidth / nchan


	for line in flagfile:
		line = newline.sub(r'', line)
		if start == 0:
			print line
		if line == "***BEGIN*PASS***":
			start = 1
		elif line == "***END*PASS***":
			print line
			start = 0
		elif start == 1:
			freq1 = float(space.sub(r'', line[87:98]))
			#print freq1, band_strt, band_stop, chan_wdth, ifstep
			if (freq1 >= band_strt) and (freq1 <= band_stop):
				ant1 = int(space.sub(r'', line[50:52]))
				if uvdata.header.naxis[3] > 1:
					ifs1 = int( 1 + ((freq1 - uvdata.header.crval[2]) / ifstep) )
					chans1 = int( (freq1 - (band_strt + (ifstep * (ifs1 - 1)))) / chan_wdth)
				else:
					ifs1 = 1
					chans1 = int( (freq1 - band_strt) / chan_wdth)

				start = 2
				print line[0:49], str(ant1).rjust(2), line[53:76], str(ifs1).rjust(2), line[81:86], str(chans1).rjust(12), line[99:140]
				pass
		elif start == 2:
			freq2 = float(space.sub(r'', line[87:98]))
			if (freq2 >= band_strt) and (freq2 <= band_stop):
				ant2 = int(space.sub(r'', line[50:52]))
				if uvdata.header.naxis[3] > 1:
					ifs2 = int( 1 + ((freq2 - uvdata.header.crval[2]) / ifstep) )
					chans2 = int( (freq2 - (band_strt + (ifstep * (ifs2 - 1)))) / chan_wdth)
					if chans2 == 0:
						chans2 = nchan
				else:
					ifs2 = 1
					chans2 = int( (freq2 - band_strt) / chan_wdth)
					
				# Deal with commands which cross a band edge
				if chans2 < chans1 :
					
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs1).rjust(2), line[81:86], str(nchan).rjust(12), line[99:140]
					print line[0:49], str(ant1).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(1).rjust(12), line[99:140]
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(chans2).rjust(12), line[99:140]

					uvflg.indata = uvdata
					uvflg.sources[1:] = ''
					uvflg.timerang[1:] = 0,0,0,0,999,23,59,59
					uvflg.bchan = chans1
					uvflg.echan = nchan
					uvflg.bif = ifs1
					uvflg.eif = ifs1
					uvflg.antennas[1] = ant1
					uvflg.baseline[1] = ant2
					uvflg.opcode = 'FLAG'
					uvflg.reason = 'e-MERLIN flag mask'
					uvflg.go()

					uvflg.indata = uvdata
					uvflg.sources[1:] = ''
					uvflg.timerang[1:] = 0,0,0,0,999,23,59,59
					uvflg.bchan = 1
					uvflg.echan = chans2
					uvflg.bif = ifs2
					uvflg.eif = ifs2
					uvflg.antennas[1] = ant1
					uvflg.baseline[1] = ant2
					uvflg.opcode = 'FLAG'
					uvflg.reason = 'e-MERLIN flag mask'
					uvflg.go()

				else:
					print line[0:49], str(ant2).rjust(2), line[53:76], str(ifs2).rjust(2), line[81:86], str(chans2).rjust(12), line[99:140]
					uvflg.indata = uvdata
					uvflg.sources[1:] = ''
					uvflg.timerang[1:] = 0,0,0,0,999,23,59,59
					uvflg.bchan = chans1
					uvflg.echan = chans2
					uvflg.bif = ifs1
					uvflg.eif = ifs2
					uvflg.antennas[1] = ant1
					uvflg.baseline[1] = ant2
					uvflg.opcode = 'FLAG'
					uvflg.reason = 'e-MERLIN flag mask'
					uvflg.go()

				start = 1

	return


########################################################################################################################################
# Wide-field tasks
########################################################################################################################################

def runsplit(uvdata,widetarget):
	# This task splits off each target from the final uvdataset
	# first need to SPLIT off each target and give it the name of the source
	print "Splitting sources off muilt-source file."
	split = AIPSTask('SPLIT')
	split.indata = uvdata
	split.docalib=100
	split.outclass='SPLIT' # set the outclass to 'SPLIT'
	split.outdisk=uvdata.disk
	split.outseq=0
	split.gainuse=0
	split.flagver=0
	split.doband=-1
	split.bchan=1
	split.echan=0
	split.bif=1
	split.eif=0
	split.bpver=1
	split.sources[1] = widetarget	#choose the source from widetarget
	split.go()


def chessboard64wide(uvdata,radius):
	# This function creates a chessboardlike set of smaller uvfiles.
	# Can be used for rapid imaging.	
	# radius is the half the half power beamwidth of the array.
	# uvdata is a list containing name,class,disk,sequence
	print "Chopping up uv data into chessboard-like uv datasets. This will take some time..."

	# set main input parameters
	uvfix = AIPSTask('UVFIX')
	uvfix.indata=uvdata

	#now shift the phase centre to make 4 quadrants using uvfix
	print "Shifting the phase centre to make quadrant 1."
	uvfix.shift[1]=radius*60/2
	uvfix.shift[2]=radius*60/2 
	uvfix.outdata=AIPSUVData(uvdata.name,'UV1000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 2."
	uvfix.shift[1]=-radius*60/2
	uvfix.shift[2]=radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV2000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 3."
	uvfix.shift[1]=radius*60/2
	uvfix.shift[2]=-radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV3000',uvdata.disk,1)
	uvfix.go()
	print "Shifting the phase centre to make quadrant 4."
	uvfix.shift[1]=-radius*60/2
	uvfix.shift[2]=-radius*60/2
	uvfix.outdata=AIPSUVData(uvdata.name,'UV4000',uvdata.disk,1)
	uvfix.go()

	# now average each quadrant in time
	uvavg=AIPSTask('UVAVG')
	uvavg.yinc=2
	for i in range(1,5):
		klass="UV"+ str(i) +"000"
		uvavg.indata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		uvavg.outdata=AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		print "Time averaging quadrant "+str(i)
		uvavg.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		uvdata.rename(uvdata.name,klass,uvdata.disk)

	# now average each quadrant in frequency
	avspc=AIPSTask('AVSPC')
	avspc.avoption='SUBS'
	avspc.channel=2 #averaging
	for i in range(1,5):
		klass="UV"+ str(i) +"000"
		avspc.indata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		avspc.outdata=AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		avspc.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,2)
		uvdata.rename(uvdata.name,klass,uvdata.disk)



	# **********************************************************************
	# Now shift the phase centre for each quadrant again making 16 files
	# **********************************************************************
	
	for i in range(1,5):
		klass='UV'+str(i)+'000'
		outkl='UV'+str(i)+'100'
		aipsuvdata=AIPSUVData(uvdata.name,klass,uvdata.disk,1)
		uvfix.indata=aipsuvdata
		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'100',uvdata.disk,1)
		uvfix.shift[1]=radius*60/4
		uvfix.shift[2]=radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'200',uvdata.disk,1)
		uvfix.shift[1]=-radius*60/4
		uvfix.shift[2]=radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'300',uvdata.disk,1)
		uvfix.shift[1]=radius*60/4
		uvfix.shift[2]=-radius*60/4
		uvfix.go()

		uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+'400',uvdata.disk,1)
		uvfix.shift[1]=-radius*60/4
		uvfix.shift[2]=-radius*60/4
		uvfix.go()
		uvdata = AIPSUVData(uvdata.name,klass,uvdata.disk,1)	# zap the input
		uvdata.zap()


	# **********************************************************************	
	# Now average these 16 files in time to 4 seconds
	uvavg.yinc=4 #desired average in seconds
	uvavg.zinc=2 #previous average in seconds
	for i in range(1,5):
		for j in range(1,5):
			uvavg.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvavg.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			print "Time averaging subquadrant "+str(i)
			uvavg.go()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvdata.zap()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk)

	# **********************************************************************
	# Now channel average by a factor of 2 again
	avspc.channel=2
	for i in range(1,5):
		for j in range(1,5):
			avspc.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			avspc.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			print "Averaging channels in subquadrant"
			avspc.go()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,1)
			uvdata.zap()
			uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk,2)
			uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+"00",uvdata.disk)


	# **********************************************************************
	# Now shift the phase centre for each quadrant again making 64 files
	# **********************************************************************
	
	for i in range(1,5):
		for j in range(1,5):
			aipsuvdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'00',uvdata.disk,1)
			uvfix.indata=aipsuvdata
			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'10',uvdata.disk,1)
			uvfix.shift[1]=radius*60/8
			uvfix.shift[2]=radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'20',uvdata.disk,1)
			uvfix.shift[1]=-radius*60/8
			uvfix.shift[2]=radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'30',uvdata.disk,1)
			uvfix.shift[1]=radius*60/8
			uvfix.shift[2]=-radius*60/8
			uvfix.go()

			uvfix.outdata=AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'40',uvdata.disk,1)
			uvfix.shift[1]=-radius*60/8
			uvfix.shift[2]=-radius*60/8
			uvfix.go()

			uvdata = AIPSUVData(uvdata.name,'UV'+str(i)+str(j)+'00',uvdata.disk,1)
			uvdata.zap()


	# **********************************************************************	
	# Now average these 64 files in time to 8 seconds
	uvavg.yinc=8 #desired average in seconds
	uvavg.zinc=4 #previous average in seconds
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				uvavg.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvavg.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				print "Time averaging subsubquadrant "+str(i)+str(j)+str(k)
				uvavg.go()
				uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvdata.zap()
				uvdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				uvdata.rename(uvdata.name,"UV"+ str(i) + str(j) + str(k) + "0", uvdata.disk)


	# **********************************************************************
	# Now channel average by a factor of 2 again
	avspc.channel=2
	for i in range(1,5):
		for j in range(1,5):
			for k in range(1,5):
				avspc.indata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				avspc.outdata=AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				print "Averaging channels in subsubquadrant"
				avspc.go()
				uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,1)
				uvdata.zap()
				uvdata = AIPSUVData(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk,2)
				uvdata.rename(uvdata.name,"UV"+ str(i) + str(j)+str(k)+"0",uvdata.disk)

	print "Completed averaging of " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq) + " into 64 uvfiles."





def widemaps(uvdata,indisk,radius,target,noiter):
	# This images each uv file created in chessboard
	# First set up imaging defaults	
	# Enter permanent parameters here:

	# The next 10 lines moved here from the main code.  Added findmaxb() call to calculate correct pixel size for any freq.  mkargo 20140515
	cellsiz = findmaxb(uvdata)
	cellpitch = cellsiz[0]

	imgsize = int((radius*2*60/(cellpitch*8)))
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = 2**count
	#print "imgsize =", imgsize, "i =", count, "radius =", radius, "cellpitch =", cellpitch

	docalib=-1
	gainuse=0
	flagver=0
	doband=-1
	bpver=0
	bchan=1
	echan=0
	nchav=128
	chinc=nchav
	dotv=-1

	imagr = AIPSTask('IMAGR')
	imagr.sources[1] = target
	imagr.outname = uvdata.name
	imagr.outdisk = indisk
	imagr.docalib = docalib
	imagr.gainuse = gainuse
	imagr.flagver = flagver
	imagr.doband = doband
	imagr.bpver = bpver
	imagr.bchan = bchan
	imagr.echan = echan
	imagr.nchav = nchav
	imagr.chinc = chinc
	imagr.cellsize[1] = cellpitch
	imagr.cellsize[2] = cellpitch
	imagr.imsize[1] = imgsize
	imagr.imsize[2] = imgsize
	imagr.niter = noiter
	imagr.dotv = dotv
	imagr.robust = 7

	indata=uvdata
	print "Now imaging field from UV file " + uvdata.name + '.' + uvdata.klass + '.' + format(uvdata.seq)
	imagr.indata = indata
	imagr.go()




def flat(imdata,fovradius,target):

	# Determine the cellsize used from the image header
	cellpitch = abs(imdata.header.cdelt[0] * 3600)

	imgsize = int((fovradius*2*60/(cellpitch*8)))
	thing = imgsize
	count=2
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch
	while thing > 1.0 :
		thing = thing / 2
		count = count + 1
	imgsize = 2**count
	#print "imgsize =", imgsize, "i =", count, "radius =", fovradius, "cellpitch =", cellpitch

	
 	# First trim the images
	print 'Trimming images'
 	subim=AIPSTask('SUBIM')
	subsize=int(fovradius*60*2/(cellpitch*8))
	subblc=int((imgsize-subsize)/2)
	subtrc=int(imgsize-subblc)
	#print "imsize, subsize, subblc, subtrc =", imgsize, subsize, subblc, subtrc
	buffe=24 #overlap in pixels to account for projection effects
 	subim.blc[1]=subblc-buffe
	subim.blc[2]=subblc-buffe
	subim.trc[1]=subtrc+buffe
	subim.trc[2]=subtrc+buffe
	subim.xinc=1
	subim.yinc=1
	subim.indata = imdata
	suboutname='SUB' + target[-4:-1]+target[-1]
	subim.outname=suboutname
	subim.outseq=0
	if imdata.seq < 10 :
		subim.outclass='ICL00' + format(imdata.seq)
	else:
		subim.outclass='ICL0' + format(imdata.seq)
	print 'trimming facet '+ format(imdata.seq)
	subim.go()




def pbcorr(imdata):
 	# beam correct using pbcor
	print 'Beam correction based typical WTMOD values...'
 	pbcor=AIPSTask('PBCOR')
	pbcor.indata=imdata
	pbcor.inclass='FLATN'
	pbcor.pbparm[1]=0.3 # minimum extent of beam to believe
	pbcor.pbparm[2]=1 # use the parameters below
	pbcor.pbparm[3]=-4.968
	pbcor.pbparm[4]=99.99
	pbcor.pbparm[5]=-72.78			# eventually we want to change these dynamically, right?
	pbcor.pbparm[6]=0.72
	pbcor.pbparm[7]=7.806
	pbcor.outclass='PBCOR'
	pbcor.outname=imdata.name
	pbcor.outseq=0
	pbcor.outdisk=imdata.disk
	pbcor.go()

def SAD(imdata):
	# Catalogue the bright sources in the field.
	print "Adjusting dimensions to create an RMS map..."
	regrd=AIPSTask('REGRD')
	regrd.imsize[1]=16384 #maximum size for RMSD to work.
	regrd.imsize[2]=16384
	regrd.inname=imdata.name
	regrd.outname=imdata.name
	regrd.inclass='PBCOR'
	regrd.outclass='REGRD'
	regrd.outseq=0
	regrd.outdisk=imdata.disk
	regrd.go()

	print "Creating RMS map..."
	rmsd=AIPSTask('RMSD')
	rmsd.imsize[1]=50 #size of area to get the rms
	rmsd.imsize[2]=50
	rmsd.inname=imdata.name
	rmsd.outname=imdata.name
	rmsd.inclass='REGRD'
	rmsd.outclass='RMSD'
	rmsd.outseq=0
	rmsd.outdisk=imdata.disk
	rmsd.go()

	print "Searching field for bright sources..."
	sad=AIPSTask('SAD')
	sad.inname=imdata.name
	sad.outname=imdata.name
	sad.inclass='REGRD' #the regridded image
	sad.in2name=imdata.name
	sad.in2class='RMSD' #the rms map
	sad.in2disk=imdata.disk
	sad.inseq=0
	sad.in2disk=imdata.disk
	sad.dparm[9]=10 #use the rms map in S/N ratio
	sad.ngauss=500 #maximum number of components to search for
	sad.stvers=0 #create an ST table
	sad.go()



print "Specific e-MERLIN procedures loaded! Running pipeline."
