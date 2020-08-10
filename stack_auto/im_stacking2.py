""" example usage:
python ~/ss/src/MysteryDecrement/stack-gen_function.py <radiomap> <stackpos.npy> <postagestampsize>

where, 

stackpos.npy is 3col: ID, RA, DEC

postagestampsize: in pixels (note: half the size - 0.5, so 50, generates stamp 101 x 101 pixels

ipython -i ~/ss/src/MysteryDecrement/stack-gen_function.py ~/ss/maps/VLA_VIDEO_I_080916_resid.fits somepos.npy 50
"""
import pyfits as pf
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import numpy as np
import sys
import os
from scipy.signal import convolve
from scipy.optimize import curve_fit
import time
import pickle
from astLib import astWCS
from astropy.io import fits
import json
#neddir = '/home/deane/scripts/general'
#sys.path.append(neddir)
#from nedwright import *
import ConfigParser
import subprocess

if __name__=='__main__':
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[-1])

    PATH = config.get('pipeline', 'data_path') 


    print('Image stacking Groundtruth image ...........')

    ###################
    ### General options  ".txt"
    ###################
    inputmap = PATH+config.get('field', 'fitsname')  # sys.argv[1] # radio map
    inputcat = PATH+config.get('pipeline', 'project_name')+'_truthcat_fixed.txt' # sys.argv[2]  # fits binary catalogue
    subsize = 50  # int(sys.argv[3]) # size of stacked subimages MINUS 1 and DIVIDED by 2.
    radiocube_subsize = 50 # int(sys.argv[4]) # radiocube to write out size (MINUS 1 DIVIDED by 2)

    if (subsize < radiocube_subsize):
        print('\n\n\n!!! ERROR: subsize must be larger than radiocube_subsize !!!!\n\n')
        sys.exit()

    plotdir =PATH+'image_plotdir_' + inputcat.split('/')[-1].replace('.npy','')
    os.system('mkdir %s'%plotdir)


    if radiocube_subsize:
        print('\n\n!!!! NOTE: saving radiocubes, which takes much longer. Estimated radiocubesize = %i GB '\
              %((radiocube_subsize*2 + 1)/61. * 33))
        time.sleep(5)
        makeradiocube = True
        radiocubename = 'stackall_' + inputmap[:-5] + '_%ipix.npy'%(2*radiocube_subsize + 1)
        makenoisecube = True
        noisecubename = 'noise_stackall_' + inputmap[:-5] + '_origcube%ipixels.npy'%(2*subsize + 1)
    else:
        makeradiocube = False
        makenoisecube = False




    randomise = False # randomise positions?
    if (randomise):
        scatter_arcsec = 20.  # ideally, one would use maps bounds to avoid edge and local rms effects

    FITthreshold = 3 # only attempt 2D Gauss fit if |SNR| > FITthreshold

    ############
    ### LOAD DATA
    #############
    radiocat = np.loadtxt(inputcat, usecols =(1,2))
    radioRA = radiocat[:, 0]
    radioDEC = radiocat[:, 1]

    vlamap = np.squeeze(pf.getdata(inputmap)) # units: Jy/beam
    hdr = pf.getheader(inputmap)
    VLAWCS=astWCS.WCS(hdr,mode='pyfits')
    ### WCS parameters
    refRA = hdr['CRVAL1']
    refDEC = hdr['CRVAL2']
    refpixRA = hdr['CRPIX1'] - 1 # to make first pixel = 0
    refpixDEC = hdr['CRPIX2'] -1 # to make first pixe = 0
    numpixRA = hdr['NAXIS1']
    numpixDEC = hdr['NAXIS2']
    try:
        pixscale = hdr['CDELT1']*(-1)  # assumes square pixels
    except KeyError:
        pixscale = hdr['CD1_1']*(-1)  # assumes square pixels

    if (len(sys.argv) > 7): # could use parser module
        print('\n\n assuming that %s is the rms map\n\n'%sys.argv[7])
        noisemap_path = sys.argv[7]
        noisemap =  np.squeeze(pf.getdata(noisemap_path))
        usenoisemap = True
    else:
        usenoisemap = False
        print('\n\n no rms map provided, will estimate from subim, size=%i, npixels = %i'%(subsize,subsize**2))








    def randomise_pos(radioRA,radioDEC,scatter_arcsec):
        radioRA += np.random.rand(len(radioRA))*(scatter_arcsec/3600.)  \
          - (scatter_arcsec/3600./2.)
        radioDEC += np.random.rand(len(radioRA))*(scatter_arcsec/3600.) \
          - (scatter_arcsec/3600./2.)

    def remove_coords_outside_radiomap(radioRA,radioDEC): # this doesn't seem to work for COSMOS map
        ra_temp, dec_temp = [],[]
        for obj in range(len(radioRA)):
            if (VLAWCS.coordsAreInImage(radioRA[obj],radioDEC[obj])):
                ra_temp.append(radioRA[obj])
                dec_temp.append(radioDEC[obj])

        radioRA = np.array(ra_temp)
        radioDEC = np.array(dec_temp)

        return radioRA, radioDEC




    ############################
    ### Plotting/fittig fuctions
    ############################

    def oneGaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, bpa, offset):
        sigma_y = sigma_x    
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(bpa)**2)/(2*sigma_x**2) + (np.sin(bpa)**2)/(2*sigma_y**2)
        b = -(np.sin(2*bpa))/(4*sigma_x**2) + (np.sin(2*bpa))/(4*sigma_y**2)
        c = (np.sin(bpa)**2)/(2*sigma_x**2) + (np.cos(bpa)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( -1.0* (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
        return g.ravel()

    def fit1g(subim):
        subsize = int((len(subim) - 1)/2.)
        peakf = subim.min()
        y=range(len(subim[0]))
        x=range(len(subim[0]))
        x, y = np.mgrid[-subsize:subsize+1, -subsize:subsize+1]
        xguess=np.where(subim==subim.min())[0][0]
        yguess=np.where(subim==subim.min())[1][0]
        sigma_guess = 5.
        initial1=[peakf,xguess-len(subim)/2 ,yguess-len(subim)/2,sigma_guess,sigma_guess,0,0] 
        wes=subim.reshape(subim.size,)

        try:
            popt1, pcov1 = curve_fit(oneGaussian,(x,y),wes, p0=initial1)
            data_fitted = oneGaussian((x, y), *popt1).reshape(len(x[0]),len(x[0]))
            resid =subim - data_fitted

        except RuntimeError:
            data_fitted = subim*0
            resid = subim - data_fitted
            popt1,pcov1 = 0,0
        return data_fitted,popt1,pcov1,resid


    def plotstack_data(infits,subsize,ind):

        hdr = pf.getheader(infits)
        subim = np.squeeze(pf.getdata(infits))     * 1e9
        try:
            pixscale = hdr['CDELT2']*3600 # arcsec
        except KeyError:
            pixscale = hdr['CD2_2']*3600 # arcsec
        npix = hdr['NAXIS2']
        ex = npix * pixscale / 2. # arcsec
        extent = [-ex,ex,-ex,ex]
        cmap = pl.cm.RdGy
        vmin = subim.min()
        vmax = subim.max()


        pl.figure()
        pl.imshow(subim,origin='lower',interpolation='nearest',extent=extent,cmap=cmap)#,vmin=vmin,vmax=vmax)
        pl.xlabel('$\Delta$RA / arcsec')
        pl.ylabel('$\Delta$Dec. / arcsec')
        pl.colorbar(label='nJy')
        pl.savefig(infits[:-5]+'-data.png')





    def fit_and_plotstack_modelresid(infits,subsize,ind):
        """
        Logical flow:
        if pixel |SNR| < 2.5, dont' fit
        else: 
             attempt fit
             if good, good. 
             elif fakefit:
                 check ints, set bestparms = 0
             elif failedfit:
                  set besparms = 0


           case a: successful fit, actual component
           case b: successful fit, spurious
           case c: unsuccessful fit, return integer
           case d: don't even try based on abs(peak)/rms

        abs_peak = max(abs(submin.min(),submin.max()))
        if (abs_peak < 2.5*median_std):
            nofit, record zeros
        else:
            plotstack_modelresid()

        """



        hdr = pf.getheader(infits)
        subim = np.squeeze(pf.getdata(infits))     * 1e9

        try:
            data_fitted, popt1,pcov1,resid = fit1g(subim)
            goodfit = True
        except IndexError: 
            goodfit = False

        if goodfit:

            pixscale = hdr['CDELT2']*3600 # arcsec
            npix = hdr['NAXIS2']
            ex = npix * pixscale / 2. # arcsec
            extent = [-ex,ex,-ex,ex]

            cmap = pl.cm.RdGy
            vmin = subim.min()
            vmax = subim.max()

            pl.figure()
            pl.imshow(subim,origin='lower',interpolation='nearest',extent=extent,cmap=cmap)#,vmin=vmin,vmax=vmax)
            pl.xlabel('$\Delta$RA / arcsec')
            pl.ylabel('$\Delta$Dec. / arcsec')
            pl.colorbar(label='nJy')
            pl.savefig(infits[:-5]+'-data.png')

            if 'std' not in infits:

                pl.figure()
                pl.imshow(data_fitted,origin='lower',interpolation='nearest',extent=extent,cmap=cmap) #,vmin=vmin,vmax=vmax)
                pl.xlabel('$\Delta$RA / arcsec')
                pl.ylabel('$\Delta$Dec. / arcsec')
                if (type(popt1)==int):
                    pl.title('Gaussian fits FAILED')
                else:
                    pl.title('Gauss FWHM = %.1f arcsec; DC offset = %.2f nJy'%(popt1[3]*2.35*pixscale,popt1[6]))
                pl.colorbar(label='nJy')
                pl.savefig(infits[:-5]+'-model.png')

                pl.figure()
                pl.imshow(resid,origin='lower',interpolation='nearest',extent=extent,cmap=cmap) #,vmin=vmin,vmax=vmax)
                pl.xlabel('$\Delta$RA / arcsec')
                pl.ylabel('$\Delta$Dec. / arcsec')
                pl.colorbar(label='nJy')
                pl.savefig(infits[:-5]+'-resid.png')


                f = open(infits[:-5]+'_fitsummary.txt','w')
                print>>f, '#fit results:'
                print>>f, '#amp_nJy,x0,y0,fwhm,DCoffset_nJy'
                if(type(popt1) != int):
                    print>>f, '%.2f \t %.2f \t %.2f \t %.2f \t %.2f'%\
                      (popt1[0],popt1[1]*pixscale,popt1[2]*pixscale,\
                        popt1[3]*2.35*pixscale,popt1[6])
                    best_fit_parms = popt1[0],popt1[1]*pixscale,popt1[2]*pixscale,\
                                     popt1[3]*2.35*pixscale,popt1[6]
                elif (type(popt1) == int):
                    best_fit_parms = 0,0,0,0,0
                    print>>f, '%.2f \t %.2f \t %.2f \t %.2f \t %.2f'%best_fit_parms

                f.close()

                return best_fit_parms

                #print>>f, '\namp = %.2f  nJy'%(popt1[0])
                #print>>f, 'x0 = %.2f arcsec; y0 = %.2f arcsec'%(popt1[1]*pixscale,popt1[2]*pixscale)
                #print>>f, 'FWHM = %.2f arcsec'%(popt1[3]*2.35*pixscale)
                #print>>f, 'DC offset = %.2f nJy'%(popt1[6])





    ########################
    ### STACKING function
    ########################

    def stack(radioRA,radioDEC):
        print 'total number of stacks = %i'%(len(radioRA))
        print 'initialising...'
        cube = np.zeros([num_stacks,2*subsize+1,2*subsize+1])
        if (usenoisemap):
            noisecube = np.zeros([num_stacks,2*subsize+1,2*subsize+1])
            noise = 'dummy'
        else:
            noise = np.zeros(num_stacks)
            noisecube = 'dummy'
        ind = np.ones(num_stacks,dtype='int') # =1 if included in median/avg calcs

        start_time = time.time()

        for gal in range(num_stacks):
            RApixoffset,DECpixoffset = astWCS.WCS.wcs2pix(VLAWCS,radioRA[gal],radioDEC[gal])
            center = np.array([int(np.round(RApixoffset)), int(np.round(DECpixoffset))]) 
            subim = vlamap[center[1]-subsize : center[1]+subsize+1,center[0]-subsize : center[0]+subsize+1]
            if (cube[0,:,:].shape == subim.shape):
                cube[gal,:,:] = subim
            else:
                cube[gal,:,:] = np.zeros([2*subsize+1,2*subsize+1])
                subim = np.zeros([2*subsize+1,2*subsize+1])
                ind[gal] = 0
            #if (subim.max() > maxsrcflux) or (subim.min() < -maxsrcflux):
            #    ind[gal] = 0

            if (usenoisemap):
                noise_subim = noisemap[center[1] - subsize:center[1] + subsize + 1,\
                                           center[0] - subsize:center[0] + subsize + 1]
                if (cube[0,:,:].shape == noise_subim.shape):
                    noisecube[gal,:,:] = noise_subim
                else:
                    noise_subim = np.zeros([2*subsize+1,2*subsize+1])
                    noisecube[gal,:,:] = noise_subim
                #if (np.nanmean(noise_subim) > 2*clipflux):
                #    ind[gal] = 0
            else:
                noise_subim_std = np.std(subim)
                noise[gal] = noise_subim_std
                #if (abs(noise[gal]) > clipflux):
                #    ind[gal] = 0

            if (np.isnan(subim.sum() ) or (np.count_nonzero(subim) == 0)): # removes any Nan values
                ind[gal] = 0

            #cube[gal,:,:] = subim   #### WTF???

            if not (gal%10000):
                print inputmap.split('/')[-1] + ':\t\t' + str(gal + 1) + ' of ' + str(num_stacks) + ' objects'
                print("avg for last 10k stacks = %s seconds \n" % (np.round((time.time() - start_time),decimals=2)))

                start_time = time.time()

        ind = np.where(ind == 1)[0]
        return cube, noisecube, noise, ind


    ###################################
    ####### 
    ###################################
    def make_stacked_maps():
        median_map = np.zeros([2*subsize + 1, 2*subsize + 1])
        for i in range(2*subsize + 1):
            for j in range(2*subsize + 1):
                median_map[i,j] = np.nanmedian(cube[ind,i,j]) 

        median_std_map = np.zeros([2*subsize + 1, 2*subsize + 1])
        for i in range(2*subsize + 1):
            for j in range(2*subsize + 1):
                median_std_map[i,j] = np.nanstd(cube[ind,i,j])/np.sqrt(len(ind)) / 1.48


        avgmap = np.zeros([2*subsize + 1, 2*subsize + 1])
        for i in range(2*subsize + 1):
            for j in range(2*subsize + 1):
                if (usenoisemap):
                    avgmap[i,j] = np.average(cube[ind,i,j],weights=1./noisecube[ind,i,j]**2)
                else:
                    avgmap[i,j] = np.average(cube[ind,i,j],weights=1./noise[ind]**2)

        cube_std_map = np.zeros([2*subsize + 1, 2*subsize + 1])
        for i in range(2*subsize + 1):
            for j in range(2*subsize + 1):
                cube_std_map[i,j] = np.nanstd(cube[ind,i,j])/np.sqrt(len(ind))


        return median_map,median_std_map,avgmap, cube_std_map




    ###################################
    ####### save file with sample sizes
    ###################################
    def save_stack_results(niter):


        #np.savetxt('indices_stacked_loopniter_%03d.txt'%niter,ind,fmt='%i')
        win = subsize - 10
        #median, median_std = median_map[subsize,subsize], median_std_map[subsize,subsize]  #units are Jy
        if (median_map[subsize,subsize] < 0):
            #median, median_std = median_map[win:-win,win:-win].min(), median_std_map[subsize,subsize]  #units are Jy
            median, median_std = median_map[subsize,subsize], median_std_map[subsize,subsize]  #units are Jy
        else:
            #median, median_std = median_map[win:-win,win:-win].max(), median_std_map[subsize,subsize]
            median, median_std = median_map[subsize,subsize], median_std_map[subsize,subsize]
        if (usenoisemap):
            wmean, wmean_std = np.average(cube[ind,subsize,subsize],weights=1./noisecube[ind,subsize,subsize]**2), cube[ind,subsize,subsize].std() /np.sqrt(len(ind)) 
        else:
            wmean, wmean_std = np.average(cube[ind,subsize,subsize],weights=1./noise[ind]**2), cube[ind,subsize,subsize].std() /np.sqrt(len(ind)) 

        median_str = '%.2e \pm %.2e (S/N = %.1f), '%(median,median_std,median/median_std)
        wmean_str = '%.2e \pm %.2e Jy/b (S/N = %.1f)'%(wmean,wmean_std,wmean/wmean_std)
        app = ' total = ' + str(len(ind))  + '\nweighted mean = ' + wmean_str + ', median = ' + median_str +' \n'

        f = open(PATH+'stats.txt','a')
        f.write('input catalogue:  ' + inputcat.split('/')[-1])
        f.write('input map:        ' + inputmap.split('/')[-1])
        f.write(app)
        f.write('____________________________________\n')
        f.close()

        hdr['NAXIS'] = 4
        hdr['NAXIS1'] = subsize*2 + 1
        hdr['NAXIS2'] = subsize*2 + 1
        hdr['NAXIS3'] = 1
        hdr['NAXIS4'] = 1
        hdr['CRVAL1'] = 0.0
        hdr['CRVAL2'] = 0.0
        hdr['CRPIX1'] = subsize+1
        hdr['CRPIX2'] = subsize+1
        hdr['BUNIT'] = 'Jy/beam'
        hdr_cube = hdr.copy()
        hdr_cube['NAXIS3'] = num_stacks #! should consider only saving ind, rather than num stacks?
        ### now update header for output fits files ###

        fitsname = os.path.join(plotdir,inputmap[:-5].split('/')[-1]+'_median_map'\
           +'-subsize%i-num_stacks%i.fits'%(subsize,len(ind)))
        pf.writeto(fitsname,median_map.reshape(1,1,2*subsize+1,2*subsize+1),header=hdr,clobber=True)
        plotstack_data(fitsname,subsize,ind)
        # determine whether you want to fit 
        abs_peak = max(abs(median_map.min()),median_map.max())
        if (abs_peak < FITthreshold * np.mean(median_std_map)):
            best_fit_parm_median = np.zeros(5) #! hard-code warning
        else:
            best_fit_parm_median = fit_and_plotstack_modelresid(fitsname,subsize,ind)

        fitsname = os.path.join(plotdir,inputmap[:-5].split('/')[-1]+'_median_std_map'\
          +'-subsize%i-num_stacks%i.fits'%(subsize,len(ind)))
        pf.writeto(fitsname,median_std_map.reshape(1,1,2*subsize+1,2*subsize+1),header=hdr,clobber=True)
        plotstack_data(fitsname,subsize,ind)
        #fit_and_plotstack_modelresid(fitsname,subsize,ind)

        fitsname = os.path.join(plotdir,inputmap[:-5].split('/')[-1]+'_avgmap'\
          +'-subsize%i-num_stacks%i.fits'%(subsize,len(ind)))
    #    pf.writeto(fitsname,(cube[ind,:,:].sum(axis=0)/len(ind)).reshape(1,1,2*subsize+1,2*subsize+1),header=hdr,clobber=True)
        pf.writeto(fitsname,avgmap.reshape(1,1,2*subsize+1,2*subsize+1),header=hdr,clobber=True)
        plotstack_data(fitsname,subsize,ind)
        #fit_and_plotstack_modelresid(fitsname,subsize,ind)

        fitsname = os.path.join(plotdir,inputmap[:-5].split('/')[-1]+'_avgmap_std'\
           +'-subsize%i-num_stacks%i.fits'%(subsize,len(ind)))
        pf.writeto(fitsname,cube_std_map.reshape(1,1,2*subsize+1,2*subsize+1),header=hdr,clobber=True)
        plotstack_data(fitsname,subsize,ind)
        #fit_and_plotstack_modelresid(fitsname,subsize,ind)

        #if (writeMasterCubes):
        #    pf.writeto('mastercube_subsize%i.fits'%\
        #                   subsize,cube.reshape(1,num_stacks,2*subsize+1,2*subsize+1),header=hdr_cube,clobber=True)
        #    pf.writeto('masterNoisecube_subsize%i.fits'%\
        #                   subsize,noisecube.reshape(1,num_stacks,2*subsize+1,2*subsize+1),header=hdr_cube,clobber=True)

        fitsnames = os.path.join(plotdir,inputmap[:-5].split('/')[-1]+'_median_map'\
           +'-subsize%i-num_stacks%i.fits'%(subsize,len(ind)))
        os.system('mv'+' '+fitsnames+' '+PATH+'Imstacked-model.fits')

        return best_fit_parm_median, median_std





    ###################################
    ##### main functions called below
    ###################################


    ### initialise cat
    #print('\n pre-loading input catalog...\n')
    #radioRA, radioDEC, preloaded_cat,preloaded_hdulist = initialise_input_cat(inputcat,cattype)





    niter = 0



    if (randomise): # add artifical positional scatter?
        radioRA, radioDEC = randomise_pos(radioRA,radioDEC,scatter_arcsec)
    # remove coords outside radio map; but doesn't work for COSMOS map...
    #radioRA, radioDEC = remove_coords_outside_radiomap(radioRA,radioDEC)

    ### and stack...
    num_stacks = len(radioRA)

    cube, noisecube, noise, ind = stack(radioRA,radioDEC)
    median_map,median_std_map,avgmap,cube_std_map = make_stacked_maps()
    best_fit_parm_median, median_std = save_stack_results(niter)
    #amp_nJy, x0, y0, fwhm, DCoffset_nJy, median_centpix, median_std, num_stacks, J-K colour

    print('\nmaking noise cube...')


    if makeradiocube:
        np.save(PATH+'%s_%ix%ipixels.npy'%(inputmap.split('/')[-1],radiocube_subsize*2+1,radiocube_subsize*2+1), \
                cube[:,subsize-radiocube_subsize:subsize+radiocube_subsize+1,\
                       subsize-radiocube_subsize:subsize+radiocube_subsize+1])
    if makenoisecube:
        noise_out = np.zeros(len(cube[:, 0, 0]))

        for subim in range(len(cube[:, 0, 0])):
            noise_out[subim] = np.nanstd(cube[subim, :, :])
        stampsize = 2 * subsize + 1
        np.save(PATH+'noise_%s_stampsize%i.npy' % (inputmap.split('/')[-1],(stampsize * 2 + 1)), noise_out)


    print('completed stacking')



       
