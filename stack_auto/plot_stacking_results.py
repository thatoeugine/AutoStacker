import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pyrap.tables as pt
import astropy.io.fits as pf
from astropy.table import Table
import sys
plt.style.use("seaborn")
import configparser
from matplotlib import rc,rcParams
rc('text', usetex=True)
# activate latex text rendering
rc('axes', linewidth=2)
rc('font', weight='bold')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']



########################################################################################
#                Setting up constants and reading data from the stacked MS file
#                and the T-RECS catalogue.
########################################################################################
if __name__=='__main__':
    config = configparser.ConfigParser()
    config.read(sys.argv[-1])


    path = config.get('pipeline', 'data_path') 
    msfilename = path+config.get('casa_imager', 'stack_msfile_name')
    tab = pt.table(msfilename.encode('ascii', 'ignore')) # open stacked MS file
    vis_data = tab.getcol('CORRECTED_DATA').T #Transposing the visibility data
    vis_data2 = tab.getcol('DATA').T 
    uvw = tab.getcol('UVW').T  # convert uvdist to units of meters,
    uvdist = np.sqrt(uvw[0,:]**2 + uvw[1,:]**2)
    flag_col = tab.getcol('FLAG') # flagged data = True (such as low elevation visibilities)

    ############### Ground values #######################
    cat_data = tab.getcol('MODEL_DATA').T


    #################################################################################################
    #                      Plotting uvdist vs. real(vis)
    #################################################################################################

    uvbins_edges = np.arange(0,9,1)*1000 #(0,0.05,0.005)*1000#uvdistance units: Kilo-lambda
    uvbins_centre = (uvbins_edges[:-1] + uvbins_edges[1:])/2.
    numuvbins = len(uvbins_centre)
    binwidths = uvbins_edges[1] - uvbins_edges[0]

    # initialise
    ampbins_mean = np.zeros([numuvbins])
    stdbins = np.zeros([numuvbins])
    Nvisperbin = np.zeros([numuvbins])
    fluxbins_mean = np.zeros([numuvbins])


    # colours to use (order for all plots)
    colorlist = ['#336699', 'cyan', 'pink', 'r' ,'green', 'k', 'yellow', 'grey','orange', 'purple', 'magenta']
    #colorlist.reverse()
    if (len(colorlist) < numuvbins):
        print ('add more colors to colorlist')
        sys.exit()
    else:
        colorlist = colorlist[:numuvbins]



    for b in range(numuvbins):
        mask = (uvdist > uvbins_edges[b])&(uvdist< uvbins_edges[b+1])&(np.logical_not(flag_col[:,0,0]))# mask of unflagged visibilities in this uvbin

        Nvisperbin[b] = mask.sum() # total number of visibilities in this uvbin

        if (Nvisperbin[b] == 0):
            ampbins_mean[b],stdbins[b] = 0,0# if no vis, set to zero
        else:
            ampbins_mean[b] = np.nanmean(vis_data.real[[0,3],:,:][:,:,mask]) #average real amplitude values in bin "b"
            stdbins[b] = np.nanstd(vis_data2.real[[0,3],:,:][:,:,mask]) / Nvisperbin[b]**0.5 # rms of that bin
            fluxbins_mean[b] = np.nanmean(cat_data.real[[0,3],:,:][:,:,mask])#average flux values in bin "b" (i.e. Theoretical flux)

    # Note [0,3] are correlations XX & YY



    #######
    ### amp vs uvdist, with uncertainties For mean
    #######
    fig = plt.figure(figsize=(10,6.8))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    yerr = stdbins/np.sqrt(Nvisperbin) #noise_per_vis/np.sqrt(np.sum(Nvisperbin,axis=0)) #yerr = noise_per_vis/np.sqrt(np.sum(allsrcs[:,2,:],axis=0))
    xerr = binwidths/2. * np.ones(numuvbins)
    for b,color in enumerate(colorlist):
        ax1.semilogy(uvbins_centre[b],ampbins_mean[b],'o',color=color,mec='none',alpha=1)
        ax1.errorbar(uvbins_centre[b],ampbins_mean[b],xerr=xerr[b],yerr=yerr[b],\
                 ecolor=color,lw=1,alpha=1,fmt='none',capsize=2)
        ax2.semilogy(uvbins_centre,fluxbins_mean,'--',color='r',mec='none',alpha=0.3)
    ax1.set_xlabel('uv distance [m]', fontsize=18)
    ax1.set_ylabel('Real(amplitude) [Jy]', fontsize=18)
    #ax1.set_ylim(0,np.nanmax(ampbins_mean)*1.2)
    ax1.set_xlim(0,uvbins_edges.max())
    ax2.set_xlabel('uv distance [m]',fontsize=18)
    ax2.set_xlim(ax1.get_xlim())
    np.savetxt(path+'uvdistplot_ampdataptsMean.txt',np.vstack([uvbins_centre,xerr,ampbins_mean,yerr]))
    plt.savefig(path+'ampMean_uvdist.pdf')
    
