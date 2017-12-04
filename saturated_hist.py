import sys
sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86\analysis')
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from readsnap import readsnap
import matplotlib.colors as colors
import pipeline as pipeline
import matplotlib.ticker as ticker

def saturated_data_2d(sdir='.',smin=0,smax=0,bin_nr=0):
    # Produce the 2d histogram of number densities in a given size bin by
    # combining data from all snapshots after reaching a saturated turbulent state
    nd_list = []
    ng_list = []
    nds_list = []
    for snum in np.arange(smin,smax+1):
            #pipeline.scratch()
            #pipeline.grab_wheeler(sdir,snum,snap_type=1)
            ss=pipeline.snap_ext(snum,four_char=1)
            #filename = './scratch/dust_snap_'+ss+'.h5'
            
            filename = sdir+'/dust_snap_'+ss+'.h5'
            
            infi=h5py.File(filename,'r')
            #nd   = np.array(infi["Number_Density_List_Of_Dust_Neighbors"])
            ng   = np.array(infi["Number_Density_List_Of_Gas_Neighbors"])
            #hd   = np.array(infi["Smoothing_Length_List_Of_Dust_Neighbors"])
            #hg   = np.array(infi["Smoothing_Length_List_Of_Gas_Neighbors"])
            #nd_list.append(nd)
            #ng_list.append(ng)
            
            # stuff specific to size bins 
            bins = np.array(infi["Size_Bins"])
            dset_name = 'Particle_Positions_Bin'+str(bin_nr)
            ok = np.array(infi[dset_name])
            nds  = np.array(infi["Number_Density_List_Of_Dust_Neighbors_Size"])
            nds_list.append(nds[bin_nr][ok])
            ng_list.append(ng[ok])
            infi.close()
            
    #nd_list = np.ravel(np.array(nd_list))
    nds_list = np.ravel(np.array(nds_list))
    ng_list = np.ravel(np.array(ng_list))
    
    x = np.log10(ng_list)
    y = np.log10(nds_list)    
    w = nds_list**(-1)
    plt.hist2d(x,y,bins=100,norm=colors.LogNorm(),weights=w)
    plt.xlabel('log10(n_gas)')
    plt.ylabel('log10(n_dust)')
    s0 = sdir.split("/")
    title = 'Volume weighted; data:'+ s0[len(s0)-1] + '\n' + 'bin nr.'+str(bin_nr)
    plt.title(title)
    plt.ylim([-2.5,4.5])
    plt.xlim([-3,3.5])
    
    ax = plt.gca()
    #ax.axvline(x=0, ls='dotted', lw=1.3, color='k')
    ax.axhline(y=0, ls='dotted', lw=1.3, color='k')
    x = np.arange(-5,5)
    plt.plot(x,x,ls='dotted',lw=1.3,color='k')
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%03.1f'))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%03.1f'))
    plt.show()

def saturated_data_1d(sdir='.',smin=0,smax=0,bin_nr=0):
    # Produce the 1d histogram of number densities in a given size bin by
    # combining data from all snapshots after reaching a saturated turbulent state
    nd_list = []
    ng_list = []
    nds_list = []
    for snum in np.arange(smin,smax+1):
            ss=pipeline.snap_ext(snum,four_char=1)
            #filename = './scratch/dust_snap_'+ss+'.h5' 
            filename = sdir+'/dust_snap_'+ss+'.h5'
            
            infi=h5py.File(filename,'r')
            #nd   = np.array(infi["Number_Density_List_Of_Dust_Neighbors"])
            ng   = np.array(infi["Number_Density_List_Of_Gas_Neighbors"])
            #hd   = np.array(infi["Smoothing_Length_List_Of_Dust_Neighbors"])
            #hg   = np.array(infi["Smoothing_Length_List_Of_Gas_Neighbors"])
            #nd_list.append(nd)
            #ng_list.append(ng)
            
            # stuff specific to size bins 
            bins = np.array(infi["Size_Bins"])
            dset_name = 'Particle_Positions_Bin'+str(bin_nr)
            ok = np.array(infi[dset_name])
            nds  = np.array(infi["Number_Density_List_Of_Dust_Neighbors_Size"])
            nds_list.append(nds[bin_nr][ok])
            ng_list.append(ng[ok])
            infi.close()
            
    #nd_list = np.ravel(np.array(nd_list))
    nds_list = np.ravel(np.array(nds_list))
    ng_list = np.ravel(np.array(ng_list))
    
    x = np.log10(nds_list/ng_list)
    w = nds_list**(-1)
    plt.hist(x,bins=100,weights=w,histtype='step',normed=True)
    ax = plt.gca()
    ax.set_yscale('log',basey=10)    
    
    plt.xlabel('log10(n_dust/n_gas)')
    s0 = sdir.split("/")
    title = 'Volume weighted; data:'+ s0[len(s0)-1] + '\n' + 'bin nr.'+str(bin_nr)
    plt.title(title)
    #plt.ylim([-2.5,3.5])
    #plt.xlim([-2,2])
    plt.show()
    
##########   MAIN : CALL YOUR FUNCTIONS HERE  ##########
plt.figure()
saturated_data_2d(sdir='../c.k05mu.01',smin=5,smax=29,bin_nr=0)


