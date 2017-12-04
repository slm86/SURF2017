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

def std(sdir='.',smin=0,smax=0,bin_nr=0):
    # Plot the standard deviation of dust-to-gas density fluctuations 
    # as a function of time
    s = np.zeros(smax+1)
    for snum in np.arange(smin,smax+1):
            print(snum)
            ss=pipeline.snap_ext(snum,four_char=1)
            #filename = './scratch/dust_snap_'+ss+'.h5' 
            filename = sdir+'/dust_snap_'+ss+'.h5'
            
            infi=h5py.File(filename,'r')
            #nd   = np.array(infi["Number_Density_List_Of_Dust_Neighbors"])
            ng   = np.array(infi["Number_Density_List_Of_Gas_Neighbors"])
            # stuff specific to size bins 
            bins = np.array(infi["Size_Bins"])
            dset_name = 'Particle_Positions_Bin'+str(bin_nr)
            ok = np.array(infi[dset_name])
            nds  = np.array(infi["Number_Density_List_Of_Dust_Neighbors_Size"])
            infi.close()
            
            x = np.log10(nds[bin_nr][ok]/ng[ok])
            
#            w = nds[bin_nr][ok]**(-1)
#            n, bins= np.histogram(x,bins=100,weights=w,normed=True)
            s[snum]=np.sqrt(np.var(x))
            
    num = np.arange(smin,smax+1)
    ax.plot(num,s-s[0],label=str(bin_nr))

sdir = '../c.k05N4/'   
smax = 50
 
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for i in [0,1,2,3]:
    std(sdir=sdir,smin=0,smax=smax,bin_nr=i)

ax.set_xlim(left=0)
ax.set_xlim(right=smax)
ax.set_ylim(bottom=0.01)
ax.set_ylim(top=0.5)

# why aren't these working on log scale?
#ax.set_yticks([0.1,0.5])
#ax.set_yticklabels(['me','you'])
#ax.set_xticks([1,10,20])
#ax.set_xticklabels(['1','10','100'])

ax.set_yscale('log')

ax.set_xlabel('snapshot number')
ax.set_ylabel('sigma(log10(n_dust/n_gas))')
plt.legend()
data = 'data '+sdir
ax.text(smax,0.1, data, rotation=-90)

plt.title('Time evolution of the standard deviation of density fluctuations')
plt.show()