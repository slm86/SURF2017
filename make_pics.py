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

# Produce pictures of all snapshots in a simulation and create a GIF 

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;

sdir = './scratch'
snum_min = 1
snum_max = snum_min
zmin = 0.510
zmax = 0.525
dt = 0.025
# don't forget to adjust colorbar min/max scale below
vmin = -0.5
vmax = 0.5

x_subtract = 0
y_subtract = 0
z_subtract = 0

xg, yg = np.mgrid[0:1:512j, 0:1:512j]
for snum in np.arange(snum_min,snum_max+1):
    print(snum)
    H  = readsnap(sdir,snum,0,header_only=1);    
    time = H['time']
    print('time',time)
    P  = readsnap(sdir,snum,0);
    Pd = readsnap(sdir,snum,3);
    x = P['p'][:,0]
    y = P['p'][:,1]
    z = P['p'][:,2]
    vx = P['v'][:,0]
    vy = P['v'][:,1]
    vz = P['v'][:,2]
    rho = P['rho']
    
    xd = Pd['p'][:,0]
    yd = Pd['p'][:,1]
    zd = Pd['p'][:,2]
    vxd = Pd['v'][:,0]
    vyd = Pd['v'][:,1]
    rd = Pd['R_d']
    md = Pd['m']
    
    ok, = np.where((z>zmin)&(z<zmax))    
    okd, = np.where((zd>zmin)&(zd<zmax))
    
    x = x[ok] - x_subtract
    y = y[ok] - y_subtract
    z = z[ok] - z_subtract
    rho = rho[ok]
    
    xd = xd[okd] - x_subtract
    yd = yd[okd] - y_subtract
    zd = zd[okd] - z_subtract
    rd = rd[okd]
    
    y_subtract += np.average(vy)*dt
    x_subtract += np.average(vx)*dt
    z_subtract += np.average(vz)*dt
    
    # this is esentially the code of dustywind_plot.plot_some_things()
    f = plt.figure(figsize=(12,6))
    ax1 = f.add_axes([0.1,0.1,0.4,0.8])
    ax2 = f.add_axes([0.5,0.1,0.4,0.8])#, sharex=ax1)
    
    # change if box is rectangular
    # xg, yg = np.mgrid[0:1:512j, 0:2:1024j]
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    ax2.set_yticklabels('')
    
    # variable of interest, to be plotted
    w = np.log10(rho)
    print(np.amin(w),np.amax(w),'sizes',np.amin(rd),np.amax(rd))
    cmap = 'seismic'#'jet'
   
    dg=interpolate.griddata((x,y),w,(xg,yg),method='linear',fill_value=np.median(w))
    im=ax1.imshow(np.flipud(np.transpose(dg)),interpolation='bicubic',\
                    vmin=vmin,vmax=vmax,cmap=cmap,extent=(0,1,0,2),aspect='auto');
                  
    ax1.set_xticks(np.linspace(0,1,5,endpoint=False))
    ax_ins = inset_axes(ax1,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.92,.90,0.05),\
                        bbox_transform=ax1.transAxes,borderpad=0.15)
    
    cb=pylab.colorbar(im,cax=ax_ins,ticks=[-2,-1,-0.5,0,0.5,1,2],orientation='horizontal')
    cb.ax.tick_params(labelsize=10)
  
    cmap_dust ='viridis'
    im2 = ax2.scatter(xd,yd,c=(rd),s=(100*rd**1.5),cmap=cmap_dust,rasterized=True); #5*rd #300*rd**2
   
    f.subplots_adjust(wspace=0,hspace=0)
    plt.show()
    
    ss=snap_ext(snum,four_char=0)
    filename = './fig'+str(ss)+'.png'
    plt.savefig(filename, format='png', dpi=100)
    plt.close('all')

