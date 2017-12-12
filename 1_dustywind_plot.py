import sys
#sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86\analysis')
#sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86\analysis\pfh_python-routines')
# change to these when working on wheeler, or appropriate for your machine
sys.path.append('/panfs/ds08/hopkins/slm86/pfh_python-routines/gadget_lib')
#sys.path.append('/panfs/ds08/sxs/jsquire/gizmo-dust-spectrum')
import numpy as np
import h5py as h5py
import scipy.interpolate as interpolate
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
#import pipeline
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from readsnap import readsnap
#from matplotlib import rc
#plt.rc('text', usetex=True)
#plt.rc('font',family='sans serif')

########### LOAD DATA ###########

sdir = '../output'
snum = 100

H  = readsnap(sdir,snum,0,header_only=1);    
print('time',H['time'])
P  = readsnap(sdir,snum,0);
Pd = readsnap(sdir,snum,3);

# read gas data
x = P['p'][:,0]
y = P['p'][:,1]
z = P['p'][:,2]
rho = P['rho']   
#h  = P['h'] # smoothing length  
vx = P['v'][:,0]
vy = P['v'][:,1]
vz = P['v'][:,2]
m = P['m']
v  = np.sqrt(vx**2+vy**2+vz**2)

# read dust data
xd = Pd['p'][:,0]
yd = Pd['p'][:,1]
zd = Pd['p'][:,2]
vxd = Pd['v'][:,0]
vyd = Pd['v'][:,1]
vzd = Pd['v'][:,2]
rd = Pd['GrainSize']
md = Pd['m']
vd  = np.sqrt(vxd**2+vyd**2+vzd**2)

# compute box-average velocities 
avg_v = np.average(v)
avg_vd = np.average(vd)
print('mach',avg_vd-avg_v)

print('mass',np.sum(md)/np.sum(m))

print('max log10(rho)',np.log10(np.amax(rho)))
print('min log10(rho)',np.log10(np.amin(rho)))

print('max(vd)-max(vg)',np.amax(vd)-np.amax(v))

print(np.amax(zd), np.amax(yd), np.amax(xd))

####################################

def plot_some_things(threeD=True,zmin=0.55,zmax=0.60):
    # This plots different variables of interest, obtained directly from reading 
    # the snapshots plus scatter of dust particles
    
    # default is for 3D data to take a slice zmin < z < zmax
    # for 2D data, call function with threeD=False
    
    #f,(ax1,ax2)=plt.subplots(1,2,sharey='row',sharex=True,squeeze=True)
    f = plt.figure(figsize=(12,6))
    ax1 = f.add_axes([0.1,0.1,0.4,0.8])
    ax2 = f.add_axes([0.5,0.1,0.4,0.8])#, sharex=ax1)
    xg, yg = np.mgrid[0:1:512j, 0:1:512j] #0:2:1024j]
    
    if(threeD==True):
        ok = np.where((z>zmin)&(z<zmax))    
        okd = np.where((zd>zmin)&(zd<zmax))
    else:
        ok, = np.where(x)
        okd, = np.where(xd)
    print(np.shape(ok))
    
    # ax1 holds the gas plot
    # w holds the variable of interest, to be plotted; here, gas density
    w = np.log10(rho[ok])
    vmin=np.min(w); vmax=np.max(w);
    cmap='jet'
    
    dg=interpolate.griddata((x[ok],y[ok]),w,(xg,yg),method='linear',fill_value=np.median(w))
    # adjust colorscale vmin and vmax according to the vmin,vmax found above
    im=ax1.imshow(np.flipud(np.transpose(dg)),interpolation='bicubic',\
                  vmin=-1,vmax=1,cmap=cmap,extent=(0,1,0,1),aspect='auto');
    #plt.contourf(xg,yg,dg)
    #im = pylab.imshow(np.flipud(np.transpose(dg)),extent=(0,1,0,2))
    title1 = "log10(rho_gas)"+' at time '+str(np.round(H['time'],decimals=2))
    ax1.set_title(title1) #,fontsize=16
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    ax1.set_xticks(np.linspace(0,1,5,endpoint=False))
    
    ax_ins = inset_axes(ax1,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.92,.90,0.05),\
                        bbox_transform=ax1.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im,cax=ax_ins,ticks=[-1.5,-1.0,-0.5,0,0.5,1.0,1.5],orientation='horizontal')
    cb.ax.tick_params(labelsize=10)
    
    # ax2 holds the dust plot
    # dust particles are color coded by grain size
    ax2.scatter(xd[okd],yd[okd],s=0.2,rasterized=True,c=np.log(rd[okd]),cmap='jet')
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    #plt.xlabel('x axis')
    ax2.set_title("Dust Particles")#,fontsize=16)
    '''ax_ins = inset_axes(ax2,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),
                        bbox_transform=ax2.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im2,cax=ax_ins,orientation='horizontal')'''
    cb.ax.tick_params(labelsize=9)
    ax2.set_yticklabels('')
    
    s = 'data '+sdir
    ax2.text(1.0,1.0, s, rotation=-90)
    f.subplots_adjust(wspace=0,hspace=0)
    
    fig_name = 'snap_plot_'+str(snum)+'.png'
    f.savefig(fig_name,bbox_inches='tight',pad_inches=0.1)
    print("Saved ",fig_name)


def dust_bin_over_gas(sdir,snum,bin_nr):
    # Plots the colormap of gas density and superposes a scatter of dust particles
    # from a specified size bin
    plt.figure()
    w = np.log10(rho)
    xg, yg = np.mgrid[0:1:512j, 0:2:1024j]
    cmap='jet'
    dg=interpolate.griddata((x,y),w,(xg,yg),method='linear',fill_value=np.median(w))
    im=pylab.imshow(np.flipud(np.transpose(dg)),interpolation='bicubic',vmin=-2.5,vmax=2.5,cmap=cmap,extent=(0,1,0,2));
    title1 = "log10(rho_gas)"+' at time '+str(np.round(H['time'],decimals=2))
    plt.title(title1,fontsize=12)
    plt.xlabel('x axis')
    plt.ylabel('y axis')
    
    filename = sdir + '/dust_snap_'+pipeline.snap_ext(snum,four_char=1)+'.h5'
    infi=h5py.File(filename,'r')
    bins = np.array(infi["Size_Bins"])
    ok_ds = []
    for i in np.arange(bins.size - 1):
        dset_name = 'Particle_Positions_Bin'+str(i)
        ok_ds.append(np.array(infi[dset_name]))
    nds  = np.array(infi["Number_Density_List_Of_Dust_Neighbors_Size"])
    ok = np.ravel(ok_ds[bin_nr])
    plt.scatter(xd[ok],yd[ok],s=0.2,c='k',rasterized=True)
    
    ax = plt.gca()
    ax_ins = inset_axes(ax,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),
                        bbox_transform=ax.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im,cax=ax_ins,orientation='horizontal')
    cb.ax.tick_params(labelsize=7)
    
    s = 'data '+sdir
    plt.text(1.0,1.0, s, rotation=-90)
    
##########   MAIN : CALL YOUR FUNCTIONS HERE  ##########
plot_some_things() 

#dust_bin_over_gas(sdir,snum,3)


