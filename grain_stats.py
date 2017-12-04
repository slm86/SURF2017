import sys
sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86')
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from readsnap import readsnap
# In this script we use new snapshot files created with grain_density_from_snapshot.py
# to plot various dust properties

# read dust h5 file
sdir = "./scratch"
snum = 100
infi=h5py.File('C:\Users\Stefi\AppData\Local\lxss\home\slm86\c.k05N4\dust_snap_0100.h5','r')

########### LOAD DATA ###########

#sdir = infi["Snapshot_Directory"].value
#snum = infi["Snapshot_Number"].value # note that for scalars above, we need to add the '.value' to pull out the actual number
time = infi["Snapshot_Time"].value 
Ndus = infi["Total_Number_Of_Dust_Particles"].value 
Ngas = infi["Total_Number_Of_Gas_Particles"].value 

nd   = np.array(infi["Number_Density_List_Of_Dust_Neighbors"])
ng   = np.array(infi["Number_Density_List_Of_Gas_Neighbors"])
hd   = np.array(infi["Smoothing_Length_List_Of_Dust_Neighbors"])
hg   = np.array(infi["Smoothing_Length_List_Of_Gas_Neighbors"])

# stuff specific to size bins
bins = np.array(infi["Size_Bins"])
ok_ds = []
for i in np.arange(bins.size - 1):
    dset_name = 'Particle_Positions_Bin'+str(i)
    ok_ds.append(np.array(infi[dset_name]))
nds  = np.array(infi["Number_Density_List_Of_Dust_Neighbors_Size"])
hds  = np.array(infi["Smoothing_Length_List_Of_Dust_Neighbors_Size"])

infi.close()

#readsnap('file_directory',snapshot_number,particle_type,cosmological=?,skip_bh=?,header_only=?)
H  = readsnap(sdir,snum,0,header_only=1)
P  = readsnap(sdir,snum,0)
Pd  = readsnap(sdir,snum,3)

x = P['p'][:,0]
y = P['p'][:,1]
m = P['m']
rho = P['rho']    
vx = P['v'][:,0]
vy = P['v'][:,1]
vz = P['v'][:,2]
v  = np.sqrt(vx**2+vy**2+vz**2)

idd= Pd['id']
xd = Pd['p'][:,0]
yd = Pd['p'][:,1]
vxd = Pd['v'][:,0]
vyd = Pd['v'][:,1]
vzd = Pd['v'][:,2]
vd  = np.sqrt(vxd**2+vyd**2+vzd**2)
md = Pd['m']
rd = Pd['R_d']

# Create grid for plotting
xg,yg = np.mgrid[0:1:2048j, 0:2:4096j]
#xg,yg = np.mgrid[0:1:512j, 0:2:1024j]

####################################
def density_histogram():
    # Dust-volume weighted histogram of dust-to-gas density fluctuations
    # All dust is considered here, no size differentiation
    ax = plt.gca()
    x = np.log10(ng/nd)
    w = nd**(-1)
    plt.hist(x,bins=50,weights=w,histtype='step',label=time)
    ax.set_yscale('log',basey=10)
    plt.xlabel('log10(n_gas/n_dust)')
    plt.xlim([-4,3])
    title ='Volume weighted' +' data:'+sdir+'\n'+'time '+str(np.round(time,decimals=3))
    plt.title(title)
    plt.show()

def density_histogram_1d(bin_nr=0):
    # Dust-volume weighted histogram of dust-to-gas density fluctuations
    # One species of dust only - those grains in size bin: bin_nr
    plt.figure()
    ax = plt.gca()
    #ok = np.transpose(np.array(ok_ds[bin_nr]))[:,0]
    ok = np.ravel(ok_ds[bin_nr])
    
    # Histogram of the ratio n_dust(bin_nr)/n_gas, at the positions of grains in bin_nr
    # Histogram is weighted by the dust volume ~ 1/number density
    x = np.log10(nds[bin_nr][ok]/ng[ok])
    w = nds[bin_nr][ok]**(-1)
    plt.hist(x,bins=100,weights=w,histtype='step')
    
    ax.set_yscale('log',basey=10)
    plt.xlabel('log10(n_dust/n_gas)')
    title = 'Bin number '+str(bin_nr)+'/'+str(bins.size-2)+'; weighted by volume'
    plt.title(title)
    s = 'data '+sdir
    plt.text(1.0,1.0, s, rotation=-90)
    plt.show()
    
    
def density_histogram_1d_all():
    # Dust-volume weighted histograms of dust-to-gas density fluctuations for
    # each dust species, displayed on the same plot
    plt.figure()
    ax=plt.gca()
    for bin_nr in np.arange(bins.size - 1):
        ok = np.transpose(np.array(ok_ds[bin_nr]))[:,0]
        #ok = np.ravel(ok_ds[bin_nr])
        x = np.log10(nds[bin_nr][ok]/ng[ok])
        w = nds[bin_nr][ok]**(-1)
        a,bb = np.histogram(x,bins=100,weights=w)
        b = np.zeros(a.size)
        for j in np.arange(a.size):
            b[j]=(bb[j]+bb[j+1])/2
        plt.plot(b,a,'+',label=str(bin_nr))
    plt.legend(title='Size Bins')
    #ax.set_yscale('log',basey=10)
    plt.xlabel('log10(n_dust/n_gas)')
    plt.ylabel('Probability density')
    title ='Volume weighted' +' data:'+sdir+'\n'+'time '+str(np.round(time,decimals=3))
    plt.title(title)
    ax.set_xlim([-2,2])
    
def density_histogram_2d(bin_nr=0):
    # Dust-volume weighted histogram of dust and gas neighbor number densities
    # One species of dust only - those grains in size bin: bin_nr
    plt.figure()
    ax = plt.gca()
    ok = np.transpose(np.array(ok_ds[bin_nr]))[:,0]
    x = np.log10(ng[ok])
    y = np.log10(nds[bin_nr][ok])    
    w = nds[bin_nr][ok]**(-1)
    plt.hist2d(x,y,bins=100,norm=colors.LogNorm(),weights=w)
    plt.xlabel('log10(n_gas)')
    plt.ylabel('log10(n_dust)')
    title = 'Bin number '+str(bin_nr)+'/'+str(bins.size-2)+'; weighted by volume'
    plt.title(title)
    s = 'data '+sdir
    plt.text(1.0,1.0, s, rotation=-90)
    ax.axvline(x=0, ls='dotted', lw=2, color='k')
    plt.show()
    
def velocity_size_pdf():
    # Dust streaming velocity dependence on grain size
    # Streaming velocity is calculated as box average dust velocity minus
    # box average gas velocity in the y-direction
    plt.figure()
    ax = plt.gca()
       
    x = np.log10(vd/rd)
    plt.xlabel('log10(v_dust/R_d)')
    
    xx = (vyd-np.average(vy))
    #plt.xlabel('v_stream/R_d')
    
    plt.hist(x,bins=20,normed=True)  
    plt.ylabel('Probability density')
    plt.show()
    
    plt.figure()
    plt.scatter(rd,np.abs(xx),s=0.1,rasterized=True)
    ax2 = plt.gca()
    ax2.set_xscale('log',basex=10)
    ax2.set_yscale('log',basey=10)
    #ax2.set_xlim([0,0.1])
    ax2.set_ylim([0.1,20])
    ax2.set_ylabel('Dust Streaming Velocity')
    ax2.set_xlabel('Grain Size')

def some_function():
    # Check mass normalization and magnitude of some velocities
    md0 = np.sum(nds[0]* hds[0]**2)
    print('mass in bin1',md0)

    md1 = np.sum(nds[1]* hds[1]**2)
    print('mass in bin2',md1)

    md2 = np.sum(nds[2]* hds[2]**2)
    print('mass in bin3',md2)

    md3 = np.sum(nds[3]* hds[3]**2)
    print('mass in bin4',md3)

    print('avg rho dust',np.average(nds))
    mg = np.sum(ng* hg**2)
    print('gas total mass ',mg)
    print('dust/gas mass ratio',(np.sum(nds* hds**2))/mg)

    print('dust density fluctuations',((nds[0]+nds[1]+nds[2]+nds[3])-nd)/nd)
    
    print('mean gas vx',np.average(vx))
    print('mean gas vy',np.average(vy))
    print('mean dust vx',np.average(vxd))
    print('mean dust vy',np.average(vyd))

def plot_routine():
    # Produce colormap plots of number density
    fig = plt.figure(figsize=(12,6))
    cmap = 'jet'
    vmin = -1.5
    vmax = 1.5
    # colormap of gas density
#    ax1 = plt.subplot(121)
#    w = np.log10(rho)
#    dg=interpolate.griddata((x,y),w,(xg,yg),method='linear',fill_value=np.median(w))
#    im=pylab.imshow(np.flipud(np.transpose(dg)),interpolation='bicubic',\
#                    vmin=vmin,vmax=vmax,cmap=cmap,extent=(0,1,0,2));
#    title1 = "log10(rho_gas)"+' at time '+str(np.round(H['time'],decimals=2))
#    plt.title(title1)
#    plt.xlabel('x axis')
#    plt.ylabel('y axis')
#    ax_ins = inset_axes(ax1,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),
#                        bbox_transform=ax1.transAxes,borderpad=0.15)
#    cb=pylab.colorbar(im,cax=ax_ins,orientation='horizontal')
#    cb.ax.tick_params(labelsize=7)
    
    # colormap of dust density
#    ax2 = plt.subplot(122)
#    dummy = np.log10(nd)
#    dg2 = interpolate.griddata((xd,yd),dummy,(xg,yg),method='linear',fill_value=np.median(dummy))
#    im2 = pylab.imshow(np.flipud(np.transpose(dg2)),interpolation='bicubic', \
#                       vmin=vmin,vmax=vmax,cmap=cmap,extent=(0,1,0,2)); 
#    ax_ins=inset_axes(ax2,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),\
#                      bbox_transform=ax2.transAxes,borderpad=0.15)
#    cb=pylab.colorbar(im2,cax=ax_ins,orientation='horizontal')
#    cb.ax.tick_params(labelsize=7)
#    ax2.set_title('log10(nd)')
    
    # scatter of all dust particles, color-coded by size
    ''' ax1 = plt.subplot(1,3,1)
    plt.scatter(xd,yd,s=0.2,c=np.log(rd),cmap='jet',rasterized=True);
    plt.xlim([0,1])
    plt.ylim([0,2])
    ax1.set_title('All Dust Particles')'''
    
    #?see how to make colorbar in log scale
    # colormap of dust density, for grains in the smallest size bin
    
    ax20 = plt.subplot(1,4,1)
    dummy = np.log10(nds[0])
    dg2 = interpolate.griddata((xd,yd),dummy,(xg,yg),method='linear',fill_value=np.median(dummy))
    im2 = pylab.imshow(np.flipud(np.transpose(dg2)),interpolation='bicubic', \
                       vmin=-2.5,vmax=2.5,cmap=cmap,extent=(0,1,0,2)); 
    ax_ins=inset_axes(ax20,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),\
                      bbox_transform=ax20.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im2,cax=ax_ins,orientation='horizontal')
    cb.ax.tick_params(labelsize=7)
    ax20.set_title('Bin 0')
    ax20.set_yticks([0,0.5,1,1.5,2])
    ax20.set_xticks([0,0.5])
    
    # colormap of dust density, for grains in the bin 1
    ax21 = plt.subplot(1,4,2)
    dummy = np.log10(nds[1])
    dg3=interpolate.griddata((xd,yd),dummy,(xg,yg),method='linear',fill_value=np.median(dummy))
    im3=pylab.imshow(np.flipud(np.transpose(dg3)),interpolation='bicubic',\
                    vmin=-2.5,vmax=2.5,cmap=cmap,extent=(0,1,0,2)); 
    ax_ins = inset_axes(ax21,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),\
                        bbox_transform=ax21.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im3,cax=ax_ins,orientation='horizontal')
    cb.ax.tick_params(labelsize=7)
    ax21.set_title('Bin 1')
    ax21.set_yticks([])
    ax21.set_xticks([0,0.5])
    
    # colormap of dust density, for grains in the largest size bin
    ax22 = plt.subplot(1,4,3)
    dummy = np.log10(nds[2])
    dg3=interpolate.griddata((xd,yd),dummy,(xg,yg),method='linear',fill_value=np.median(dummy))
    im3=pylab.imshow(np.flipud(np.transpose(dg3)),interpolation='bicubic',\
                    vmin=-2.5,vmax=2.5,cmap=cmap,extent=(0,1,0,2)); 
    ax_ins = inset_axes(ax22,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),\
                        bbox_transform=ax22.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im3,cax=ax_ins,orientation='horizontal')
    cb.ax.tick_params(labelsize=7)
    ax22.set_title('Bin 2')
    ax22.set_yticks([])
    ax22.set_xticks([0,0.5])
    
    # colormap of dust density, for grains in the largest size bin
    ax23 = plt.subplot(1,4,4)
    dummy = np.log10(nds[3])
    dg3=interpolate.griddata((xd,yd),dummy,(xg,yg),method='linear',fill_value=np.median(dummy))
    im3=pylab.imshow(np.flipud(np.transpose(dg3)),interpolation='bicubic',\
                    vmin=-2.5,vmax=2.5,cmap=cmap,extent=(0,1,0,2)); 
    ax_ins = inset_axes(ax23,width="100%",height="100%",loc=1,bbox_to_anchor=(0.05,0.91,0.90,0.05),\
                        bbox_transform=ax23.transAxes,borderpad=0.15)
    cb=pylab.colorbar(im3,cax=ax_ins,orientation='horizontal')
    cb.ax.tick_params(labelsize=7)
    ax23.set_title('Bin 3')
    ax23.set_yticks([])
    ax23.set_xticks([0,0.5,1])
    
    ax21.set_xlabel('x axis')
    ax20.set_ylabel('y axis')
    s = 'data '+sdir
    plt.text(1.1,0.5, s, rotation=-90)
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)
    fig.suptitle('Colormap of log10(n_dust)')

def dust_density_custom(a_bin):
    # colormap of dust desnity of grains in one bin at positions of grains in a different bin size
    # a_bin gives the species who's density we plot
    # b_bin gives the positions at which we evaluate this density
    
    plt.figure()
    for i in np.arange(4):
        b_bin = i
        ok = np.ravel(ok_ds[b_bin])
        n  = nds[a_bin][ok]
        print(np.average(n))
        returns = plt.hist(np.log10(n), bins=300, range=(-3,3), histtype='step', label=str(b_bin), normed=False)
        
    print(np.sum(returns[0]), np.size(ok))
    
    s = 'data '+sdir
    plt.text(3.1,2500,s,rotation=-90)
    title = "Number density of grains in bin "+str(a_bin)+" around grains in other bin sizes" + \
                "\n"+"snapshot "+str(snum)
    plt.title(title)
    plt.xlabel('log10(n)')
    #plt.yscale('log')
    plt.xlim([-3,3])
    plt.ylim([0,8000])
    plt.legend()
    plt.show()
    return plt.gca()

##########   MAIN : CALL YOUR FUNCTIONS HERE  ##########

vmin = -2
vmax = 2

#some_function()
#plot_routine()

#density_histogram_1d_all()
#for i in np.arange(bins.size - 1):
    #density_histogram_1d(bin_nr=i)
    #density_histogram_2d(bin_nr=i)

velocity_size_pdf()
#density_histogram_1d_all()
#density_histogram_2d(bin_nr=3)

#dust_density_custom(0)

#g_r, r, i = pairCorrelationFunction_2D(xd,yd,[1,2],0.4,0.1)
#plt.plot(r,g_r)
#plt.show()


'''
plt.contourf(xg,yg,dg)
imd = pylab.imshow(np.flipud(np.transpose(dg)),extent=(0,1,0,2))
'''