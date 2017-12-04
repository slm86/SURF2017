import sys
sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86\analysis\scratch')
sys.path.append('C:\Users\Stefi\Desktop\pfh_python-routines\gadget_lib')
# change to these when working on wheeler, or appropriate for your machine
#sys.path.append('/panfs/ds08/hopkins/slm86/pfh_python-routines/gadget_lib')
#sys.path.append('/panfs/ds08/sxs/jsquire/gizmo-dust')
import numpy as np
from readsnap import readsnap
import time
import h5py as h5py

# these are some random 'utility' routines needed by the master routine below
def snap_ext(snum,four_char=0):
    ext='00'+str(snum);
    if (snum>=10): ext='0'+str(snum)
    if (snum>=100): ext=str(snum)
    if (four_char==1): ext='0'+ext
    if (snum>=1000): ext=str(snum)
    return ext;

# here is the actual master routine
def compute_rdf_wr(x, y, z, x2, y2, z2, vx, vy, vz, vx2, vy2, vz2, S, radii, dr, rMax, Nt=10):
    """Compute the three-dimensional radial distribution function and the radial 
    velocities between a pair of dust particle species contained in a box with
    side lengths specified by vector S. 
    Arguments:
        x, x2               an array of x positions of centers of particles
        y, y2               an array of y positions of centers of particles
        z, z2               an array of z positions of centers of particles
        vx, vx2             corresponding x-component velocities
        vy, vy2             corresponding y-component velocities
        vz, vz2             corresponding z-component velocities
        S                   3-component vector with dimensions of the box in space
        radii               vector specifying the radii at which we want to work out
                            the RDF and radial velocities. This is done by considering 
                            particles inside spherical shells at these radii.
        dr                  vector specifying the width of spherical shells
        
        Nt              parameter that affects the speed of the algorithm 
                        -- don't change default unless for good reasons
                        ## Nt is the number of small tiles that could fit in the x direction in the box
                        ## it should be at least 2, which is the case with one small tile and one large tile
                        ## if you try 1, the large tile will not fit the box, and because of the way this is 
                        ## written you will get a divide by zero error
                        
    Returns a tuple: (g_average, wr_average)
        g(r)            a numpy array containing the correlation function (RDF) g(r)
        wr(r)           a numpy array containing the radial velocity function wr(r)
    """
    from numpy import zeros, sqrt, where, pi, mean, arange
    twoD = False
    box_volume = S[0]*S[1]*S[2]
        
    if((np.amin(z)==0) & (np.amax(z)==0)):
        twoD = True
        box_volume = S[0]*S[1]
        S[2] = 0
    
    numberDensity = len(x) / box_volume
    num_radii = radii.size
    shell_volume = zeros(num_radii)
    for i in range(num_radii):
        rOuter = radii[i] + dr[i]/2
        rInner = radii[i] - dr[i]/2
        print(rInner,radii[i],rOuter)
        if twoD: 
            shell_volume[i]= pi * (rOuter**2 - rInner**2)
        else: 
            shell_volume[i]= (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))
            
    g_average = zeros(num_radii)
    wr_average = zeros(num_radii)
    
    # Tile the box
    # We do this in order to speed up the algorithm; rather than centering on
    # every particle of species 1 in the box and searching the particles of
    # species 2 at the given distances, we divide the box into tiles and loop
    # over these. Since the radii we are interested in are smaller than 1/10 of 
    # the box length, it's only necessary to look at the particles in a tile.
    Nt_x = Nt
    st = np.float(S[0])/np.float(Nt_x) # size of a tile
    Nt_y = Nt*S[1]/S[0]
    Nt_z = Nt*S[2]/S[0]
    
    tiles_x = np.linspace(st/2,S[0]-st/2,Nt_x)
    tiles_y = np.linspace(st/2,S[1]-st/2,Nt_y)
    if(twoD==True): 
        Nt_z = 2
        tiles_z = zeros(Nt_z+1)
    else:           
        tiles_z = np.linspace(st/2,S[2]-st/2,Nt_z)
    Nt_total = (Nt_x-1)*(Nt_y-1)*(Nt_z-1) # total number of tiles
    
    # Main loop over tiles to compute RDF
    for tx in np.arange(Nt_x-1):
        for ty in np.arange(Nt_y-1):
            for tz in np.arange(Nt_z-1):
                # select species 1 particles inside current tile
                index_inner, = np.where((tiles_x[tx]<=x)&(x<=tiles_x[tx+1])&\
                                        (tiles_y[ty]<=y)&(y<=tiles_y[ty+1])&\
                                        (tiles_z[tz]<=z)&(z<=tiles_z[tz+1]))
                # select species 2 particles inside extended tile
                index_outer, = np.where((tiles_x[tx]-st/2<=x2)&(x2<=tiles_x[tx+1]+st/2)&\
                                        (tiles_y[ty]-st/2<=y2)&(y2<=tiles_y[ty+1]+st/2)&\
                                        (tiles_z[tz]-st/2<=z2)&(z2<=tiles_z[tz+1]+st/2))

                num_interior_particles = len(index_inner)   
                num_outer_particles = len(index_outer)
                g = zeros([num_interior_particles, num_radii])
                wr = zeros([num_interior_particles, num_radii])
                r = zeros((num_outer_particles,3))
                w = zeros((num_outer_particles,3))
                
                # Compute pairwise correlation for each interior particle
                for p in range(num_interior_particles):
                    index = index_inner[p]
                    # compute distances and relative velocities between the current
                    # particle p and all particles labeled by index_outer
                    r[:,0] = x[index] - x2[index_outer]
                    r[:,1] = y[index] - y2[index_outer]
                    r[:,2] = z[index] - z2[index_outer]
                    w[:,0] = vx[index] - vx2[index_outer]
                    w[:,1] = vy[index] - vy2[index_outer]
                    w[:,2] = vz[index] - vz2[index_outer]
                    d = sqrt(r[:,0]**2 + r[:,1]**2 + r[:,2]**2)
                    for i in arange(num_radii):
                        # count the particles that are inside the shell ar radii[i]
                        dummy, = np.where((radii[i]-dr[i]/2 < d) & (d < radii[i]+dr[i]/2))
                        result = dummy.size
                        g[p,i] = np.float(result) / numberDensity
                        # calculate the radial component of the relative velocity
                        # w between particle p and all those in the shell
                        # wr = (r.w)/|r| where r.w is the dot product
                        wr_dummy = np.sum(r[dummy,:]*w[dummy,:])/d[dummy]
                        # select only those particles approaching
                        neg, = np.where(wr_dummy<0)
                        if(len(neg)>0): 
                            wr[p,i] = mean(wr_dummy[neg])                       
                        else: wr[p,i] = 0
    
                # Average g(r) and wr(r) over all particles p inside this tile
                # we are summing the results over all tiles, and will divide by
                # the tile number at the end
                for i in range(num_radii):
                    g_average[i] += mean(g[:, i]) / shell_volume[i]
                    wr_average[i] += mean(wr[:,i]) 
                    
    g_average = g_average/Nt_total
    wr_average = wr_average/Nt_total
    return (g_average, wr_average)
#### END OF FUNCTION ####

def loop_rdf_wr(sdir,snum,num_bins,box_size=[1,1,1]):
    """This computes the radial distribution function(g(r)) and radial velocities
    between grains of size bin i and those of size bin j, at separations given by
    vector 'radii'. The results are stored as 3 dimensional matrices; 
    the first two axes are just spanning the possible i-j pairs, so it will be 
    lower triangular; the third direction runs over the separations at which g 
    and wr are computed"""
    Pd = readsnap(sdir,snum,3);
    xd = Pd['p'][:,0]
    yd = Pd['p'][:,1]
    zd = Pd['p'][:,2]
    vxd = Pd['v'][:,0]
    vyd = Pd['v'][:,1]
    vzd = Pd['v'][:,2]
    
    # Shells
    # to be modified such that these numbers are not hard coded, but passed on 
    # as a variable or something..
    radii = np.array([0.0025,0.005,0.01,0.025,0.05])
    dr = np.array([0.0025,0.0025,0.0075,0.0225,0.0275])
    num_radii = radii.size
    
    # Construct the dust size bins and the vectors with their indices, ok_ds
    # The size bins are equally spaced in log-space
    # ok_ds is a list of vectors, each vector containing the indices of dust 
    # particles in a certain size bin
    e_min=np.amin(Pd['R_d'])
    e_max=np.amax(Pd['R_d'])
    size_bins = np.linspace(np.log10(e_min*0.99),np.log10(e_max*1.01),num=num_bins+1,endpoint=True)
    ok_ds=[]
    for i in np.arange(num_bins):
        dummy, = np.where((np.log10(Pd['R_d'])>=size_bins[i])&(np.log10(Pd['R_d'])<size_bins[i+1]))
        ok_ds.append(dummy)
    
    g_list = np.zeros((num_bins,num_bins,num_radii))
    wr_list = np.zeros((num_bins,num_bins,num_radii))
    for i in np.arange(num_bins):
       for j in np.arange(i+1):
           print('Computing bins',i,j)
           # This is the call to the main routine to compute g(r) and wr(r)
           # between dust species i and dust species j
           g,wr = compute_rdf_wr(xd[ok_ds[i]],yd[ok_ds[i]],zd[ok_ds[i]],\
                                 xd[ok_ds[j]],yd[ok_ds[j]],zd[ok_ds[j]],\
                                 vxd[ok_ds[i]],vyd[ok_ds[i]],vzd[ok_ds[i]],\
                                 vxd[ok_ds[j]],vyd[ok_ds[j]],vzd[ok_ds[j]],\
                                 box_size, radii, dr, 0.1, Nt=10)
           g_list[i,j,:] = g
           wr_list[i,j,:] = wr
    return (g_list,wr_list,radii)

# MAIN program
# Specify snapshot directory, snapshot number, and number of dust species
sdir = './scratch/'
num_bins = 10

snum_min = 0
snum_max = 100
snum_step = 150

for snum in np.arange(snum_min,snum_max + snum_step,snum_step):
    start_time = time.time()
    # Call the main routine
    g_result, wr_result, r = loop_rdf_wr(sdir,snum,num_bins)

    # Write results to h5 file
    ss=snap_ext(snum,four_char=1)
    filename=sdir+'/rdf_snap_'+ss+'.h5'
    outfi = h5py.File(filename,'w')
    dset_radii = outfi.create_dataset('Radii',data=r)
    dset___rdf = outfi.create_dataset('RDF',data=g_result)
    dset____wr = outfi.create_dataset('wr',data=wr_result)
    outfi.close()
    print('Created h5 file for snap ',snum)
    print('Duration=',time.time()-start_time(),' s')
