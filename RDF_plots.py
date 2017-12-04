''' RDF PLOT '''
import sys
sys.path.append('C:\Users\Stefi\AppData\Local\lxss\home\slm86\analysis')
sys.path.append('C:\Users\Stefi\Desktop\pfh_python-routines')
import numpy as np
#from readsnap import readsnap
import matplotlib.pyplot as plt
import h5py as h5py

sdir = './scratch'
num_bins = 10
# Read rdf_snap h5 file
name = 'c.k05N4/rdf_wr_snap_0200'
path = sdir+'/'+name+'.h5'
infi = h5py.File(path,'r') 
title = name
g_result = np.array(infi["RDF"])
wr_result = np.abs(np.array(infi["wr"]))
r = np.array(infi["Radii"])
infi.close()

# Plot results
cs = ['b','r','g','c','m']

def plot_1():
    plt.figure()
    ax = plt.gca()
    x_data = np.arange(num_bins)
    for j in np.arange(r.size):
        y_data = np.zeros(num_bins)
        for i in np.arange(num_bins):
            y_data[i] = g_result[i,i,j]
        plt.plot(x_data,y_data,'-',label=str(r[j]))

    plt.xlabel('grain size')
    plt.ylabel('autocorrelation g(r)')
#    plt.yscale('log')
    plt.ylim([0,5])
    plt.xlim([0,9.1])
    plt.xticks(np.arange(num_bins))
    plt.legend(title="r, radius=",loc=2)
    plt.title(title)
    ax.axhline(y=1,ls='dashed',lw=1,c='k') # horizontal dashed line to mark g=1
    
def plot_2():
    plt.figure()
    ax = plt.gca()
    for s in [0,1,2,3,4]:
        j = 1
        x_data1 = np.arange(s,num_bins,1)
        y_data1 = np.zeros(num_bins-s)    
        for ii in np.arange(len(x_data1)):
            i = x_data1[ii]
            y_data1[ii] = g_result[i,i-s,j]
        plt.plot(x_data1,y_data1,'-',lw=0.75,c=cs[s],label=str(s))
    #textstr = "separation radius="+str(r[j])
    #ax.text(0.60,0.95,textstr,transform=ax.transAxes, verticalalignment='top')

    for s in [0,1,2,3,4]:
        j = 2
        x_data1 = np.arange(s,num_bins,1)
        y_data1 = np.zeros(num_bins-s)    
        for ii in np.arange(len(x_data1)):
            i = x_data1[ii]
            y_data1[ii] = g_result[i,i-s,j]
        plt.plot(x_data1,y_data1,'--',lw=0.75,c=cs[s])

    for s in [0,1,2,3,4]:
        j = 4
        x_data1 = np.arange(s,num_bins,1)
        y_data1 = np.zeros(num_bins-s)    
        for ii in np.arange(len(x_data1)):
            i = x_data1[ii]
            y_data1[ii] = g_result[i,i-s,j]
        plt.plot(x_data1,y_data1,':',lw=1.2,c=cs[s])
    
    textstr2 = "separation radii"+'\n'+'solid='+str(r[1])+'\n'+'dashed='+str(r[2])+\
               '\n'+"dotted="+str(r[4])
    ax.text(0.75,0.95,textstr2,transform=ax.transAxes, verticalalignment='top')
    
    plt.xlabel('Smax, largest grain size bin')
    plt.ylabel('g(r)')
#    plt.yscale('log')
    plt.ylim([0,3])
    plt.xlim([0,9.1])
    plt.xticks(np.arange(num_bins))
    plt.legend(title="Smax - Smin=",loc=2)
    plt.title(title)
    ax.axhline(y=1,ls='dashed',lw=1,c='k') # horizontal dashed line to mark g=1
    
def plot_2_one(x=0):
    ''' Does what plot_2() does, but selecting only one value for Smax-Smin=x
    Default is x=0, i.e. the autocorrelation. plot_1() does the autocorrelation
    as well, with different colors representing different radii. For x=0 the
    result should coincide with plot_1'''
    plt.figure()
    ax = plt.gca()
    s = x

    j = 1
    x_data1 = np.arange(s,num_bins,1)
    y_data1 = np.zeros(num_bins-s)    
    for ii in np.arange(len(x_data1)):
        i = x_data1[ii]
        y_data1[ii] = g_result[i,i-s,j]
    plt.plot(x_data1,y_data1,'-',lw=0.75,c=cs[s],label=str(s))
    
    j = 2
    x_data1 = np.arange(s,num_bins,1)
    y_data1 = np.zeros(num_bins-s)    
    for ii in np.arange(len(x_data1)):
        i = x_data1[ii]
        y_data1[ii] = g_result[i,i-s,j]
    plt.plot(x_data1,y_data1,'--',lw=0.75,c=cs[s])

    j = 3
    x_data1 = np.arange(s,num_bins,1)
    y_data1 = np.zeros(num_bins-s)    
    for ii in np.arange(len(x_data1)):
        i = x_data1[ii]
        y_data1[ii] = g_result[i,i-s,j]
    plt.plot(x_data1,y_data1,':',lw=1,c=cs[s])

    j = 4
    x_data1 = np.arange(s,num_bins,1)
    y_data1 = np.zeros(num_bins-s)    
    for ii in np.arange(len(x_data1)):
        i = x_data1[ii]
        y_data1[ii] = g_result[i,i-s,j]
    plt.plot(x_data1,y_data1,'-.',lw=1,c=cs[s])

    j = 0
    x_data1 = np.arange(s,num_bins,1)
    y_data1 = np.zeros(num_bins-s)    
    for ii in np.arange(len(x_data1)):
        i = x_data1[ii]
        y_data1[ii] = g_result[i,i-s,j]
    plt.plot(x_data1,y_data1,'k-',lw=0.75)
    
    textstr2 = "separation radii"+'\n'+ "solid black="+str(r[0])+\
               '\n'+'solid='+str(r[1])+'\n'+'dashed='+str(r[2])+ \
               '\n'+"dotted="+str(r[3])+'\n'+"dash-dots="+str(r[4])
    ax.text(0.65,0.95,textstr2,transform=ax.transAxes, verticalalignment='top')
    plt.xlabel('Smax, largest grain size bin')
    plt.ylabel('RDF')
    #plt.yscale('log')
    plt.ylim([0,5])
    plt.xlim([0,9.1])
    plt.xticks(np.arange(num_bins))
    plt.legend(title="Smax - Smin=",loc=2)
    plt.title(title)
    ax.axhline(y=1,ls='dashed',lw=1,c='k') # horizontal dashed line to mark g=1

def plot_3():
    #WR plot
    plt.figure()
    ax = plt.gca()
    for s in [0,1,2,3,4]:
        j = 4
        x_data1 = np.arange(s,num_bins,1)
        y_data1 = np.zeros(num_bins-s)    
        for ii in np.arange(len(x_data1)):
            i = x_data1[ii]
            y_data1[ii] = wr_result[i,i-s,j]
            print(y_data1)
        plt.plot(x_data1,y_data1,'-',lw=1,c=cs[s],label=str(s))

    textstr2 = "separation radius="+str(r[j])
    ax.text(0.55,0.95,textstr2,transform=ax.transAxes, verticalalignment='top')
    
    plt.xlabel('Smax, largest grain size bin')
    plt.ylabel('wr(r)')
    plt.yscale('log')
    #plt.ylim([1,100])
    plt.xlim([0,9.1])
    plt.xticks(np.arange(num_bins))
    plt.legend(title="Smax - Smin=",loc=2)
    plt.title(title)

# Call the plot routines
plot_1()
#plot_2()
plot_2_one()
#plot_3()





