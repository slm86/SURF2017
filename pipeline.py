''' Pipeline to Wheeler '''
import sys
sys.path.append('C:\Users\Stefi\Desktop\pfh_python-routines\grains')
import paramiko
from scp import SCPClient
import subprocess

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;

def grab_wheeler(sdir='.',snum=0,snap_type=0,destination='./scratch',other_filename=''):
    ''' snap_type=0 grabs the snapshot produced by GIZMO
        snap_type=1 grabs the snapshot created with grain_density_from_sanpshot.py
        other_filename if specified, grabs that file'''
    ssh = paramiko.SSHClient()

    passphrase='0206Teodora205'
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('wheeler.caltech.edu', username='slm86', password=passphrase, \
                 key_filename='C:\Users\Stefi\AppData\Local\lxss\home\slm86\.ssh\id_rsa')

    scp = SCPClient(ssh.get_transport())
    
    if(snap_type==0):
        ss=snap_ext(snum,four_char=0)
        filename=sdir+'/snapshot_'+ss+'.hdf5'     
    else:
        ss=snap_ext(snum,four_char=1)
        filename=sdir+'/dust_snap_'+ss+'.h5'
    
    if(len(other_filename)>0): filename=sdir + '/' + other_filename
    scp.get(filename,destination)
    print('File ',snum,' copied')
    ssh.close()

def scratch(): 
    ''' Clears the /scratch folder on local machine '''
    # print(subprocess.check_output('dir',shell=True))
    subprocess.call('del scratch\dust*',shell=True)
    subprocess.call('del scratch\snap*',shell=True)
    subprocess.call('del .\scratch\rdf*',shell=True)
    print('Scratch folder empty')

### USE THIS TO COPY OVER FILES QUICKLY ###
#for i in np.arange(0,30):
#    grab_wheeler(sdir='/panfs/ds08/hopkins/slm86/c.k05N4.2/output/',\
#                 snum=i,snap_type=1,destination='../c.k05N4.2/')
    
grab_wheeler(sdir='/panfs/ds08/sxs/jsquire/gizmo-dust/s.kneg05N4-M1.5/output',snum=0,snap_type=0,\
             destination='./scratch',other_filename='rdf_wr_snap_0450.h5')