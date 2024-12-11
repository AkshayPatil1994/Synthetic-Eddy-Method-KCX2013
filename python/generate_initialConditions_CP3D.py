import numpy as np
import matplotlib.pyplot as plt
import sys
#
# This function writes the restart file for CP3D
#
def n2carray(udata,vdata,wdata,fileout,time=np.array(0.0),iteration=np.array(1,dtype=np.int64),scaling=1.0):
    # Now flatten the array
    udata = scaling*udata.flatten(order='F')	        # Flatten it to make FORTRAN binary compatible
    vdata = scaling*vdata.flatten(order='F')	        # Flatten it to make FORTRAN binary compatible
    wdata = scaling*wdata.flatten(order='F')	        # Flatten it to make FORTRAN binary compatible
    pdata = np.zeros_like(wdata)
    # Open file and write to file
    with open(fileout, 'wb') as f:
        udata.tofile(f)
        vdata.tofile(f)
        wdata.tofile(f)
        pdata.tofile(f)
    f.close()
# START MAIN PROGRAM
writeRestart = True
checkPlot = False
zloc = 40
retau = 1000.0
H = 1.0
kvisc = 1.0e-6
nslices = 1048
Ny = 764
Nz = 128
# Define the datastructure
udata = np.zeros((nslices,Ny,Nz))
vdata = np.zeros((nslices,Ny,Nz))
wdata = np.zeros((nslices,Ny,Nz))
utau = retau*kvisc/H
# Load original z
umean = np.loadtxt('mean_re350.dat',skiprows=1)
z = umean[:,1]
# Specify the filename
for myslice in range(0,nslices):
    filename = 'slices/uslicedata_'+str(myslice+1)+'.dat'
    udata[myslice,:,:] = np.loadtxt(filename)
    filename = 'slices/vslicedata_'+str(myslice+1)+'.dat'
    vdata[myslice,:,:] = np.loadtxt(filename)
    filename = 'slices/wslicedata_'+str(myslice+1)+'.dat'
    wdata[myslice,:,:] = np.loadtxt(filename)

# Write the binary restart file
print("U inf:",np.any(np.isinf(udata) == True))
print("V inf:",np.any(np.isinf(vdata) == True))
print("W inf:",np.any(np.isinf(wdata) == True))
print("U nan:",np.any(np.isnan(udata) == True))
print("V nan:",np.any(np.isnan(vdata) == True))
print("W nan:",np.any(np.isnan(wdata) == True))

# Swap Y and Z directions
temp = np.transpose(udata,axis=(0,2,1))
del udata
udata = temp
temp = np.transpose(vdata,axis=(0,2,1))
del vdata
vdata = temp
temp = np.transpose(wdata,axis=(0,2,1))
del wdata
wdata = temp
# Write to file
if(writeRestart):
    n2carray(udata,vdata,wdata,'fld.bin')
# Plot if requested
if(checkPlot):
    plt.pcolor(udata[:,:,10])
    plt.colorbar()
    #Compute the u'
    uprime = udata - np.nanmean(udata,axis=(0,1))
    # Check a slice
    plt.subplot(2,2,1)
    plt.pcolor((uprime[:,:,zloc].T),cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.pcolor((vdata[:,:,zloc].T),cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,3)
    plt.pcolor((wdata[:,:,zloc].T),cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.plot(np.sqrt(np.nanmean(uprime**2,axis=(0,1))),z)
    plt.plot(np.sqrt(np.nanmean(vdata**2,axis=(0,1))),z)
    plt.plot(np.sqrt(np.nanmean(wdata**2,axis=(0,1))),z)

    plt.figure(2)
    plt.semilogx(z*retau/H,np.nanmean(udata,axis=(0,1)))
    plt.semilogx(z[1:10]*retau/H,z[1:10]*retau/H,'r--')
    plt.show()
