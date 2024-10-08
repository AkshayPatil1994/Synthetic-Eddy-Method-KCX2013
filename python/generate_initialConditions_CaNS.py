import numpy as np
import matplotlib.pyplot as plt
import sys 

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
        time.tofile(f)
        iteration.tofile(f)
    f.close()
# START MAIN PROGRAM
writeRestart = True
checkPlot = False
zloc = 40
nslices = 1048
Ny = 764
Nz = 128
# Define the datastructure
udata = np.zeros((nslices,Ny,Nz))
vdata = np.zeros((nslices,Ny,Nz))
wdata = np.zeros((nslices,Ny,Nz))
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
    plt.pcolor((uprime[:,:,zloc].T)*5.0e-2,cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.pcolor((vdata[:,:,zloc].T)*5.0e-2,cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,3)
    plt.pcolor((wdata[:,:,zloc].T)*5.0e-2,cmap='turbo')
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.plot(np.sqrt(np.nanmean(uprime**2,axis=(0,1))),z)
    plt.plot(np.sqrt(np.nanmean(vdata**2,axis=(0,1))),z)
    plt.plot(np.sqrt(np.nanmean(wdata**2,axis=(0,1))),z)

    plt.figure(2)
    plt.semilogx(z*1000.0/1.0,np.nanmean(udata,axis=(0,1)))
    plt.semilogx(z[1:10]*1000.0/1.0,z[1:10]*1000.0/1.0,'r--')
    plt.show()
