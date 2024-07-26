import numpy as np
import matplotlib.pyplot as plt

# -
#
# SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
def read_restart_file(filenamei,ng):
    import numpy as np
    offset     = 0
    iprecision = 8
    disp       = np.prod(ng)
    fldinfo    = np.zeros([2])
    data       = np.zeros([ng[0],ng[1],ng[2],4])
    with open(filenamei,'rb') as f:
        for p in range(4):
            f.seek(offset)
            fld = np.fromfile(f,dtype=np.float64,count=disp)
            data[:,:,:,p] = np.reshape(fld,(ng[0],ng[1],ng[2]),order='F')
            offset += iprecision*disp
        f.seek(offset)
        fldinfo[:] = np.fromfile(f,dtype=np.float64,count=2)
    f.close
    u=data[:,:,:,0]
    v=data[:,:,:,1]
    w=data[:,:,:,2]
    p=data[:,:,:,3]
    time  =     fldinfo[0]
    istep = int(fldinfo[1])

    return u, v, w, p, time, istep
# START THE MAIN PROGRAM
N = [1048,764,128]
targetRe = 350.0
targetH = 1.0
targetNu = 1.0e-6
utau = targetRe*targetNu/targetH
# Read the binary snapshot
[u,v,w,p,time,istep] = read_restart_file('fld.bin',N)
# Load reference grid
umean = np.loadtxt('mean_re350.dat',skiprows=1)
z = umean[:,1]
# Compute some quantities to verify
uprime = u - np.mean(u,axis=(0,1))
stress = np.mean(w*uprime,axis=(0,1)) - targetNu*np.gradient(np.mean(u,axis=(0,1)),z)
# Plotting
plt.figure(1)
plt.subplot(1,3,1)
plt.pcolor(u[:,:,60]/utau,cmap='turbo')
plt.subplot(1,3,2)
plt.pcolor(v[:,:,60]/utau,cmap='turbo')
plt.subplot(1,3,3)
plt.pcolor(w[:,:,60]/utau,cmap='turbo')


plt.figure(2)
plt.plot(np.sqrt(np.mean(uprime**2,axis=(0,1)))/utau,z*targetRe/targetH,'r')
plt.plot(np.sqrt(np.mean(v**2,axis=(0,1)))/utau,z*targetRe/targetH,'b')
plt.plot(np.sqrt(np.mean(w**2,axis=(0,1)))/utau,z*targetRe/targetH,'k')
plt.plot(stress/utau**2,z*targetRe/targetH,'mo',markerfacecolor='None')
plt.plot(np.mean(w*uprime,axis=(0,1))/utau**2,z*targetRe/targetH,'g')
plt.plot(z/targetH-1,z*targetRe/targetH,'--')
plt.xlim([-1,3.0])
plt.show()
