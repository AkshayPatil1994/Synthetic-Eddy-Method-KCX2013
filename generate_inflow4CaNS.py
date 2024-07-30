import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# USER INPUT
plotData = True
saveData = False
saveSep = False
# FLOW SCALING PARAMETERS
Re = 392.24
H = 1.0
kvisc = 1e-6
# Name of the output files
utau = Re*kvisc/H
inflowOutFilename='inflow_new.dat'
meanoutfile = 'mean_interpolated_re395.dat'
covaroutfile = 'covar_interpolated_re395.dat'
# First read the fine data
covar_org = np.loadtxt('covar_re395.dat',skiprows=1)
mean_org = np.loadtxt('mean_re395.dat',skiprows=1)
# Load the new grid to interpolate fine data
zgrid = np.loadtxt('grid.out')
myz = zgrid[:,1]         # Set the right z location
# Define the data structure to save the values
mean = np.zeros((len(myz),np.shape(mean_org)[1]))
covar = np.zeros((len(myz),np.shape(covar_org)[1]))
# Interpolate all data to coarse grid
for i in range(2,np.shape(mean_org)[1]):
    mean[:,i] = np.interp(myz,mean_org[:,1]*(1.0/(Re/H)),mean_org[:,i]*utau)

for j in range(2,np.shape(covar_org)[1]):
    covar[:,j] = np.interp(myz,mean_org[:,1]*(1.0/(Re/H)),covar_org[:,j]*utau**2)
# Setup the z axis
mean[:,1] = myz
covar[:,1] = myz
# Plotting check
if(plotData):
    plt.figure(1)
    plt.title('Mean velocity check')
    plt.plot(mean_org[:,2]*utau,mean_org[:,1]*(1.0/(Re/H)),'kd',label='reference',markerfacecolor='None')
    plt.plot(mean[:,2],mean[:,1],'r+',label='coarse',markerfacecolor='None')
    plt.legend()
    plt.figure(2)
    plt.title('Covariance check')
    for j in range(2,np.shape(covar_org)[1]):
        plt.plot(covar_org[:,j]*utau**2,covar_org[:,1]*(1.0/(Re/H)),'kd',label='reference',markerfacecolor='None')
        plt.plot(covar[:,j],covar[:,1],'r+',label='coarse',markerfacecolor='None')
        if(j==2):
            plt.legend()
    plt.show()
# Save data to file
if(saveData):
    ds = np.shape(covar)
    # Assign the right values
    dataout = np.zeros((ds[0], 11))
    dataout[:,0] = mean[:,1]    # Z
    dataout[:,1] = mean[:,2]    # Umean
    dataout[:,2] = covar[:,2]   # u'u'
    dataout[:,3] = covar[:,6]   # u'v'
    dataout[:,4] = covar[:,5]   # u'w'
    dataout[:,5] = covar[:,6]   # v'u'
    dataout[:,6] = covar[:,4]   # v'v'
    dataout[:,7] = covar[:,7]   # v'w'
    dataout[:,8] = covar[:,5]   # w'u'
    dataout[:,9] = covar[:,7]   # w'v'
    dataout[:,10] = covar[:,3]  # w'w'
    # Write data to file
    myheader = ['z', 'Umean', 'Ruu', 'Ruv', 'Ruw', 'Rvu', 'Rvv', 'Rvw', 'Rwu', 'Rwv', 'Rww']
    df = pd.DataFrame(dataout, columns=myheader)
    # Save dataframe as double precision csv with whitespace
    df.to_csv(inflowOutFilename, index=False, header=True, sep=' ', float_format='%.15e')
# Save separate files for generating initial conditions
if(saveSep):
    mymeanheader = ['y','y+','Umean','dUmean/dy' ,'Wmean','dWmean/dy' ,'Pmean']
    mycovarheader = ['y','y+','R_uu' ,'R_vv' ,'R_ww' ,'R_uv' ,'R_uw' ,'R_vw']
    df2 = pd.DataFrame(mean, columns=mymeanheader)
    df2.to_csv(meanoutfile,index=False, header=True, sep=' ', float_format='%.15e')
    df3 = pd.DataFrame(covar, columns=mycovarheader)
    df3.to_csv(covaroutfile,index=False, header=True, sep=' ', float_format='%.15e')
