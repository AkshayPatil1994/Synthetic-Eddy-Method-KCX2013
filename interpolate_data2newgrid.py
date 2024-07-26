import numpy as np
import matplotlib.pyplot as plt

# First read the fine data
covar_org = np.loadtxt('covar_org.dat',skiprows=1)
mean_org = np.loadtxt('mean_org.dat',skiprows=1)
# Load the new grid to interpolate fine data
zgrid = np.loadtxt('grid_re350.dat')
# Define the data structure to save the values
mean = np.zeros((len(zgrid),np.shape(mean_org)[1]))
covar = np.zeros((len(zgrid),np.shape(covar_org)[1]))
# Interpolate all data to coarse grid
for i in range(2,np.shape(mean_org)[1]):
    mean[:,i] = np.interp(zgrid[:,1]/zgrid[-1,1],mean_org[:,1]/1000.0,mean_org[:,i])
for j in range(2,np.shape(covar_org)[1]):
    covar[:,j] = np.interp(zgrid[:,1]/zgrid[-1,1],covar_org[:,1]/1000.0,covar_org[:,j])
# Setup the z axis 
mean[:,1] = zgrid[:,1]
covar[:,1] = zgrid[:,1]
# Save data to file
np.savetxt('covar_re350.dat',covar[1:-1,:])
np.savetxt('mean_re350.dat',mean[1:-1,:])
