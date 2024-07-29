import numpy as np
import matplotlib.pyplot as plt

# Load reference data
covar = np.loadtxt('covar_re350.dat',skiprows=1)
mean = np.loadtxt('mean_re350.dat',skiprows=1)
ds=np.shape(covar)
# Assign the right values
dataout = np.zeros((ds[0],11))
print(np.shape(dataout))
dataout[:,0] = mean[:,1]	# Z
dataout[:,1] = mean[:,2]	# Umean
dataout[:,2] = covar[:,2]	# u'u'
dataout[:,3] = covar[:,6]	# u'v'
dataout[:,4] = covar[:,5]	# u'w'
dataout[:,5] = covar[:,6]	# v'u'
dataout[:,6] = covar[:,4]	# v'v'
dataout[:,7] = covar[:,7]	# v'w'
dataout[:,8] = covar[:,5]	# w'u'
dataout[:,9] = covar[:,7]	# w'v'
dataout[:,10] = covar[:,3]	# w'w'
# Write data to file
np.savetxt('inflow.dat',dataout)
