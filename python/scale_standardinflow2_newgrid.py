import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# USER INPUT
plotData = True
saveData = False
# First read the fine data
data_input = np.loadtxt('reference_inflow.dat',skiprows=1)
# Load the new grid to interpolate fine data
zgrid = np.loadtxt('grid.out')
myz = zgrid[1:-1,1]         # Set the right z location
# Define the data structure to save the values
inflowdata = np.zeros((len(myz),11))
inflowdata[:,0] = myz
# Interpolate all data to coarse grid
for i in range(1,11):
    inflowdata[:,i] = np.interp(myz/myz[-1],data_input[:,0]/data_input[-1,0],data_input[:,i])
# Plotting check
if(plotData):
    plt.figure(1)
    plt.plot(inflowdata[:,1],inflowdata[:,0]/inflowdata[-1,0],'ro')
    plt.plot(data_input[:,1],data_input[:,0]/data_input[-1,0],'k-')    
    # Plot rms
    plt.figure(2)
    plt.plot(inflowdata[:,2],inflowdata[:,0]/inflowdata[-1,0],'ro')
    plt.plot(inflowdata[:,3],inflowdata[:,0]/inflowdata[-1,0],'ro')
    plt.plot(inflowdata[:,4],inflowdata[:,0]/inflowdata[-1,0],'ro')
    plt.plot(data_input[:,2],data_input[:,0]/data_input[-1,0],'k-')
    plt.plot(data_input[:,3],data_input[:,0]/data_input[-1,0],'k-')
    plt.plot(data_input[:,4],data_input[:,0]/data_input[-1,0],'k-')
    plt.show()
# Save data to file
if(saveData):
    # Write data to file
    myheader = ['z', 'Umean', 'Ruu', 'Ruv', 'Ruw', 'Rvu', 'Rvv', 'Rvw', 'Rwu', 'Rwv', 'Rww']
    df = pd.DataFrame(inflowdata, columns=myheader)
    # Save dataframe as double precision csv with whitespace
    df.to_csv('inflow.dat', index=False, header=True, sep=' ', float_format='%.15e')
