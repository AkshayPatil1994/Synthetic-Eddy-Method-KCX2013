# Serial Implementation of the Synthetic Eddy Method based on Kim, Castro, Xie (2013)

## How to compile the code
```
cd src
make
```

## How to change parameters
Follow the instructions in `parameters.in` to change the input parameters

## Setup the `inflow.dat` input file
A reference `inflow.dat` has been listed in the `python/inflow.dat` directory of the code. The file contains the **dimensional** vertical coordinate `z` or `$x_3$` and **non-dimensional** mean velocity and Reynolds stress tensor components used by the code to generate the initial conditions. 

1. Be sure to have the profiles collocated with the vertical grid used in the simulation
2. The first line of the file is ignored as it is the header file containing the information about the column data

How to obtain the `inflow.dat` file from scratch?

1. Download the reference data from a datasource of your choice. For standard datasets for channel flows, the Texas university fileserver is a good options (https://turbulence.oden.utexas.edu/)
2. Be sure to download the mean profiles and the covariance profiles
3. Consult the python script `interpolate_data2newgrid.py` to see how the data is arranged. 
**WARNING** This code assumes that the z direction is vertical!
4. Use the `interpolate_data2newgrid.py` script to interpolate the data on your simulation grid.
5. Once that is completed, you can use the `inflow.dat` to run the fortran code that generates individual slices.
6. These slices can then be stitched together using the `generate_inflow4CaNS.py` script which will output the `fld.bin` file used in CaNS as a checkpoint. Note that the checkpoint starts at `t = 0`.

## How to run the code
Copy the `parameters.in` file to the same location as the executable `genic`
```
# Copy the parameters file (assuming you are in src directory already)
cp ../parameters.in src/   
# Run the compiled program 
./genic
```
Example output

```
  ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓██████▓▒░  
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
 ░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░        
 ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
  ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓██████▓▒░  
 -------------------------------------------------------
 Utau input:   6.0000000000000001E-003
 Integral length scale used   1.6666666666666666E-002
 Input bulk velocity is  0.12155821700951484     
 Lagrangian time scale:   0.13710851538207491     
 Time step is    2.0000000000000000E-002
 -------------------------------------------------------
 WARNING: The digital filter has isotropic kernel both in y and z!
 Value of alpha is:  0.89175315170062130       0.79522368356799134     
||                                                               |  1.60%   16/1000
```

## How to visualise the results
Use one of the python scripts to load the dataset and plot slices
