# Serial Implementation of the Synthetic Eddy Method based on Kim, Castro, Xie (2013)

## How to compile the code
```
cd src
gfortran -O3 -fcheck=bounds utils.f90 main.f90 -o genflow
```

## How to change parameters
Follow the instructions in `parameters.in` to change the input parameters

## How to visualise the results
Use one of the python scripts to load the dataset and plot slices
