# Kinematic SLCM
Deterministic and Stochastic Superdroplet Models in a 2D kinematic flow

This simulation code in Fotran uses a stochastic Lagrangian (particle-based) cloud model (SLCM) to simulate warm and mixed-phase clouds. 
The microphysics is currently passive and running on a prescribed flow provided via an HDF5 input file with the velocity field data.
Parallel computing is implemented with OpenMP.

## Compilation
  Simply run the "make" command on the project's root folder. 
  
## Running
  Execute the run.sh script with two input parameters. The first one is the input file (examples in the Input/ folder) and the 
  other one is the number of threads. 
  Example: 
      ./run.sh input.in 16
      
## Results
  Output data is recorded into hierarchy data files (HDF5) with .h5 extension. You will need a library or a software tool to handle them and extract the data.
  I personally use Matlab, but feel free to choose whatever suits you best. 
  For more information, please refer to https://www.hdfgroup.org/solutions/hdf5/
      
