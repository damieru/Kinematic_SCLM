clear
rm ksdh5
rm *mod
rm *o
h5fc -pg -fdefault-real-8 -fdefault-double-8 -fopenmp -O3 variables.f90 functions.f90 -o ksdh5.out kinematic_SD.f90 mpdata_2d.f90

