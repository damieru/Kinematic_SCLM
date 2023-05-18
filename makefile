CC=h5fc
CFLAGS=-fopenmp -O3 -J$(ODIR) -fdefault-real-8 -fdefault-double-8

ODIR = src/obj
SDIR = src
_OBJ = profiler.o singleton.o numerical_methods.o variables.o functions.o HDF5_functions.o kinematic_SD.o mpdata_2d.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

ksd: $(OBJ)  
	$(CC) -o $@.out $^ $(CFLAGS)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(CC) -c -o $@ $^ $(CFLAGS)

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod *.out

