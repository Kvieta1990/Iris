OBJ	= main.o fileio.o pre_feff.o read_xmu.o clc_average.o dev_check.o remove_dups.o
FC	= gfortran
EXEC	= MD2FEFF
FLAGS	= -fbounds-check -o
LIB	= -L /data/home/apw247/Packages/Fortran/lib -linterp

$(EXEC): $(OBJ)
	$(FC) -fopenmp $(FLAGS) $@ $^ $(LIB)

main.o: main.f90
	$(FC) -fopenmp -fbounds-check -c $<

fileio.o: fileio.f90
	$(FC) -fbounds-check -c $<

pre_feff.o: pre_feff.f90 remove_dups.f90
	$(FC) -fbounds-check -c $^

read_xmu.o: read_xmu.f90
	$(FC) -fbounds-check -c $<

clc_average.o: clc_average.f90
	$(FC) -fbounds-check -c $<

dev_check.o: dev_check.f90
	$(FC) -fbounds-check -c $<

clean:
	rm -rf *.o $(EXEC)