all: wave_mpi

wave_mpi.o: wave_mpi.f90
	mpif90 -c -Ofast wave_mpi.f90

wave_mpi: wave_mpi.o
	mpif90 -o wave_mpi wave_mpi.o -lmpi -Ofast 

clean:
	rm wave_mpi *.o
