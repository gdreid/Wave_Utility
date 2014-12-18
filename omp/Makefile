all: wave_omp

wave_omp: wave_omp.f90
	pgfortran -mp -fast wave_omp.f90 -o wave_omp

clean:
	rm wave_omp *.o
