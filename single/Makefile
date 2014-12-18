all: wave_single

wave_single: wave_single.f90
	pgfortran -fast -Mipa=fast,inline wave_single.f90 -o wave_single

clean:
	rm wave_single *.o
