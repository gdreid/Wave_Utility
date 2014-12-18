all: wave_oacc

wave_oacc: wave_oacc.f90
	pgfortran -acc -ta=nvidia -Minfo=accel wave_oacc.f90 -o wave_oacc

clean:
	rm wave_oacc *.o
