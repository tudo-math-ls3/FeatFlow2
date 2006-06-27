LAPACK=/home/user/mmoelle1/nobackup/Featflow2/object/libraries/lib-pc-pentium4-linux-ifc/liblapack.a
BLAS=/home/user/mmoelle1/work/codes/libgoto/libgoto_opteron-r1.02.a

all:
	/usr/local/ifc-9.0-32/bin/ifort -O3 -xW -cm -fpe0 -vec_report0 -ipo -c invert.f90
	/usr/local/ifc-9.0-32/bin/ifort -O3 -xW -cm -fpe0 -vec_report0 -ipo -c test_invert.f90
	/usr/local/ifc-9.0-32/bin/ifort -O3 -xW -cm -fpe0 -vec_report0 -ipo -c invertold.f
	/usr/local/ifc-9.0-32/bin/ifort -O3 -xW -cm -fpe0 -vec_report0 -ipo test_invert.o invert.o invertold.f $(LAPACK) $(BLAS) -o test_invert

prof:
	/usr/local/ifc-9.0-32/bin/ifort -p -O3 -xW -cm -fpe0 -vec_report0 -ipo -c invert.f90
	/usr/local/ifc-9.0-32/bin/ifort -p -O3 -xW -cm -fpe0 -vec_report0 -ipo -c test_invert.f90
	/usr/local/ifc-9.0-32/bin/ifort -p -O3 -xW -cm -fpe0 -vec_report0 -ipo -c invertold.f
	/usr/local/ifc-9.0-32/bin/ifort -p -O3 -xW -cm -fpe0 -vec_report0 -ipo test_invert.o invert.o invertold.f $(LAPACK) $(BLAS) -o test_invert

debug:
	ifort -warn -check -c invert.f90
	ifort -warn -check -c test_invert.f90
	ifort -warn -check -c invertold.f
	ifort -warn -check test_invert.o invert.o invertold.f $(LAPACK) $(BLAS) -o test_invert

