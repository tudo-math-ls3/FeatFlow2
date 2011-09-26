shallow_ifort: sw.o vartypes.o exp_sw_schemes.o imp_sw_schemes.o sw_flux.o
	ifort -fast vartypes.o exp_sw_schemes.o imp_sw_schemes.o sw_flux.o sw.o -o sw

vartypes.o: vartypes.f90
	ifort -fast -c vartypes.f90

exp_sw_schemes.o: exp_sw_schemes.f90\
	vartypes.o
	ifort -fast -c exp_sw_schemes.f90

imp_sw_schemes.o: imp_sw_schemes.f90\
	vartypes.o
	ifort -fast -c imp_sw_schemes.f90

sw_flux.o: sw_flux.f90\
	vartypes.o
	ifort -fast -c sw_flux.f90

sw.o: sw.f90\
	vartypes.o exp_sw_schemes.o imp_sw_schemes.o sw_flux.o
	ifort -fast -c sw.f90


#shallow_g95: sw.f90
#	g95 -O3 vartypes.f90 exp_sw_schemes.f90 imp_sw_schemes.f90 sw_flux.f90 sw.f90 -o sw
	


#shallow_ifort: sw.f90
#	ifort -fast vartypes.f90 exp_sw_schemes.f90 imp_sw_schemes.f90 sw_flux.f90 sw.f90 -o sw
