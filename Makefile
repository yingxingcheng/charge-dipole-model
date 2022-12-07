
FC=gfortran
FFLAGS =  -O3 -llapack

install :
	$(FC) $(FFLAGS) helper.f90 gauss_charge_dipole_model.f90 -o main.x


clean:
	rm -f *.o  *.mod *.x
