FC = mpif90 
LD = $(FC)
FFLAGS = 
LDFLAGS = $(FFLAGS)
MPIRUN = mpirun

mpi_teams_microbenchmarks: mpi_teams_microbenchmarks.o 
		$(LD) $(LDFLAGS) -o  $@ $^

mpi_teams_microbenchmarks.o: mpi_teams_microbenchmarks.f90
		$(FC) -c $(FFLAGS) $< -o $@

run:mpi_teams_microbenchmarks
	$(MPIRUN) -np 4 ./$<

clean:
		rm -rf *.o mpi_teams_microbenchmarks *.mod
