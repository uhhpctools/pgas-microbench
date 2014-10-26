FC = caft
LD = $(FC)
FFLAGS = 
LDFLAGS = $(FFLAGS)
GASNETRUN = gasnetrun_ibv

caf2.0_teams_microbenchmarks: caf2.0_teams_microbenchmarks.caf my_timer.o
		$(FC) $(LDFLAGS) $^ -o  $@ 

run:caf2.0_teams_microbenchmarks
	$(GASNETRUN) -np 4 ./$<

clean:
		rm -rf *.o caf2.0_teams_microbenchmarks *.mod *.rmod *.f90 *.a
