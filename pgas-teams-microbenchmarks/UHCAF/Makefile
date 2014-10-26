FC = uhcaf
LD = $(FC)
FFLAGS = 
LDFLAGS = $(FFLAGS)
CAFRUN = cafrun

caf_teams_microbenchmarks: caf_teams_microbenchmarks.o my_timer.o
		$(LD) $(LDFLAGS) -o  $@ $^

caf_teams_microbenchmarks.o: caf_teams_microbenchmarks.f90
		$(FC) -c $(FFLAGS) $< -o $@

run:caf_teams_microbenchmarks
	$(CAFRUN) -np 4 ./$<

clean:
		rm -rf *.o caf_teams_microbenchmarks *.mod
