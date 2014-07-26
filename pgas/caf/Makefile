FC = uhcaf
CC = uhcc
#FC = ftn
LD = $(FC)
FFLAGS =
LDFLAGS = $(FFLAGS)

caf_tests: caf_tests.o wtime.o
	$(LD) $(LDFLAGS) -o  $@ $^

caf_tests.o: caf_tests.F90
	$(FC) -c $(FFLAGS) $< -o $@

.PHONY:clean
clean:
	rm -rf *.o caf_tests *.mod
