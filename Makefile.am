CC= $(MPICC)
bin_PROGRAMS = TauBench 
TauBench_SOURCES = flux_fake.c main_fake.c lim_fake.c smooth_fake.c \
      flop.h lim_fake.h nodep.h types.h \
      expand.h flux_fake.h main_fake.h smooth_fake.h util.h
EXTRA_DIST = SUPPORTED_MPI_VERSIONS SUPPORTED_BATCH_SYSTEMS RESULTS
run: TauBench
	./submit -v -bm TauBench
