.POSIX:
.SUFFIXES:
.SUFFIXES: .c
.SUFFIXES: .o

MPICC = mpicc
LINK = $(MPICC)
CFLAGS = -g -O2

M = \
flux.o\
main.o\
lim.o\
smooth.o\

taubench: $M
	$(LINK) $M $(LDFLAGS) -lm -o $@
.c.o:
	$(MPICC) $(CFLAGS) $< -c
clean:
	rm -f $M taubench

flux.o: flux.c util.h types.h main.h flux.h flop.h expand.inc nodep.inc
lim.o: lim.c types.h util.h main.h lim.h flop.h nodep.inc expand.inc
main.o: main.c types.h main.h lim.h flux.h smooth.h util.h flop.h
smooth.o: smooth.c types.h main.h smooth.h flop.h util.h nodep.inc
