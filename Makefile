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

include dep.mk
