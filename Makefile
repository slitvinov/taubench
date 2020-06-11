.POSIX:
.PHONY: clean
CC = cc
CFLAGS = -g -Ofast
LINK = $(CC)

M = \
flux.o\
main.o\
lim.o\
smooth.o\

taubench: $M
	$(LINK) $(LDFLAGS) -lm -o $@ $M
.c.o:
	$(CC) $(CFLAGS) $< -c

clean:
	rm -f $M taubench

include dep.mk
