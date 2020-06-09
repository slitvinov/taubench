.POSIX:
.PHONY: clean
CC = cc
CFLAGS = -g -Ofast
LINK = $(CC)

M = \
flux_fake.o\
main_fake.o\
lim_fake.o\
smooth_fake.o\

taubench: $M
	$(LINK) $(LDFLAGS) -lm -o $@ $M
.c.o:
	$(CC) $(CFLAGS) $(CO_CFLAGS) $< -c

clean:
	rm -f $M taubench

include dep.mk
