CC = cc
CFLAGS = -g -Ofast
LINK = $(CC)

M = \
flux_fake.o\
main_fake.o\
lim_fake.o\
smooth_fake.o\

TauBench: $M
	$(LINK) $(LDFLAGS) -lm -o $@ $M

.c.o:
	$(CC) $(CFLAGS) $(CO_CFLAGS) $< -c
