flux.o: flux.c util.h types.h main.h flux.h flop.h expand.h nodep.h
lim.o: lim.c types.h util.h main.h lim.h flop.h nodep.h expand.h
main.o: main.c types.h lim.h flux.h smooth.h util.h flop.h
smooth.o: smooth.c types.h main.h smooth.h flop.h util.h nodep.h
