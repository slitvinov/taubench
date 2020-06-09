flux_fake.o: flux_fake.c util.h types.h main_fake.h flux_fake.h flop.h \
 expand.h nodep.h
lim_fake.o: lim_fake.c types.h util.h main_fake.h lim_fake.h flop.h \
 nodep.h expand.h
main_fake.o: main_fake.c types.h lim_fake.h flux_fake.h smooth_fake.h \
 util.h flop.h
smooth_fake.o: smooth_fake.c types.h main_fake.h smooth_fake.h flop.h \
 util.h nodep.h
