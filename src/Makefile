# make file for JML

OBJS = jml.o MCMCcoal.o seq-gen.o global.o evolve.o model.o nucmodels.o\
	aamodels.o eigen.o treefile.o progress.o gamma.o twister.o

LIBS = -lstdc++ -lm
#CC = cc -g
CC = gcc

jml: $(OBJS)
	$(CC) -o jml $(OBJS) $(LIBS)

jml.o: jml.h MCMCcoal.c seq-gen.c jml.cpp MCMCcoal.h seq-gen.h treeinfo.h nexusdata.h GenFunctions.h
	$(CC) -o jml.o -c jml.cpp

MCMCcoal.o: MCMCcoal.h MCMCcoal.c
	$(CC) -o MCMCcoal.o -c MCMCcoal.c

seq-gen.o: seq-gen.c seq-gen.h global.h treefile.h evolve.h model.h nucmodels.h aamodels.h progress.h twister.h
	$(CC) -o seq-gen.o -c seq-gen.c

aamodels.o: aamodels.c aamodels.h
	$(CC) -o aamodels.o -c aamodels.c

eigen.o: eigen.c eigen.h
	$(CC) -o eigen.o -c eigen.c

evolve.o: evolve.c evolve.h global.h tree.h
	$(CC) -o evolve.o -c evolve.c

gamma.o: gamma.c gamma.h
	$(CC) -o gamma.o -c gamma.c

global.o: global.c global.h
	$(CC) -o global.o -c global.c

model.o: model.c model.h
	$(CC) -o model.o -c model.c

nucmodels.o: nucmodels.c nucmodels.h
	$(CC) -o nucmodels.o -c nucmodels.c

progress.o: progress.c progress.h
	$(CC) -o progress.o -c progress.c

treefile.o: treefile.c treefile.h
	$(CC) -o treefile.o -c treefile.c

twister.o: twister.c twister.h
	$(CC) -o twister.o -c twister.c
	
clean:
	rm -rf *.o *~

