GCC = gcc
ICCFLAGS =
#OPENMPICCFLAGS = -qopenmp
#PTHREADSICCFLAGS = -pthread

#GCC = gcc
#GCCFLAGS = -fcilkplus -O3
#OPENMPGCCCFLAGS = -fopenmp

CMAIN = ex3
#NPROCS = 4


all: main.o
	$(GCC) $(GCCFLAGS) $^ -o $(CMAIN)

main.o: main.c
	$(GCC) -c $(ICCFLAGS) $^


clean:
	rm -f *.o *~ $(CMAIN)
