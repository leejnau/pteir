# Makefile for OU MS Thesis code

#CFLAGS =-g -c -std=c99 -pg   # debug flags
#CFLAGS = -O2 -funroll-loops -c -std=c99 -pg -fopenmp   # profile flags
CFLAGS = -O3 -funroll-loops -c -std=c99 -fopenmp  # optimized flags

all: seqs.o divsufsort.o sa.o rt.o teir.o teir

teir: seqs.o divsufsort.o sa.o rt.o teir.o main.c
	#gcc $(CFLAGS) main.c
	#gcc *.o -g main.c -o main -fopenmp # debug
	#gcc -O2 -funroll-loops -pg *.o main.c -o main -std=c99 -fopenmp # profile
	gcc *.o -O3 -funroll-loops main.c -o teir -std=c99 -fopenmp  # optimized

teir.o: seqs.h seqs.c sa.h sa.c teir.h teir.c
	gcc $(CFLAGS) teir.c

rt.o: rt.h rt.c
	gcc $(CFLAGS) rt.c

sa.o: sa.h sa.c
	gcc $(CFLAGS) sa.c 

divsufsort.o: divsufsort.h divsufsort.c
	gcc $(CFLAGS) divsufsort.c

seqs.o: seqs.h seqs.c
	gcc $(CFLAGS) seqs.c

clean: 
	rm *.o teir
