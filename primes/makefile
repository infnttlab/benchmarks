HOSTNAME=$(shell hostname)
GCCVERSION=$(shell gcc -dumpversion)
CC=g++
CFLAGS=-std=c++11

all: compile
	mkdir -p bin
	mv primes_omp_$(HOSTNAME)_$(GCCVERSION) primes_omp_$(HOSTNAME)_$(GCCVERSION)_O0 primes_omp_$(HOSTNAME)_$(GCCVERSION)_O1 primes_omp_$(HOSTNAME)_$(GCCVERSION)_O2 primes_omp_$(HOSTNAME)_$(GCCVERSION)_O3 primes_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast bin

compile: primes_omp primes_omp_O0 primes_omp_O1 primes_omp_O2 primes_omp_O3 primes_omp_Ofast

primes_omp: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -o primes_omp_$(HOSTNAME)_$(GCCVERSION)

primes_omp_O0: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -O0 -o primes_omp_$(HOSTNAME)_$(GCCVERSION)_O0

primes_omp_O1: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -O1 -o primes_omp_$(HOSTNAME)_$(GCCVERSION)_O1

primes_omp_O2: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -O2 -o primes_omp_$(HOSTNAME)_$(GCCVERSION)_O2

primes_omp_O3: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -O3 -o primes_omp_$(HOSTNAME)_$(GCCVERSION)_O3

primes_omp_Ofast: primes_omp.cpp
	$(CC) $(CFLAGS) -fopenmp primes_omp.cpp -Ofast -o primes_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast

clean:
	rm -fr bin
