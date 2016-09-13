HOSTNAME=$(shell hostname)
GCCVERSION=$(shell gcc -dumpversion)
CC=g++
CFLAGS=-std=c++11

all: compile
	mkdir -p bin
	mv pi_loop_omp_$(HOSTNAME)_$(GCCVERSION) pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O0 pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O1 pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O2 pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O3 pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION) pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O0 pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O1 pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O2 pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O3 pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast bin

compile: pi_loop_omp pi_loop_omp_O0 pi_loop_omp_O1 pi_loop_omp_O2 pi_loop_omp_O3 pi_loop_omp_Ofast pi_loop_no_omp pi_loop_no_omp_O0 pi_loop_no_omp_O1 pi_loop_no_omp_O2 pi_loop_no_omp_O3 pi_loop_no_omp_Ofast


pi_loop_omp: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)

pi_loop_omp_O0: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp -O0 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O0

pi_loop_omp_O1: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp -O1 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O1

pi_loop_omp_O2: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp -O2 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O2

pi_loop_omp_O3: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp -O3 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_O3

pi_loop_omp_Ofast: pi_loop_omp.cpp
	$(CC) $(CFLAGS) -fopenmp -Ofast pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast

pi_loop_no_omp: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)

pi_loop_no_omp_O0: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) pi_loop_no_omp.cpp -O0 -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O0

pi_loop_no_omp_O1: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) -O1 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O1

pi_loop_no_omp_O2: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) -O2 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O2

pi_loop_no_omp_O3: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) -O3 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_O3

pi_loop_no_omp_Ofast: pi_loop_no_omp.cpp
	$(CC) $(CFLAGS) -O3 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(GCCVERSION)_Ofast

clean:
	rm -fr bin

