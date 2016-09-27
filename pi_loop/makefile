#COMP=intel
COMP=gcc
HOSTNAME = $(shell hostname -s)
CXXVERSION = $(shell gcc -dumpversion)
CXXFLAGS = -std=c++11
OMPFLAG = -fopenmp

ifeq ($(COMP),intel)
CXX=icpc
CXXVERSION = $(shell icpc -dumpversion)
CXXFLAGS = -std=c++0x
OMPFLAG = -openmp
endif

#.PHONY: print_vars
#print_vars:
#	echo $(CXX)

all: compile
	mkdir -p bin
	mv pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION) pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O0 pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O1 pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O2 pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O3 pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_Ofast pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)  pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O0 pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O1 pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O2 pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O3 pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_Ofast bin         

compile: pi_loop_omp pi_loop_omp_O0 pi_loop_omp_O1 pi_loop_omp_O2 pi_loop_omp_O3 pi_loop_omp_Ofast pi_loop_no_omp pi_loop_no_omp_O0 pi_loop_no_omp_O1 pi_loop_no_omp_O2 pi_loop_no_omp_O3 pi_loop_no_omp_Ofast


pi_loop_omp: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)

pi_loop_omp_O0: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -O0 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O0

pi_loop_omp_O1: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -O1 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O1

pi_loop_omp_O2: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -O2 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O2

pi_loop_omp_O3: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -O3 pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O3

pi_loop_omp_Ofast: pi_loop_omp.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -Ofast pi_loop_omp.cpp -o pi_loop_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_Ofast

pi_loop_no_omp: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)

pi_loop_no_omp_O0: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) pi_loop_no_omp.cpp -O0 -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O0

pi_loop_no_omp_O1: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) -O1 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O1

pi_loop_no_omp_O2: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) -O2 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O2

pi_loop_no_omp_O3: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) -O3 pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_O3

pi_loop_no_omp_Ofast: pi_loop_no_omp.cpp
	$(CXX) $(CXXFLAGS) -Ofast pi_loop_no_omp.cpp -o pi_loop_no_omp_$(HOSTNAME)_$(CXX)_$(CXXVERSION)_Ofast

clean:
	rm -fr bin

