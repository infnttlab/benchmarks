#!/bin/bash
/opt/openmpi14/bin/mpicc -Wall -g -o pardago -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
