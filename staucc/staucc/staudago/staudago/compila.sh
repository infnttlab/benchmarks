#!/bin/bash
gcc -Wall -g -o staudago -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
