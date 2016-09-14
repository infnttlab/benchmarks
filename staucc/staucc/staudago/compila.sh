#!/bin/bash
#icc -Wall -g -o staudago -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
#Advisor
icc -Wall -g -O2 -fno-inline-functions -o staudago -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm 

#nedit main.c manage_prob.c io_file.c p_sim.h manage_prob.h io_file.h &
