#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>

#ifndef _P_SIM_H
   #include "p_sim.h"
   #define _P_SIM_H
#endif


void uliRows_uliColumns(FILE *, ULI *, ULI *);

unsigned **load_table(ULI *, ULI *, FILE *);
int **iLoad_table(FILE *, int, struct globale*);

void add_vector_in_matrix(LD *, LD *, ULI , LD, int, struct globale*);
void write_matrix(FILE * fPtr, LD *matrix, ULI cols, int id, struct globale*);

int load_files(int , int , struct globale*, unsigned*** , unsigned***, int***);
