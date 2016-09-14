#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>

#ifndef _P_SIM_H
	#include "p_sim.h"
	#define _P_SIM_H
#endif
#ifndef _IO_FILE_H
	#include "io_file.h"
	#define _IO_FILE_H
#endif

void cleanbuffers(int, struct globale*);
void feeding(int, struct globale*);
void rule_max_appl(LD **, int, struct globale*);
void rule_prob(int, LD, int, struct globale*);
void dpp_tossing(LD, gsl_rng *, int, struct globale*);

int  dpp_step1(int, gsl_rng *, int, LD*, LD*, LD, int, struct globale*);
void dpp_step2(int, gsl_rng *, int, int, LD*, LD*, struct globale*);
void dpp_step3(int, gsl_rng *, int, struct globale*);
LD   dpp_step4(int, LD, LD, int, struct globale*);

void dpp_update_tgt_vector(int, gsl_rng *, int, struct globale*);
void dpp_update_tgt_matrix(int, gsl_rng *, int, struct globale*);

void receive(int, int, int, struct globale*);

void get_tau1(LD, gsl_rng *, LD, int, LD*, LD*, struct globale*);
void get_tau2(gsl_rng *, int, LD*, LD*, struct globale*);

void sum_prob(struct globale*);
LD get_min_tau(LD);
void get_HOR(int, struct globale*);
int checkNegFS(int, int, struct globale*);
LD getFS(int , struct globale*);

