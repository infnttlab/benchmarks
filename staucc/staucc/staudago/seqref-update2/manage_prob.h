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

int  dpp_step1(gsl_rng *, int, LD*, LD*, LD, int, struct globale*);
void dpp_step2(gsl_rng *, int, int, LD*, LD*, struct globale*, int);
void dpp_step3(gsl_rng *, int, struct globale*, int);
LD   dpp_step4(int, LD, LD, int, struct globale*);

void dpp_update_tgt_vector(gsl_rng *, int, struct globale*, int);
void dpp_update_tgt_matrix(gsl_rng *, int, struct globale*, int);

void receive(int, int, int, struct globale*);

void get_tau1(LD, gsl_rng *, LD, int, LD*, LD*, struct globale*, LD*);
void get_tau2(gsl_rng *, int, LD*, LD*, struct globale*);

void sum_prob(struct globale*);
void get_HOR(int, struct globale*, LD*);

int checkNegFS(int, int, struct globale*);
LD getFS(int , struct globale*);
LD get_min_tau(LD);
