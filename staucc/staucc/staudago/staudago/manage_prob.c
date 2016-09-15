#include "manage_prob.h"
#include <time.h>
//Reset the send/recv buffer
void cleanbuffers(int nMembranes, struct globale *G) {
  int k;
  
  for(k = 0; k<nMembranes*G->uliColumns[0]; k++) {
	G->ldM_send[k] = 0;
	//G->ldM_recv[k] = 0; //Inutile se sequenziale, usata solo in receive
  }
}



// Compute the highest order of reactions 
void get_HOR(int membid, struct globale *G) {
  ULI i, j, tmp=0, tmp2, current;

  for (i=0; i < G->uliColumns[membid]; i++) G->ldHOR[membid][i] = 0;
  for (i=0; i < G->uliRows[membid];    i++) G->uliOrder[membid][i] = 0;

  for (i = 0; i < G->uliRows[membid]; i++)
  	for (j = 0; j < G->uliColumns[membid]; j++)
  		G->uliOrder[membid][i] += G->unLeft_side[membid][i][j];				//18%% di 46%%

  for (j = 0; j < G->uliColumns[membid]; j++) {
	current = 0;
  	tmp2 = 0;

  	for (i = 0; i < G->uliRows[membid]; i++) {
		if (G->unLeft_side[membid][i][j] != 0) {
  	        	switch (G->uliOrder[membid][i]) {					//23%% di 46%%
  	        		case 1:
  	        			tmp  = 1;
  	        			tmp2 = 1;
  	        			break;
  	        		case 2:
  	        			if (G->unLeft_side[membid][i][j] == 1) {
  	        				tmp  = 2;
  	        				tmp2 = 2;
  	        			} else {
  	        				tmp  = 2.0 + (1.0 / (G->M[membid][j] - 1.0));
  	        				tmp2 = 3;
  	        			}
  	        			break;
  	        		case 3:
  	        			if (G->unLeft_side[membid][i][j] == 1) {
  	        				tmp  = 3;
  	        				tmp2 = 4;
  	        			} else if (G->unLeft_side[membid][i][j] == 2) {
  	        				tmp  = 1.5 * (2.0 + (1.0 / (G->M[membid][j] - 1.0)));
  	        				tmp2 = 5;
  	        			} else {
  	        				tmp  = 3.0 + (1.0 / (G->M[membid][j] - 1.0)) + (2.0 / (G->M[membid][j] - 2.0));
  	        				tmp2 = 6;
  	        			}
  	        			break;
  	        		default:
  	        			break;
  	        	}
  	        }
  	        if (tmp2 > current) {
  	        	current = tmp2;
  	        	G->ldHOR[membid][j] = tmp;
  	        }
  	}
  }
}



// Reinitialize the multisets
void feeding (int membid, struct globale *G) {
  int i;

  for (i = 0; i < G->uliColumns[membid]; i++)
  	  if (G->ldM_feed[membid][i] > 0) G->M[membid][i] = G->ldM_feed[membid][i];
}



// Identify critical reactions 
void rule_max_appl(LD **M, int membid, struct globale *G) {
  unsigned long c=0,d=0;
  LD tmp=0;
  ULI uliCheck;
  struct entry *curr_node=NULL;

  /* The number of possible application is equal to the lower
   * amount of molecules of the reactants
   */

  curr_node = G->LSlist[membid];

  for(c=0;c<G->uliRows[membid];c++) {
  	uliCheck = ULONG_MAX;
  	  
  	for(d=0; d<curr_node->iRows; d++) {
  		tmp = G->M[membid][curr_node->uliMatr[d][0]] / ( curr_node->uliMatr[d][1] * 1.0);
  		if(tmp < uliCheck) uliCheck = tmp;
  	}
  	if ((uliCheck < 10) && (uliCheck > 0)) G->uliCritical[membid][c] = 1;
  	else G->uliCritical[membid][c] = 0;
  	curr_node = curr_node->next;
  }
}



// Free space in the membrane 
LD getFS(int membid, struct globale *G) {
  int i;
  
  G->ldFS[membid] = G->ldMembSize[membid];

  for(i = 0; i<G->uliColumns[membid]; i++) G->ldFS[membid] -= (G->M[membid][i] * G->ldMolSize[membid][i]);
  
  return(G->ldFS[membid]);
}



// Compute the propensity functions: stochastic constant * species multiplicity 
void rule_prob(int membid, LD t, int uliStep, struct globale *G) {
  unsigned long c=0, d=0, e=0;
  struct entry *curr_node=NULL;

  getFS(membid,G);

  curr_node = G->LSlist[membid];
  memcpy(G->ldR_prob[membid], G->ldC_vector[membid], G->uliRows[membid] * sizeof(LD));

  for(c=0; c<G->uliRows[membid]; c++) {
  	//iTgtcheck		
	
  	//divide the stochastich constant for the available volume (only for reactions)
  	//divide the stochastic constant by the volume size
  	if (G->iTgt_vector[membid][c] == -1) {
  	        G->ldR_prob[membid][c] /= G->ldFS[membid];
  	        G->ldR_prob[membid][c] /= G->ldMembSize[membid];
  	}
  	        
  	for(d=0; d<curr_node->iRows; d++){
  	        switch(curr_node->uliMatr[d][1]){
  	        	case 1:
  	        		G->ldR_prob[membid][c] *= G->M[membid][curr_node->uliMatr[d][0]];
  	        		break;
  	        	case 2:
  	        		G->ldR_prob[membid][c] *= (G->M[membid][curr_node->uliMatr[d][0]] * (G->M[membid][curr_node->uliMatr[d][0]] - 1)) * 0.5;
  	        		break;
  	        	default:
  	        		for(e=1; e<=curr_node->uliMatr[d][1];e++)
  	        			G->ldR_prob[membid][c] *= (G->M[membid][curr_node->uliMatr[d][0]] - e + 1.0) / (1.0*e);
  	        		break;
  	        }
  	}
  	curr_node = curr_node->next;
  }


  //print the propensity rule in the log file
  if((uliStep % (ULI) G->every) == 0) {
  	fprintf(G->rprob_fPtr[membid],"%LG ",t);
  	for(c=0; c<G->uliRows[membid]-1; c++) fprintf(G->rprob_fPtr[membid],"%LG ",G->ldR_prob[membid][c]);	//58,1%%%
  	fprintf(G->rprob_fPtr[membid],"%LG\n",G->ldR_prob[membid][G->uliRows[membid]-1]);
  }
}



//Extract tau and mu for a SSA step
void dpp_tossing(LD a_0, gsl_rng *r, int membid, struct globale *G) {
  unsigned rule=0;
  LD rnd=0.0, rnd_t=0.0,alpha=0.0;
  
  alpha = G->ldR_prob[membid][rule];

  rnd = gsl_rng_uniform_pos (r) * a_0;
  rnd_t = gsl_rng_uniform_pos (r);

  while(rnd > alpha){
  	  rule++;
  	  alpha += G->ldR_prob[membid][rule];
  }

  G->uliR_app[membid][rule]++;
  G->uliK_rule[membid][rule] = 1;

  G->tau[membid] = (1.0/ a_0) * log(1.0/rnd_t);
}



// Dynamic step of the process 
int dpp_step1(gsl_rng *r, int membid, LD* a0c, LD* t2, LD t, int uliStep, struct globale *G) {
  ULI j;
  LD a_0=0.0, eps = 0.03;
  int stepkind;

  //Reinitialize the multisets (for open systems)
  feeding(membid,G);

  //Set to zero the vectors
  for(j=0;j<G->uliRows[membid];j++) { G->uliR_app[membid][j]=0; G->ldR_prob[membid][j]=0.0; }

  //Compute the rules number of applications
  rule_max_appl(G->M,membid,G);

  //Compute the propensity functions
  rule_prob(membid, t, uliStep,G);										//64%%%

  for(j=0;j<G->uliRows[membid];j++) a_0 += G->ldR_prob[membid][j];


  if (G->iFlag_SSA[membid] == 0) { //The first time at least, as set in main.c
  	if(a_0 > 0) {
		stepkind = 1;
		get_tau1(eps, r, a_0, membid, a0c, t2,G);
  	} else {
  	        /* If no rule can be applied, then set tau to infinite and update 
  	         * the system according to the smallest tau and the received data
  	         */
  	        stepkind = 2;
  	        G->tau[membid] = LONG_MAX;
  	        //get_min_tau();
	}
  }
  else {
	/*This branch is needed to complete a SSA step began during the previous iteration*/
  	/*Tau is equal to the old value*/
  	stepkind = 34;
  	G->tau[membid] = G->ldTau_ssa[membid];
  }
  
  return(stepkind);
}


void dpp_step2(gsl_rng *r, int membid, int stepkind, LD* a0c, LD* t2, struct globale *G) {	

  if(stepkind == 1) {
  	get_tau2(r, membid, a0c, t2,G);
			
        if ((G->iFlag_SSA[membid] == 0) || (G->ldTau_ssa[membid] == G->tau[membid])) {
        	G->iFlag_SSA[membid] = 0;
        	//Update the system state according to the selected rules
        	if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G);
		else dpp_update_tgt_matrix(r, membid,G);
        }
        else G->ldTau_ssa[membid] -= G->tau[membid];
  }
  else if(stepkind == 34) {
  	if (G->ldTau_ssa[membid] == G->tau[membid]) {
  	        stepkind = 3;
  	        G->iFlag_SSA[membid] = 0;
				
  	        //Complete the SSA step: execute a rule and receive data
		if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G);
		else dpp_update_tgt_matrix(r, membid,G);
  	}
  	else {
  	        stepkind = 4;
  	        G->ldTau_ssa[membid] -= G->tau[membid];
  	}
  }
}



void dpp_step3(gsl_rng *r, int membid, struct globale *G) {	
  int j;
  fprintf(G->log_fPtr[membid],"FS negative\n");
  
  for (j=0; j< G->uliRows[membid]; j++) G->uliK_rule[membid][j] = 0;

  G->tau[membid] /= 2.0;
  G->ldTau_ssa[membid] += G->tau[membid];
  //Restore data
  memcpy(G->M[membid], G->Mbkp[membid], G->uliColumns[membid] * sizeof(LD));
  G->ldFS[membid] = G->ldFSbkp[membid];

  if(G->iFlagStep[membid] == 2 || G->iFlagStep[membid] == 3)
	for (j = 0; j < G->uliRows[membid]; j++)
		if (G->uliCritical[membid][j] == 0) G->uliK_rule[membid][j] = gsl_ran_poisson(r, G->tau[membid] * G->ldR_prob[membid][j]);

  //dpp_update(size, r, membid,G);
  if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G);
  else dpp_update_tgt_matrix(r, membid,G);
}



LD dpp_step4(int membid, LD t, LD uliStep_max, int uliStep, struct globale *G) {	
  if (G->tau[membid] == LONG_MAX) G->tau[membid] = uliStep_max - t + 1;
  t +=  G->tau[membid];
  //uliStep++;

  // Add a vector to the output buffer
  if (uliStep % (ULI) G->every == 0) {
  	  add_vector_in_matrix(G->ldM_matrix[membid],G->M[membid],G->uliColumns[membid],t,membid,G);
  	  G->uliRow_counter[membid]++;
  }

  //If the buffer is full, then write it in the output file
  if (G->uliRow_counter[membid] == (ULI) G->buffer_lines) {
  	  write_matrix(G->m_fPtr[membid],G->ldM_matrix[membid],G->uliIndexCol,membid,G);
  	  G->uliRow_counter[membid] = 0;
  }

  //printf("%d: %LG %LG\n",membid,t,tau[membid]);  
  return(t);
}



//Update the state and send data according to communication rules - target vector
void dpp_update_tgt_vector(gsl_rng *r, int id, struct globale *G) {
  ULI i,j,k;
  int iTgt;
  LD* sendV;

  sendV = (LD*) malloc(G->uliColumns[id] * sizeof(LD));
  
  for(i=0; i<G->uliRows[id]; i++) {
	if(G->uliK_rule[id][i] > 0){
	
  		//printf("Eseguo %ld in %d (%d)\n",i,id,G->iTgt_vector[id][i]);
  	
  		//Execute an internal rule
  		if (G->iTgt_vector[id][i] == -1) {
			for(j=0;j<G->uliColumns[id];j++)
  		        	G->M[id][j] += (1.0*G->iVar[id][i][j]) * G->uliK_rule[id][i];
  		}
  		//Nondeterministic send
  		else if (G->iTgt_vector[id][i] == -2) {
  		        for(k=0;k<G->uliColumns[id];k++) {
  		        	//ldM_sendV[id][k] = unRight_side[id][i][k];
				sendV[k] = G->unRight_side[id][i][k];
  		        	G->M[id][k] -= G->uliK_rule[id][i] * G->unLeft_side[id][i][k];
  		        }
  		        for(j=0; j<G->uliK_rule[id][i]; j++) {
  		        	iTgt = gsl_rng_uniform_int(r, 0)+1; //era size-1
				for(k=0;k<G->uliColumns[id];k++)
					G->ldM_send[(iTgt*G->uliColumns[id])+k] += sendV[k];
  		        }
  		}
  		else {
  		        //Execute a deterministic communication rule
  		        for(k=0;k<G->uliColumns[id];k++) {
  		        	//ldM_sendV[id][k] = uliK_rule[id][i] * (1.0 *unRight_side[id][i][k]);
				G->ldM_send[(G->iTgt_vector[id][i]*G->uliColumns[id])+k] += G->uliK_rule[id][i] * (1.0 *G->unRight_side[id][i][k]);
  		        	G->M[id][k]    -= G->uliK_rule[id][i] * (1.0 *G->unLeft_side[id][i][k]);
  		        }
  		}
  	  }
  }
  free(sendV);
}



//Update the state and send data according to communication rules - target matrix
void dpp_update_tgt_matrix(gsl_rng *r, int id, struct globale *G) {
  ULI i,j,k;
  int *iMflag=NULL;

  iMflag = calloc(1, sizeof(int));//era size

  //Fill the buffer matrix
  for(i=0; i<G->uliRows[id]; i++) {
  	if(G->uliK_rule[id][i] > 0) {
		for(j=0; j<G->uliColumns[id]; j++) {
  			iMflag[G->iTgt_matrix[id][i][j]] = 1;
			//Se togli cast LD sbaglia moltiplicazione*/
			G->ldM_send[(G->iTgt_matrix[id][i][j]*G->uliColumns[id])+j] += 1.0 * ((LD) G->uliK_rule[id][i] * (LD)G->iVarSend[id][i][j]); 
			//se togli cast LD sbaglia moltiplicazione
  		}
	}
  }

  //Update the internal state of the process*/
  for(i=0;i<G->uliRows[id];i++)
  	for(k=0;k<G->uliColumns[id];k++)
  		//if(iTgt_matrix[id][i][k] != id) // IN OGNI CASO DEVE TOGLIERE LA PARTE SINISTRA
  		G->M[id][k] -= (LD)(G->uliK_rule[id][i] * G->unLeft_side[id][i][k]);

  free(iMflag);
}





void receive(int nMstart, int nMstop, int nMembranes, struct globale *G) {
  int k, kk, j;

  //MPI_Allreduce(G->ldM_send, G->ldM_recv, nMembranes*G->uliColumns[0], MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for(k=nMstart;k<=nMstop;k++) {
  	kk = k-nMstart;
 	for(j=0;j<G->uliColumns[kk];j++)
  		G->M[kk][j] +=  G->ldM_send[(kk*G->uliColumns[kk])+j];
		//G->M[kk][j] +=  G->ldM_recv[(kk*G->uliColumns[kk])+j];
  }
}



//tauleaping-like step
void get_tau1(LD eps, gsl_rng *r, LD a_0, int membid, LD* a0c, LD* t2, struct globale *G){
	LD  mu=0.0, sigma=0.0, max, rnd_t, tmp = LONG_MAX, a_0_c = 0.0, tau1 = LONG_MAX, tau2, tmp2, tmp3, tmp4; //, alpha
	ULI j,k;

	for (j=0; j< G->uliRows[membid]; j++) G->uliK_rule[membid][j] = 0;

	for (j = 0; j < G->uliColumns[membid]; j++) {
		max = 1.0;
		mu = 0.0;
		sigma = 0.0;
		
		/* For each rule, calculate the mean and standard deviation 
		 * of the changes in the propensity functions
		 */
		for (k = 0; k < G->uliRows[membid]; k++) {
			if (G->uliCritical[membid][k] == 0) {
				mu += (LD)(G->iVar[membid][k][j]) * G->ldR_prob[membid][k];
				sigma += pow((LD)(G->iVar[membid][k][j]),2) * G->ldR_prob[membid][k];
			}
		}

		//Compute the first candidate for the length of the time step
		tmp2 = ((eps * G->M[membid][j]) / G->ldHOR[membid][j]);
		if (tmp2 >= 1) max = tmp2;

		if (mu != 0 && sigma != 0) {
			tmp3 = max/fabs(mu);
			tmp4 = pow(max,2)/sigma;
			if (tmp3 < tmp4) tmp = tmp3;
			else tmp = tmp4;
		}
		if (tmp < tau1) tau1 = tmp;
	}

	//If the first candidate is too short, then execute a SSA step
	if (tau1 < 10.0 * (1.0/ a_0)) {
		G->iFlagStep[membid] = 1;
		/*Extract the rule to execute and get the tau value*/
		dpp_tossing(a_0,r,membid,G);
		G->ldTau_ssa[membid] = G->tau[membid];
		G->iFlag_SSA[membid] = 1;
		/*Get the smallest tau among the processes*/
	}
	else { //else 1
		//Compute the second candidate for the length of the time step
		for (k = 0; k < G->uliRows[membid]; k++) if (G->uliCritical[membid][k] == 1) a_0_c += G->ldR_prob[membid][k];

		rnd_t = gsl_rng_uniform_pos(r);
		tau2 = (1.0/ a_0_c) * log(1.0/rnd_t);

		//Extract the rules according to the chosen time increment
		if (tau1 < tau2) {
			G->iFlagStep[membid] = 2;
			G->tau[membid] = tau1;
		}
		else { //else 2
			G->iFlagStep[membid]=3;
			G->tau[membid] = tau2;

			a0c[membid] = a_0_c;
			t2[membid] = tau2;
		} //else 2
	} //else 1
}


void get_tau2(gsl_rng *r, int membid, LD* a_0_c, LD* tau2, struct globale *G){
	ULI i;
	LD rnd_t, alpha;
	
	if(G->iFlagStep[membid] == 2) {
		for (i = 0; i < G->uliRows[membid]; i++)
			if (G->uliCritical[membid][i] == 0)
				G->uliK_rule[membid][i] = gsl_ran_poisson(r, G->tau[membid] * G->ldR_prob[membid][i]);
	}
	else if(G->iFlagStep[membid] == 3){ //else 2
		for (i = 0; i < G->uliRows[membid]; i++) {
			if (G->uliCritical[membid][i] == 0) {
				G->uliK_rule[membid][i] = gsl_ran_poisson(r, G->tau[membid] * G->ldR_prob[membid][i]);
				G->ldR_prob_c[membid][i] = 0.0;
			}
			else G->ldR_prob_c[membid][i] = G->ldR_prob[membid][i];
		}
		i=0;
			
		// Select the critical reaction to execute in a SSA-like manner 
		if(tau2[membid] == G->tau[membid]) {
			rnd_t = gsl_rng_uniform_pos (r) * a_0_c[membid];
			alpha = G->ldR_prob_c[membid][i];
			while(rnd_t > alpha) {
				i++;
				alpha += G->ldR_prob_c[membid][i];
			}
			G->uliK_rule[membid][i] = 1;
		}
	} 
}



int checkNegFS(int nMstart, int nMstop, struct globale *G) {
  int k;
  LD minl=1000000.0, mint;
  
  for(k=nMstart;k<=nMstop;k++) {
	if(G->ldFS[k-nMstart] < minl) minl = G->ldFS[k-nMstart];
  }

  //MPI_Allreduce(&minl, &mint, 1, MPI_LONG_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  
  mint = minl;

  if (mint < 0) for(k=nMstart;k<=nMstop;k++) G->checkFS[k-nMstart] = 1;
  else		for(k=nMstart;k<=nMstop;k++) G->checkFS[k-nMstart] = 0;
  
  if(mint < 0) k = 1; else k = 0;
  
  return(k);
}



// Get the current smallest time increment value 
LD get_min_tau(LD mymin) {
  LD min;
  
  //MPI_Allreduce(&mymin, &min, 1, MPI_LONG_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  min = mymin;

  return(min);
}



