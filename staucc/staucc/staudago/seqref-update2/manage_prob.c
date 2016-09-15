#include "manage_prob.h"

//Reset the send/recv buffer before the call to dpp_update (vector or matrix)
void cleanbuffers(int nMembranes, struct globale *G) {
  int k;
  
  for(k = 0; k<nMembranes*G->uliColumns[0]; k++) {
	G->ldM_send[k] = 0;
	//G->ldM_recv[k] = 0; //Inutile se sequenziale, usata solo in receive
  }
}



//Compute the highest order of reactions 
void get_HOR(int membid, struct globale *G, LD* ldHOR) {
  ULI i, j, tmp, tmp2;
  ULI* uliOrder, *current;
  unsigned ls;
  
  uliOrder = (ULI*) malloc(sizeof(ULI) * G->uliRows[membid]);
  current =  (ULI*) malloc(sizeof(ULI) * G->uliColumns[membid]);
  
  for(i=0; i < G->uliColumns[membid]; i++) { ldHOR[i] = 0; current[i] = 0; }
  for(i=0; i < G->uliRows[membid];    i++) uliOrder[i] = 0;

  for(i = 0; i < G->uliRows[membid]; i++)
  	for(j = 0; j < G->uliColumns[membid]; j++)
  		//G->uliOrder[membid][i] += G->unLeft_side[membid][i][j];				
		uliOrder[i] += G->unLeft_side[G->side_index[membid]+(i*G->uliColumns[membid])+j];
/*
  //questo for va GIRATO
  for(j = 0; j < G->uliColumns[membid]; j++) {
	current = 0; tmp2 = 0;
  	for(i = 0; i < G->uliRows[membid]; i++) {
		indaux = G->side_index[membid]+(i*G->uliColumns[membid])+j;
		
		if (G->unLeft_side[indaux] != 0) {
  	        	switch (uliOrder[i]) {					
  	        		case 1:
  	        			tmp  = 1;
  	        			tmp2 = 1;
  	        			break;
  	        		case 2:
					if (G->unLeft_side[indaux] == 1) {
  	        				tmp  = 2;
  	        				tmp2 = 2;
  	        			} else {
  	        				tmp  = 2.0 + (1.0 / (G->M[membid][j] - 1.0));
  	        				tmp2 = 3;
  	        			}
  	        			break;
  	        		case 3:
					if (G->unLeft_side[indaux] == 1) {
  	        				tmp  = 3;
  	        				tmp2 = 4;
					} else if (G->unLeft_side[indaux] == 2) {
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
  	        	ldHOR[j] = tmp;
  	        }
  	}
  }
*/

  for(i = 0; i < G->uliRows[membid]; i++) {
	for(j = 0; j < G->uliColumns[membid]; j++) {
		tmp  = 0;
		tmp2 = 0;
		ls = G->unLeft_side[G->side_index[membid]+(i*G->uliColumns[membid])+j];
		
		if (ls != 0) {
  	        	switch (uliOrder[i]) {					
  	        		case 1:
  	        			tmp  = 1;
  	        			tmp2 = 1;
  	        			break;
  	        		case 2:
					if (ls == 1) {
  	        				tmp  = 2;
  	        				tmp2 = 2;
  	        			} else {
  	        				tmp  = 2.0 + (1.0 / (G->M[membid][j] - 1.0));
  	        				tmp2 = 3;
  	        			}
  	        			break;
  	        		case 3:
					if (ls == 1) {
  	        				tmp  = 3;
  	        				tmp2 = 4;
					} else if (ls == 2) {
  	        				tmp  = 1.5 * (2.0 + (1.0 / (G->M[membid][j] - 1.0)));
  	        				tmp2 = 5;
  	        			} else {
  	        				tmp  = 3.0 + (1.0 / (G->M[membid][j] - 1.0)) + (2.0 / (G->M[membid][j] - 2.0));
  	        				tmp2 = 6;
  	        			}
  	        			break;
  	        		default: //qui non dovrei arrivare
  	        			break;
  	        	}
  	        }
  	        if (tmp2 > current[j]) {
  	        	current[j] = tmp2;
  	        	ldHOR[j]   = tmp;
  	        }
  	}
  }
  
  free(uliOrder);
  free(current);
}



// Reinitialize the multisets
void feeding (int membid, struct globale *G) {
  int i;

  for(i = 0; i < G->uliColumns[membid]; i++)
  	  if (G->ldM_feed[membid][i] > 0) G->M[membid][i] = G->ldM_feed[membid][i];
}



// Identify critical reactions 
void rule_max_appl(LD **M, int membid, struct globale *G) {
  unsigned long c=0,d=0;
  LD tmp=0;
  ULI uliCheck;

  for(c=0;c<G->uliRows[membid];c++) {
  	uliCheck = ULONG_MAX;
	
	for(d=0; d<G->uliColumns[membid]; d++) {
		//if(G->unLeft_side[membid][c][d]!=0)
		if(G->unLeft_side[G->side_index[membid]+(c*G->uliColumns[membid])+d] != 0)
			//quante se ne possono fare di questo tipo = quante presenti /  quante richieste nella reazione
			//tmp = G->M[membid][d] / ( G->unLeft_side[membid][c][d] * 1.0); 
			tmp = G->M[membid][d] / ( G->unLeft_side[G->side_index[membid]+(c*G->uliColumns[membid])+d] * 1.0);
  		if(tmp < uliCheck) uliCheck = tmp;
	}
	
  	if ((uliCheck < 10) && (uliCheck > 0)) G->uliCritical[membid][c] = 1;
  	else G->uliCritical[membid][c] = 0;
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

  getFS(membid,G);

  //curr_node = G->LSlist[membid];
  memcpy(G->ldR_prob[membid], G->ldC_vector[membid], G->uliRows[membid] * sizeof(LD));

  for(c=0; c<G->uliRows[membid]; c++) {
  	//iTgtcheck		
	
  	//divide the stochastich constant for the available volume (only for reactions)
  	//divide the stochastic constant by the volume size
  	if (G->iTgt_vector[membid][c] == -1) {
  	        G->ldR_prob[membid][c] /= G->ldFS[membid];
  	        G->ldR_prob[membid][c] /= G->ldMembSize[membid];
  	}

	for(d=0; d<G->uliColumns[membid]; d++) {
			//switch(G->unLeft_side[membid][c][d]){
			switch(G->unLeft_side[G->side_index[membid]+(c*G->uliColumns[membid])+d]){
			case 0: 
				break;
  	        	case 1:
  	        		G->ldR_prob[membid][c] *= G->M[membid][d];
  	        		break;
  	        	case 2:
  	        		G->ldR_prob[membid][c] *= (G->M[membid][d] * (G->M[membid][d] - 1)) * 0.5;
  	        		break;
  	        	default:
  	        		//for(e=1; e<=G->unLeft_side[membid][c][d];e++)
				for(e=1; e<=G->unLeft_side[G->side_index[membid]+(c*G->uliColumns[membid])+d];e++)
  	        			G->ldR_prob[membid][c] *= (G->M[membid][d] - e + 1.0) / (1.0*e);
  	        		break;
  	        }	
	}
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

  //G->uliR_app[membid][rule]++;
  G->uliK_rule[membid][rule] = 1;

  G->tau[membid] = (1.0/ a_0) * log(1.0/rnd_t);
}



// Dynamic step of the process 
int dpp_step1(gsl_rng *r, int membid, LD* a0c, LD* t2, LD t, int uliStep, struct globale *G) {
  ULI j;
  LD a_0=0.0, eps = 0.03;
  int stepkind;
  LD* ldHOR;

  //Reinitialize the multisets (for open systems)
  feeding(membid,G);

  //Set to zero the vectors
  //for(j=0;j<G->uliRows[membid];j++) { G->uliR_app[membid][j]=0; G->ldR_prob[membid][j]=0.0; }
  for(j=0;j<G->uliRows[membid];j++) G->ldR_prob[membid][j]=0.0; 
  
  //Compute the rules number of applications
  rule_max_appl(G->M,membid,G);

  //Compute the propensity functions
  rule_prob(membid, t, uliStep,G);										//64%%%

  for(j=0;j<G->uliRows[membid];j++) a_0 += G->ldR_prob[membid][j];


  if (G->iFlag_SSA[membid] == 0) { //The first time at least, as set in main.c
  	if(a_0 > 0) {
		stepkind = 1;
		ldHOR = (LD *) malloc(sizeof(LD) * G->uliColumns[membid]);	//Highest order of reaction vector
		
		get_HOR(membid,G,ldHOR);
		get_tau1(eps, r, a_0, membid, a0c, t2,G,ldHOR);
		
		free(ldHOR);
  	} else {
  	        /* If no rule can be applied, then set tau to infinite and update 
  	         * the system according to the smallest tau and the received data
  	         */
  	        stepkind = 2;
  	        G->tau[membid] = LONG_MAX;
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




void dpp_step2(gsl_rng *r, int membid, int stepkind, LD* a0c, LD* t2, struct globale *G, int nMembranes) {	

  if(stepkind == 1) {
  	get_tau2(r, membid, a0c, t2,G);
			
        if ((G->iFlag_SSA[membid] == 0) || (G->ldTau_ssa[membid] == G->tau[membid])) {
        	G->iFlag_SSA[membid] = 0;
        	//Update the system state according to the selected rules
        	if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G,nMembranes);
		else 		      dpp_update_tgt_matrix(r, membid,G,nMembranes);
        }
        else G->ldTau_ssa[membid] -= G->tau[membid];
  }
  else if(stepkind == 34) {
  	if (G->ldTau_ssa[membid] == G->tau[membid]) {
  	        stepkind = 3;
  	        G->iFlag_SSA[membid] = 0;
				
  	        //Complete the SSA step: execute a rule and receive data
		if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G,nMembranes);
		else 		      dpp_update_tgt_matrix(r, membid,G,nMembranes);
  	}
  	else {
  	        stepkind = 4;
  	        G->ldTau_ssa[membid] -= G->tau[membid];
  	}
  }
}



void dpp_step3(gsl_rng *r, int membid, struct globale *G, int nMembranes) {	
  int j;
  fprintf(G->log_fPtr[membid],"FS negative\n");
  
  for(j=0; j< G->uliRows[membid]; j++) G->uliK_rule[membid][j] = 0;

  G->tau[membid] /= 2.0;
  G->ldTau_ssa[membid] += G->tau[membid];
  //Restore data
  memcpy(G->M[membid], G->Mbkp[membid], G->uliColumns[membid] * sizeof(LD));
  G->ldFS[membid] = G->ldFSbkp[membid];

  if(G->iFlagStep[membid] == 2 || G->iFlagStep[membid] == 3)
	for(j = 0; j < G->uliRows[membid]; j++)
		if (G->uliCritical[membid][j] == 0) G->uliK_rule[membid][j] = gsl_ran_poisson(r, G->tau[membid] * G->ldR_prob[membid][j]);

  //dpp_update(size, r, membid,G);
  if(G->iTgtcheck == 1) dpp_update_tgt_vector(r, membid,G, nMembranes);
  else 			dpp_update_tgt_matrix(r, membid,G, nMembranes);
}



LD dpp_step4(int membid, LD t, LD uliStep_max, int uliStep, struct globale *G) {	
  if (G->tau[membid] == LONG_MAX) G->tau[membid] = uliStep_max - t + 1;
  t +=  G->tau[membid];

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

  return(t);
}



//Update the state and send data according to communication rules - target vector
void dpp_update_tgt_vector(gsl_rng *r, int id, struct globale *G, int nMembranes) {
  ULI i,j,k;
  int iTgt;
  LD* sendV;
  int indaux;

  sendV = (LD*) malloc(G->uliColumns[id] * sizeof(LD));
  
  for(i=0; i<G->uliRows[id]; i++) {
	if(G->uliK_rule[id][i] > 0){
  		//Execute an internal rule
  		if (G->iTgt_vector[id][i] == -1) {
			for(k=0;k<G->uliColumns[id];k++)
  		        	//G->M[id][j] += (1.0*G->iVar[id][i][j]) * G->uliK_rule[id][i];
				G->M[id][k] += (1.0*G->iVar[G->side_index[id]+(i*G->uliColumns[id])+k]) * G->uliK_rule[id][i];
  		}
  		//Nondeterministic send
  		else if (G->iTgt_vector[id][i] == -2) {
  		        for(k=0;k<G->uliColumns[id];k++) {
				indaux = G->side_index[id]+(i*G->uliColumns[id])+k;
				
  		        	//ldM_sendV[id][k] = unRight_side[id][i][k];
				//sendV[k] = G->unRight_side[id][i][k];
				sendV[k] = G->unRight_side[indaux];
				
  		        	//G->M[id][k] -= G->uliK_rule[id][i] * G->unLeft_side[id][i][k];
				G->M[id][k] -= G->uliK_rule[id][i] * G->unLeft_side[indaux];
  		        }
  		        for(j=0; j<G->uliK_rule[id][i]; j++) {
  		        	iTgt = gsl_rng_uniform_int(r, nMembranes-1)+1; //era size-1
				
				for(k=0;k<G->uliColumns[id];k++)
					G->ldM_send[(iTgt*G->uliColumns[id])+k] += sendV[k];
  		        }
  		}
  		else {
  		        //Execute a deterministic communication rule
  		        for(k=0;k<G->uliColumns[id];k++) {
				indaux = G->side_index[id]+(i*G->uliColumns[id])+k;
				
  		        	//ldM_sendV[id][k] = uliK_rule[id][i] * (1.0 *unRight_side[id][i][k]);
				//G->ldM_send[(G->iTgt_vector[id][i]*G->uliColumns[id])+k] += G->uliK_rule[id][i] * (1.0 *G->unRight_side[id][i][k]);
				G->ldM_send[(G->iTgt_vector[id][i]*G->uliColumns[id])+k] += G->uliK_rule[id][i] * (1.0*G->unRight_side[indaux]);
  		        	//G->M[id][k]    -= G->uliK_rule[id][i] * (1.0 *G->unLeft_side[id][i][k]);
				G->M[id][k] -= G->uliK_rule[id][i] * (1.0 *G->unLeft_side[indaux]);
  		        }
  		}
  	  }
  }
  free(sendV);
}



//Update the state and send data according to communication rules - target matrix
void dpp_update_tgt_matrix(gsl_rng *r, int id, struct globale *G, int nMembranes) {
  ULI i,j;
  int *iMflag=NULL;
  int indaux;
  
  iMflag = malloc(nMembranes*sizeof(int));//era size
  for(i=0;i<nMembranes; i++) iMflag[i] = 0;

  //Fill the buffer matrix
  for(i=0; i<G->uliRows[id]; i++) {
  	if(G->uliK_rule[id][i] > 0) {
		for(j=0; j<G->uliColumns[id]; j++) {
			indaux = G->side_index[id]+(i*G->uliColumns[id])+j;
			
  			iMflag[G->iTgt_matrix[indaux]] = 1;
			//G->ldM_send[(G->iTgt_matrix[id][i][j]*G->uliColumns[id])+j] += 1.0 * ((LD) G->uliK_rule[id][i] * (LD)G->iVarSend[id][i][j]); 
			G->ldM_send[(G->iTgt_matrix[indaux]*G->uliColumns[id])+j] += 1.0 * ((LD) G->uliK_rule[id][i] * (LD)G->iVarSend[indaux]);
			//se togli cast LD sbaglia moltiplicazione
  		}
	}
  }

  //Update the internal state of the process*/
  for(i=0;i<G->uliRows[id];i++)
  	for(j=0;j<G->uliColumns[id];j++)
  		//G->M[id][k] -= (LD)(G->uliK_rule[id][i] * G->unLeft_side[id][i][k]);
		G->M[id][j] -= (LD)(G->uliK_rule[id][i] * G->unLeft_side[G->side_index[id]+(i*G->uliColumns[id])+j]);

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
void get_tau1(LD eps, gsl_rng *r, LD a_0, int membid, LD* a0c, LD* t2, struct globale *G, LD* ldHOR){
  //LD  mu=0.0, sigma=0.0
  LD *mu, *sigma;
  LD max, rnd_t, tmp = LONG_MAX, a_0_c = 0.0, tau1 = LONG_MAX, tau2, tmp2, tmp3, tmp4; //, alpha
  ULI j,k;
  int ivaux;

  mu    = (LD*) malloc(G->uliColumns[membid] * sizeof(LD));
  sigma = (LD*) malloc(G->uliColumns[membid] * sizeof(LD));
  
  for(j = 0; j < G->uliColumns[membid]; j++) { mu[j] = 0.0; sigma[j] = 0.0; }
  for(j=0; j< G->uliRows[membid]; j++) G->uliK_rule[membid][j] = 0;
/*
  //questo for va girato
  for(j = 0; j < G->uliColumns[membid]; j++) {
  	  max = 1.0;
  	  mu = 0.0;
  	  sigma = 0.0;
  	  
  	  // For each rule, calculate the mean and standard deviation 
  	  // of the changes in the propensity functions
  	  for(k = 0; k < G->uliRows[membid]; k++) {
  		  if (G->uliCritical[membid][k] == 0) {
			  ivaux = G->iVar[G->side_index[membid]+(k*G->uliColumns[membid])+j];
  			  mu    +=     (LD)(ivaux)    * G->ldR_prob[membid][k];
  			  sigma += pow((LD)(ivaux),2) * G->ldR_prob[membid][k];
  		  }
  	  }

  	  //Compute the first candidate for the length of the time step
  	  tmp2 = ((eps * G->M[membid][j]) / ldHOR[j]);
  	  if (tmp2 >= 1) max = tmp2;

  	  if (mu != 0 && sigma != 0) {
  		  tmp3 = max/fabs(mu);
  		  tmp4 = pow(max,2)/sigma;
  		  if (tmp3 < tmp4) tmp = tmp3;
  		  else tmp = tmp4;
  	  }
  	  if (tmp < tau1) tau1 = tmp;
  }
*/
  // For each rule, calculate the mean and standard deviation 
  // of the changes in the propensity functions
  for(k = 0; k < G->uliRows[membid]; k++) {
  	if (G->uliCritical[membid][k] == 0) {
		for(j = 0; j < G->uliColumns[membid]; j++) {
			ivaux = G->iVar[G->side_index[membid]+(k*G->uliColumns[membid])+j];
  			mu[j]    +=     (LD)(ivaux)    * G->ldR_prob[membid][k];
  			sigma[j] += pow((LD)(ivaux),2) * G->ldR_prob[membid][k];
		}
  	}
  }

  for(j = 0; j < G->uliColumns[membid]; j++) {
  	max = 1.0;
  	//Compute the first candidate for the length of the time step
  	tmp2 = ((eps * G->M[membid][j]) / ldHOR[j]);
  	if (tmp2 >= 1) max = tmp2;

  	if (mu[j] != 0 && sigma[j] != 0) {
		tmp3 = max/fabs(mu[j]);
  		tmp4 = pow(max,2)/sigma[j];
  		if (tmp3 < tmp4) tmp = tmp3; else tmp = tmp4;
  	}
  	if (tmp < tau1) tau1 = tmp;
  }

  free(mu);
  free(sigma);

  //If the first candidate is too short, then execute a SSA step
  if (tau1 < 10.0 * (1.0/ a_0)) {
  	G->iFlagStep[membid] = 1;
  	//Extract the rule to execute and get the tau value
  	dpp_tossing(a_0,r,membid,G);
  	G->ldTau_ssa[membid] = G->tau[membid];
  	G->iFlag_SSA[membid] = 1;
  	/*Get the smallest tau among the processes*/
  }
  else { //else 1
       //Compute the second candidate for the length of the time step
       for(k = 0; k < G->uliRows[membid]; k++) if (G->uliCritical[membid][k] == 1) a_0_c += G->ldR_prob[membid][k];

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
  	  for(i = 0; i < G->uliRows[membid]; i++)
  		  if (G->uliCritical[membid][i] == 0)
  			  G->uliK_rule[membid][i] = gsl_ran_poisson(r, G->tau[membid] * G->ldR_prob[membid][i]);
  }
  else if(G->iFlagStep[membid] == 3){ //else 2
  	  for(i = 0; i < G->uliRows[membid]; i++) {
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
/*
  if (mint < 0) for(k=nMstart;k<=nMstop;k++) G->checkFS[k-nMstart] = 1;
  else		for(k=nMstart;k<=nMstop;k++) G->checkFS[k-nMstart] = 0;
*/  
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



