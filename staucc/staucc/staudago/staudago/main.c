//Occhio, potrebbe servire ulimit -s unlimited  (di default 10 M)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>


#include <sys/resource.h>
#include <sys/time.h> 

#include "io_file.h"
#include "manage_prob.h"


void initialize(LD uliStep_max, gsl_rng *r, struct globale *G) {	
  ULI i=0, j=0;
  LD rnd=0, mymintau, globmintau;
  char file[20],app[20];
  int count, iRows_index=0;
  struct entry *curr_node=NULL;
  FILE* fPtr;
  int nMembranes, nMstart, nMstop, nMnum, k, kk;
  int* assignments;
  struct timeval tstart, tstop;
  float elapsed, ctime, rtime;
  int* stepkind;
  LD *a0c, *t2;
  int cn;
  LD t1, t3;
  LD t;
  int uliStep=0;
  int ee;
  int buffer_size =0;

  /** NOTE: Size is the number of the parallel processes, NOT 
   * the number of membranes
   */

  gsl_rng_set(r, time(NULL));
  
  rnd = gsl_rng_uniform_int (r, G->uliMAXSEED);


  if ((fPtr = fopen ("input/numMembranes.txt","r")) == NULL)	  {
  	  printf ("\nCannot obtain the numbers of membranes, the file input/numMembranes.txt is required\n");
  	  exit(0);
  }
  else {
  	  fscanf(fPtr,"%d",&nMembranes);
  	  fclose(fPtr);
  }

  assignments = (int *) malloc(nMembranes*sizeof(int));

  
  nMnum = nMembranes;
  nMstart = 0;
  nMstop = nMembranes-1;

  G->unLeft_side	= (unsigned***) malloc(nMnum * sizeof(void*));
  G->unRight_side    	= (unsigned***) malloc(nMnum * sizeof(void*));
  G->ldC_vector		= (LD **) malloc(nMnum * sizeof(void*));
  G->M			= (LD **) malloc(nMnum * sizeof(void*));
  G->Mbkp  		= (LD **) malloc(nMnum * sizeof(void*));
  G->ldM_feed		= (LD **) malloc(nMnum * sizeof(void*));
  G->ldMolSize		= (LD **) malloc(nMnum * sizeof(void*));
  G->ldMembSize		= (LD *)  malloc(nMnum * sizeof(LD));
  G->iTgt_vector	= (int **)  malloc(nMnum * sizeof(void*));
  G->iTgt_matrix	= (int ***) malloc(nMnum * sizeof(void*));
  G->indexes		= (int **) malloc(nMnum * sizeof(void*));
  G->LSlist		= (struct entry **) malloc(nMnum *sizeof(void*));
  G->RSlist		= (struct entry **) malloc(nMnum *sizeof(void*));

  G->uliR_app		= (ULI**) malloc(nMnum * sizeof(void*));
  G->uliK_rule		= (ULI**) malloc(nMnum * sizeof(void*)); 
  G->uliOrder		= (ULI**) malloc(nMnum * sizeof(void*)); //resettato in get_HOR
  G->uliCritical	= (ULI**) malloc(nMnum * sizeof(void*));
  G->uliRow_counter	= (ULI*) malloc(nMnum * sizeof(ULI));
  G->ldR_prob		= (LD **) malloc(nMnum * sizeof(void*));  
  G->ldHOR 		= (LD **) malloc(nMnum * sizeof(void*));  //resettato in get_HOR
  G->ldR_prob_c		= (LD **) malloc(nMnum * sizeof(void*));  

  G->iVar  		= (int ***) malloc(nMnum * sizeof(void*));	
  G->iVarSend		= (int ***) malloc(nMnum * sizeof(void*));
  
  G->ldM_matrix		= (LD **) malloc(nMnum * sizeof(void*));

  G->ldFS  		= (LD *) malloc(nMnum * sizeof(LD));
  G->ldFSbkp		= (LD *) malloc(nMnum * sizeof(LD));
  G->tau		= (LD *) malloc(nMnum * sizeof(LD));

  G->m_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  G->log_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  //G->fs_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  G->rprob_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  
  G->uliRows		= (ULI*) malloc(nMnum * sizeof(ULI));
  G->uliColumns		= (ULI*) malloc(nMnum * sizeof(ULI));

  stepkind		= (int*) malloc(nMnum * sizeof(int));
  a0c  			= (LD *) malloc(nMnum * sizeof(LD));
  t2  			= (LD *) malloc(nMnum * sizeof(LD));
  G->iFlag_SSA		= (int*) malloc(nMnum * sizeof(int)); 		
  G->iFlagStep		= (int*) malloc(nMnum * sizeof(int));
  G->ldTau_ssa		= (LD *) malloc(nMnum * sizeof(LD));
  G->checkFS		= (int*) malloc(nMnum * sizeof(int));
  
  for(i=0;i<nMnum;i++) { G->tau[i] = 0.0; G->iFlag_SSA[i] = 0; G->ldTau_ssa[i] = 0.0; G->uliRow_counter[i] = 0;}

  gettimeofday(&tstart,NULL);

  t = 0.0;
  
  G->iotime = 0.0;
  
  for(k=nMstart;k<=nMstop;k++) {  
  	kk = k-nMstart;
  	
  	/************
  	 ************
  	 OUTPUT FILES
  	 ************
  	 ************/
  	
  	//Create a log file name labelled with the process id
  	strcpy(file,"output/log_");
  	itoa(k,app);
  	strcat(file,app);
  	
  	//Create a new log file into the 'output' folder
  	if ((G->log_fPtr[kk] = fopen (file,"w")) == NULL) {
  	        printf ("\nCannot write in the current directory - Log file\nAbort\n");
  	        exit(0);
  	}
  	
  	//Write process information into the log file
  	fprintf(G->log_fPtr[kk],"uliStep_max %LG\n",uliStep_max);
	
  	//Create an output file name labelled with the process id
  	strcpy(file,"output/multi_");
  	strcat(file,app);
  	if ((G->m_fPtr[kk] = fopen (file,"w")) == NULL) {
  	        printf ("\nCannot write in the current directory - Output file\nAbort\n");
  	        exit(0);
  	}
  	
	/*create the file to store the free space for every iteration
  	strcpy(file,"output/fs_");
  	strcat(file,app);
  	if ((G->fs_fPtr[kk] = fopen (file,"w")) == NULL) {
  	        printf ("\nCannot write in the current directory - fs file\nAbort\n");
  	        exit(0);
  	}
  	*/
	
  	//create the file to store the propensity function values
  	strcpy(file,"output/rprob_");
  	strcat(file,app);
  	if ((G->rprob_fPtr[kk] = fopen (file,"w")) == NULL){
  	        printf ("\nCannot write in the current directory - Log file\nAbort\n");
  	        exit(0);
  	}
	
  	  
  	/***********
  	 ***********
  	 INPUT FILES
  	 ***********
  	 ***********/
  	
  	/* Read the input files */
  	G->iTgtcheck = load_files(k,nMstart,G);
	
	///Debug
	//printf("%d: %lu %lu\n",k, G->uliRows[kk], G->uliColumns[kk]);
	
	//uliColumns e' uguale per tutte, sono le specie, ma le righe cambiano


  	//Compile the lists with the reactions
  	G->LSlist[kk] = (struct entry *)malloc(sizeof(struct entry));
  	curr_node = G->LSlist[kk];
  	
  	for(i= 0; i<G->uliRows[kk]-1; i++) {
  	        curr_node->next = (struct entry *)malloc(sizeof(struct entry));
  	        curr_node = curr_node->next;
  	}
  	curr_node->next = NULL;
	
  	G->RSlist[kk] = (struct entry *)malloc(sizeof(struct entry));
  	curr_node = G->RSlist[kk];
  	
  	for(i= 0; i<G->uliRows[kk]-1; i++){
  	        curr_node->next = (struct entry *)malloc(sizeof(struct entry));
  	        curr_node = curr_node->next;
  	}
	curr_node->next = NULL;
  	        	
  	/* To optimize the space and speedup the computation only reactants are stored in these data structure */
  	
  	/* NOTA: non sarebe meglio buttare unLeft e unRght ed usare solo queste?
  	 * La ista non si puo' usare se andiamo su CUDA e poi effettivamente velocizzo:
  	 * Una reazione max ha 3 reactant e 3 risulatati, al min 1 e 1. Verificare e poi decidere
  	 */

  	// Left-hand side 
  	curr_node = G->LSlist[kk];	
  	for(i=0; i<G->uliRows[kk]; i++) {
  	        count = 0;
  	        for(j=0; j<G->uliColumns[kk]; j++) if(G->unLeft_side[kk][i][j]!=0) count++;
  	        
  	        curr_node->iRows = count;
  	        curr_node->uliMatr = malloc(sizeof(ULI) * count);
  	        for(j=0; j<count; j++) curr_node->uliMatr[j] = malloc(sizeof(ULI) * 2);
  	        
  	        iRows_index = 0;
  	        for(j=0; j<G->uliColumns[kk]; j++) {
  	        	if(G->unLeft_side[kk][i][j] != 0){
  	        		curr_node->uliMatr[iRows_index][0] = j;
  	        		curr_node->uliMatr[iRows_index][1] = G->unLeft_side[kk][i][j];
  	        		iRows_index++;
  	        	}
  	        }
  	        curr_node = curr_node->next;
  	}
  	
  	/* Right-hand side */
  	curr_node = G->RSlist[kk];
  	for(i=0; i<G->uliRows[kk]; i++){
  	        count = 0;
  	        for(j=0; j<G->uliColumns[kk]; j++) if(G->unRight_side[kk][i][j]!=0) count++;
  	        	
  	        curr_node->iRows = count;
  	        curr_node->uliMatr = malloc(sizeof(ULI) * count);
  	        for(j=0; j<count; j++) curr_node->uliMatr[j] = malloc(sizeof(ULI) * 2);
  	        
  	        iRows_index = 0;
  	        for(j=0; j<G->uliColumns[kk]; j++) {
  	        	if(G->unRight_side[kk][i][j]!=0){
  	        		curr_node->uliMatr[iRows_index][0] = j;
  	        		curr_node->uliMatr[iRows_index][1] = G->unRight_side[kk][i][j];
  	        		iRows_index++;
  	        	}
  	        }
  	        curr_node = curr_node->next;
  	}
  	
  	//Print the lists into the log file
  	curr_node = G->LSlist[kk];
  	fprintf(G->log_fPtr[kk], "\nLEFT SIDE\n");
  	for (i=0; i < G->uliRows[kk]; i++) {
  	        fprintf(G->log_fPtr[kk], "\niRows: %d\n", curr_node->iRows);
  	        for(j=0; j<curr_node->iRows; j++) fprintf(G->log_fPtr[kk], "%lu\t%lu\n", curr_node->uliMatr[j][0],curr_node->uliMatr[j][1]);
  	        curr_node = curr_node->next;
  	}
  	
  	fprintf(G->log_fPtr[kk], "\nRIGHT SIDE\n");
  	curr_node = G->RSlist[kk];
  	for (i=0; i < G->uliRows[kk]; i++) {
  	        fprintf(G->log_fPtr[kk], "\niRows: %d\n", curr_node->iRows);
  	        for(j=0; j<curr_node->iRows; j++) fprintf(G->log_fPtr[kk], "%lu\t%lu\n", curr_node->uliMatr[j][0],curr_node->uliMatr[j][1]);
  	        curr_node = curr_node->next;
  	}	
  	
  	G->uliR_app[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Application vector*/
  	G->uliK_rule[kk]   = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Tau leaping application vector*/
  	G->uliOrder[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Order of reactions vector*/
  	G->uliCritical[kk] = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Critical reactions vector*/
  	G->ldR_prob[kk]    = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	/* Propensity functions vector*/
  	G->ldHOR[kk]	   = (LD *) malloc(sizeof(LD) * G->uliColumns[kk]);	/* Highest order of reaction vector*/
  	G->ldR_prob_c[kk]  = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	/* Critical reactions propensity functions vector*/
  	G->iVar[kk]	   = (int **) malloc(G->uliRows[kk] * sizeof(void *) ); /* Variations matrix*/
  	G->iVarSend[kk]    = (int **) malloc(G->uliRows[kk] * sizeof(void *) );
  	
  	for(i=0;i<G->uliRows[kk];i++){
  	        G->iVar[kk][i]  	= malloc(G->uliColumns[kk]* sizeof(int) );
  	        G->iVarSend[kk][i]	= malloc(G->uliColumns[kk]* sizeof(int) );
  	
	        for(j=0; j<G->uliColumns[kk];j++) {
	        	G->iVar[kk][i][j]     = 0;
	        	G->iVarSend[kk][i][j] = 0;
	        }
	}
  	// Fill the variation matrix 
  	if(G->iTgtcheck == 1){
  	        for(i=0; i<G->uliRows[kk]; i++){
  	        	for(j=0; j<G->uliColumns[kk]; j++){
  	        		if(G->iTgt_vector[kk][i] == -1){
  	        			G->iVar[kk][i][j] -= G->unLeft_side[kk][i][j];
  	        			G->iVar[kk][i][j] += G->unRight_side[kk][i][j];
  	        		} //IF I have to send the output to another process...
  	        		else G->iVar[kk][i][j] -= G->unLeft_side[kk][i][j];
  	        	}
  	        }
  	}
  	else {
  	        //dpp_update = dpp_update_tgt_matrix;
  	        for(i=0; i<G->uliRows[kk]; i++){
  	        	for(j=0; j<G->uliColumns[kk]; j++){
  	        		if(G->iTgt_matrix[kk][i][j] == 0){ //era id
  	        			G->iVar[kk][i][j] -= G->unLeft_side[kk][i][j];
  	        			G->iVar[kk][i][j] += G->unRight_side[kk][i][j];
  	        			G->iVarSend[kk][i][j] = G->iVar[kk][i][j];
  	        		}
  	        		else{
  	        			G->iVar[kk][i][j]    -= G->unLeft_side[kk][i][j];
  	        			G->iVarSend[kk][i][j] = G->unRight_side[kk][i][j];
  	        		}
  	        	}
  	        }
  	}
  	//E da qui iVar e iVarSend restano costanti
	//iVarSend solo se tgt_matrix e usato in dpp_update_tgt_matrix
	
	
  	//Print the variations matrix into the log file
	fprintf(G->log_fPtr[kk], "\nVariations Matrix\n");
  	for(i=0;i<G->uliRows[kk];i++){
  	        for(j=0;j<G->uliColumns[kk];j++) fprintf(G->log_fPtr[kk],"%d\t",G->iVar[kk][i][j]);
  	        fprintf(G->log_fPtr[kk],"\n");
  	}
  	
  	
  	
  	//Allocate the buffer for the output
/*  	
	G->ldM_matrix[kk] = (LD **) malloc(sizeof(void*) * (G->buffer_lines));
  	for(i=0;i<G->buffer_lines;i++) G->ldM_matrix[kk][i] = (LD *) malloc (sizeof(LD) * (G->uliColumns[kk]+1));
*/	
	G->ldM_matrix[kk]   = (LD *) malloc(sizeof(LD) * (G->uliColumns[kk]+1)*(G->buffer_lines));

	buffer_size += sizeof(LD)*(G->uliColumns[kk]+1)*(G->buffer_lines);
		
	if(kk == 0) {
	      G->ldM_send = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
	      //G->ldM_recv = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
	}

  	//Write on the output file

	t1 = (LD) G->uliIndexCol+1;
	fwrite(&t1,sizeof(long double),1,G->m_fPtr[kk]);

	add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
	G->uliRow_counter[kk]++;
  	
	//Qui non dovrei avere problemi di buffer pieno, pero'...
	if(G->uliRow_counter[kk] == (ULI) G->buffer_lines) {
  		write_matrix(G->m_fPtr[kk],G->ldM_matrix[kk],G->uliIndexCol,kk,G);
  		G->uliRow_counter[kk] = 0;
  	}
  	
  	getFS(kk,G);
	//fprintf(G->fs_fPtr[kk], "%LG\n", G->ldFS[kk]);
  }

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  rtime = elapsed/1000000;
  
  gettimeofday(&tstart,NULL);
  
  /********************************************/
  /********************************************/
  /********************************************/
  /*Iterative loop of the stochastic algorithm*/
  /********************************************/
  /********************************************/
  /********************************************/
  while (t < uliStep_max) {									//98,9%
  //printf("Tempo %10LG\n",t);									

  	if(uliStep % (ULI) G->every == 0) {	
		printf("%10LG ",t);
		if(ee == 8) { printf("\n"); ee=0; }
		else ee++;
		//fflush(NULL);									
	}

  	mymintau = LONG_MAX;
	
	//EACH membrane
  	for(k=nMstart;k<=nMstop;k++) {								//94%							
  	        kk = k-nMstart; 
  	        
  	        get_HOR(kk,G);	//Compute the highest order of reaction				//37% - 				
  	        //Backup M and FS 
  	        memcpy(G->Mbkp[kk], G->M[kk], G->uliColumns[kk] * sizeof(LD));			
  	        G->ldFSbkp[kk] = G->ldFS[kk];
  	        
		//printf("%d: %d\n",k,G->iFlag_SSA[kk]); fflush(NULL);
		//Dynamic step
  	        stepkind[kk] = dpp_step1(r, kk, a0c, t2, t, uliStep, G);			//56% 			
  	        if(G->tau[kk] < mymintau) mymintau = G->tau[kk];
		
		//printf("%d: %LG %d\n",k,G->tau[kk],stepkind[kk]); fflush(NULL);
		
  	}

	//ALL membranes
	globmintau = get_min_tau(mymintau);

	cleanbuffers(nMembranes,G);

	//printf("%d: step1\n",id); fflush(NULL);

	//EACH membrane	
  	for(k=nMstart;k<=nMstop;k++) {
		kk = k-nMstart;
  	        G->tau[kk] = globmintau;
  	        dpp_step2(r, kk, stepkind[kk], a0c, t2, G);
	}

	//printf("%d: step2\n",id); fflush(NULL);	

	//ALL membranes
	receive(nMstart,nMstop,nMembranes, G);
	do {
		for(k=nMstart;k<=nMstop;k++) getFS(k-nMstart, G);
		
		cn = checkNegFS(nMstart, nMstop, G);
		
		cleanbuffers(nMembranes, G);
				
		if(cn == 1) {
			for(k=nMstart;k<=nMstop;k++) dpp_step3(r,k-nMstart, G);

			//ALL membranes
			receive(nMstart,nMstop,nMembranes, G);
		}
	} while(cn == 1);

	//printf("%d: step3\n",id); fflush(NULL);	

	t3 = (LD)0.0;

	uliStep++;

	for(k=nMstart;k<=nMstop;k++) {								//6.5 %%
		kk = k-nMstart;
		
		t1 = dpp_step4(kk, t,uliStep_max, uliStep, G);					
		
		//printf("%d %LG %LG\n",k,t1,t3);
		
		if(t3 < t1) t3 = t1;
		
  		//fprintf(G->fs_fPtr[kk], "%LG\n", G->ldFS[kk]);					//6.4%
	}
	t = t3;
  }

  /********************************************************
   ********************************************************
   ***				FINE			***
   ********************************************************
   ********************************************************/
  //printf("%d: step4\n",id); fflush(NULL);

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  ctime = (elapsed/1000000)-G->iotime;
  
  printf ("\n\nData acquisition %d-%d in %.2f\n",nMstart,nMstop,rtime);
  printf ("Data processing  %d-%d in %.2f (%.2e per step)\n",nMstart,nMstop,ctime, ctime/uliStep);

  ctime = G->iotime;
  
  gettimeofday(&tstart,NULL);
  
  for(k=nMstart;k<=nMstop;k++) {
  	kk = k-nMstart;
/*
  	// Add a vector to the output buffer 
  	if (uliStep % (ULI) G->every != 0){
		add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
  		G->uliRow_counter[kk]++;
  	}
*/
	//Write always the last status
	add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
	
  	// Write in the output file what is currently in the buffer 
  	write_matrix(G->m_fPtr[kk],G->ldM_matrix[kk],G->uliIndexCol,kk,G);

	
  	//fclose(G->fs_fPtr[kk]);	
	fclose(G->m_fPtr[kk]); fclose(G->log_fPtr[kk]); 
	fclose(G->rprob_fPtr[kk]);
  } 
  
  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  ctime += (elapsed/1000000);
  
  printf("Output production in     %.2f\n",ctime);
  
  printf("\nEnd: %d steps, %d write operations with a buffer of %4.1f MB\n",\
  	uliStep,(uliStep/G->every)/G->buffer_lines+1,(float)buffer_size/1048576);
  
}



int main( int argc, char *argv[]) {
  LD step_max=10;
  FILE *fPtr;

  /*Random Number Generator definition*/
  gsl_rng *r;
  const gsl_rng_type *T_gsl;
  //struct rlimit rll;
  struct globale G;
  struct timeval tstart, tstop;
  float elapsed;
  
  
  gettimeofday(&tstart,NULL);
  
  G.uliMAXSEED = pow(2,32)-1;
  G.uliIndexCol=0;
  G.every=0.0;
  G.iTgtcheck=0;
  
  /*Random Number Generator initialization*/
  gsl_rng_env_setup();
  T_gsl = gsl_rng_default;
  r = gsl_rng_alloc(T_gsl);


  /* Read the data common to al the processes */
  system("mkdir -p output/");
  
  /*Read the maximal time of the simulation*/
  if ((fPtr = fopen ("input/time_max","r")) == NULL)	  {
  	  printf ("\nCannot open time max\n");
  	  exit(0);
  }
  else {
  	  fscanf(fPtr,"%LG",&step_max);
  	  printf("Time max: %LG ",step_max);
  	  fclose(fPtr);
  }
  
  
  /*Read the every input file*/
  if ((fPtr = fopen ("input/every", "r")) == NULL) {
  	  G.every = 1;
  	  G.buffer_lines=10000;
  	  printf ("\nCannot open every - loading default %d and %d\n",G.every,G.buffer_lines);
  }
  else {
  	  fscanf(fPtr,"%d", &G.every);
  	  fscanf(fPtr,"%d", &G.buffer_lines);
  	  printf ("Every: %d , buffer_lines: %d\n\n",G.every, G.buffer_lines);
  	  fclose (fPtr);
  }


  /* Do the job */
  initialize(step_max, r, &G);

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  elapsed /= 1000000;



  printf("End of the program in %.2f sec.\n", elapsed); fflush(NULL);
  return (0);
}