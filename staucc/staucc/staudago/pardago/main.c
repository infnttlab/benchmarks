//Occhio, potrebbe servire ulimit -s unlimited  (di default 10 M)

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


#include <sys/resource.h>
#include <sys/time.h> 

#include "io_file.h"
#include "manage_prob.h"

/***** GLOBALS *****
//FILE
FILE **m_fPtr, **log_fPtr, **fs_fPtr, **rprob_fPtr;

ULI *uliRows, *uliColumns;
ULI **uliK_rule, **uliR_app;
ULI **uliOrder, **uliCritical, *uliRow_counter;

LD **ldC_vector, *tau, **M, **Mbkp, **ldM_feed, **ldR_prob;
LD **ldR_prob_c, *ldTau_ssa;
LD **ldHOR;

LD ***ldM_matrix;
LD * ldM_send, *ldM_recv;

LD **ldMolSize, *ldMembSize, *ldFS, *ldFSbkp;
int  **iTgt_vector, ***iTgt_matrix;
int  ***iVar, **indexes, *iFlag_SSA, *iFlagStep, ***iVarSend;
int* checkFS;
unsigned ***unLeft_side, ***unRight_side;
struct entry **LSlist, **RSlist;

ULI uliIndexCol=0, uliMAXSEED;
LD every=0.0, buffer_lines;
*******************/

//void dpp_update_tgt_vector(int, gsl_rng *, int, struct globale*);
//void dpp_update_tgt_matrix(int, gsl_rng *, int, struct globale*);


//void (*dpp_update)(int, gsl_rng *, int, struct globale*);

/*Initialize the MPI process*/
void initialize(LD uliStep_max, int size, gsl_rng *r, struct globale *G) {	
  ULI i=0, j=0;
  LD rnd=0, mymintau, globmintau;
  char file[20],app[20];
  int id, count, iRows_index=0;
  MPI_Status status[size];
  struct entry *curr_node=NULL;
  FILE* fPtr;
  int nMembranes, nMstart, nMstop, nMnum, k, kk;
  int* assignments;
  struct timeval tstart, tstop;
  float elapsed;
  int* stepkind;
  LD *a0c, *t2;
  int cn;
  LD t1, t3;
  LD t;
  int uliStep=0;
  //int iTgtcheck=0;

  /*Get the rank of the process*/
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  /** NOTE: Size is the number of the parallel processes, NOT 
   * the number of membranes
   */

  /*Process zero set the seed of the pseudo-random chain and send a random number to every process*/
  if (id == 0) {
  	  gsl_rng_set(r, time(NULL));
  	  ///		  
  	  //		  printf("%d received %e\n",id,rnd);
  	  
  	  for (i=1; i < size; i++){
  		  //rnd = gsl_rng_uniform_int (r, ULONG_MAX);
  		  rnd = gsl_rng_uniform_int (r, G->uliMAXSEED);
  		  MPI_Send(&rnd, 1, MPI_LONG_DOUBLE, i, 1, MPI_COMM_WORLD);
  	  }
  }
  /*The other processes receive a random number that will be used to set their own pseudo-random chain*/
  else {
  	  MPI_Recv(&rnd, 1, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD, status);
  	  gsl_rng_set(r, rnd);
  	  ///	  
  	  //		  printf("%d received %e\n",id,rnd);
  }


  if ((fPtr = fopen ("input/numMembranes.txt","r")) == NULL)	  {
  	  printf ("\nCannot obtain the numbers of membranes, the file input/numMembranes.txt is required\n");
  	  exit(0);
  }
  else {
  	  fscanf(fPtr,"%d",&nMembranes);
  	  fclose(fPtr);
  }

  assignments = (int *) malloc(nMembranes*sizeof(int));
  nMnum = nMembranes%size;

  for(i=0;i<size;i++) {
  	  nMstart=(nMembranes/size*i)+((i>=nMnum) ? nMnum : i);
  	  nMstop=(nMembranes/size*(i+1))+((i>=nMnum) ? nMnum-1 : i);	  
  	  for(j=nMstart;j<=nMstop;j++) assignments[j]=i;
  }


  nMstart=(nMembranes/size*id)+((id>=nMnum) ? nMnum : id);
  nMstop=(nMembranes/size*(id+1))+((id>=nMnum) ? nMnum-1 : id);   
  nMnum = nMstop-nMstart+1;


  //	  if(!id) for(i=0;i<nMembranes;i++) printf("%d : %d\n",i,assignments[i]);
  //	  MPI_Barrier(MPI_COMM_WORLD); exit(1);
  //	  printf("%3d : from %3d to %3d (%3d) of %3d\n",id,nMstart,nMstop,nMnum,nMembranes);
  	  
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
  G->uliCritical		= (ULI**) malloc(nMnum * sizeof(void*));
  G->uliRow_counter	= (ULI*) malloc(nMnum * sizeof(ULI));
  G->ldR_prob		= (LD **) malloc(nMnum * sizeof(void*));  
  G->ldHOR 		= (LD **) malloc(nMnum * sizeof(void*));  //resettato in get_HOR
  G->ldR_prob_c		= (LD **) malloc(nMnum * sizeof(void*));  

  G->iVar  		= (int ***) malloc(nMnum * sizeof(void*));	
  G->iVarSend		= (int ***) malloc(nMnum * sizeof(void*));
  
  G->ldM_matrix		= (LD ***) malloc(nMnum * sizeof(void*));

  G->ldFS  		= (LD *) malloc(nMnum * sizeof(LD));
  G->ldFSbkp		= (LD *) malloc(nMnum * sizeof(LD));
  G->tau			= (LD *) malloc(nMnum * sizeof(LD));

  G->m_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  G->log_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  //G->fs_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));
  //G->rprob_fPtr	= (FILE **) malloc(nMnum * sizeof(FILE*));
  
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
  
  for(k=nMstart;k<=nMstop;k++) {  
  	  kk = k-nMstart;
  	  
  	  /************
  	   ************
  	   OUTPUT FILES
  	   ************
  	   ************/
  	  
  	  /*Create a log file name labelled with the process id*/
  	  strcpy(file,"output/log_");
  	  itoa(k,app);
  	  strcat(file,app);
  	  
  	  /*Create a new log file into the 'output' folder*/
  	  if ((G->log_fPtr[kk] = fopen (file,"w")) == NULL) {
  		  printf ("\nCannot write in the current directory - Log file\nAbort\n");
  		  exit(0);
  	  }
  	  
  	  /*Write process information into the log file*/
  	  fprintf(G->log_fPtr[kk],"My id is %d\n",id);
  	  fprintf(G->log_fPtr[kk],"uliStep_max %LG\n",uliStep_max);
  	  fprintf(G->log_fPtr[kk],"Size %d\n",size);
  	  
  	  /*Create an output file name labelled with the process id*/
  	  strcpy(file,"output/multi_");
  	  strcat(file,app);
  	  if ((G->m_fPtr[kk] = fopen (file,"w")) == NULL) {
  		  printf ("\nCannot write in the current directory - Output file\nAbort\n");
  		  exit(0);
  	  }
  	  
	/*
  	  strcpy(file,"output/fs_");
  	  strcat(file,app);
  	  if ((G->fs_fPtr[kk] = fopen (file,"w")) == NULL) {
  		  printf ("\nCannot write in the current directory - fs file\nAbort\n");
  		  exit(0);
  	  }
  	  
  	  //create the file to store the propensity function values
  	  strcpy(file,"output/rprob_");
  	  strcat(file,app);
  	  if ((G->rprob_fPtr[kk] = fopen (file,"w")) == NULL){
  		  printf ("\nCannot write in the current directory - Log file\nAbort\n");
  		  exit(0);
  	  }

  	  */

  	  /***********
  	   ***********
  	   INPUT FILES
  	   ***********
  	   ***********/
  	  
  	  /* Read the input files */
  	  G->iTgtcheck = load_files(k,nMstart,G);
	  

  	  /*Compile the lists with the reactions*/
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
  	  
  	  /* Left-hand side */
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
  	  
  	  /*Print the lists into the log file*/
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
  	  
  	  G->uliR_app[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	  /* Application vector*/
  	  G->uliK_rule[kk]   = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	  /* Tau leaping application vector*/
  	  G->uliOrder[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	  /* Order of reactions vector*/
  	  G->uliCritical[kk] = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	  /* Critical reactions vector*/
  	  G->ldR_prob[kk]    = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	  /* Propensity functions vector*/
  	  G->ldHOR[kk]	     = (LD *) malloc(sizeof(LD) * G->uliColumns[kk]);	  /* Highest order of reaction vector*/
  	  G->ldR_prob_c[kk]  = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	  /* Critical reactions propensity functions vector*/
  	  G->iVar[kk]	     = (int **) malloc(G->uliRows[kk] * sizeof(void *) );/* Variations matrix*/
  	  G->iVarSend[kk]    = (int **) malloc(G->uliRows[kk] * sizeof(void *) );
  	  
  	  for(i=0;i<G->uliRows[kk];i++){
  		  G->iVar[kk][i]	  = malloc(G->uliColumns[kk]* sizeof(int) );
  		  G->iVarSend[kk][i]      = malloc(G->uliColumns[kk]* sizeof(int) );
  	  
		  for(j=0; j<G->uliColumns[kk];j++) {
			  G->iVar[kk][i][j]     = 0;
			  G->iVarSend[kk][i][j] = 0;
		  }
	  }
  	  /* Fill the variation matrix */
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
  				  if(G->iTgt_matrix[kk][i][j] == id){
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
  	  
  	  /*Print the variations matrix into the log file*/
	  fprintf(G->log_fPtr[kk], "\nVariations Matrix\n");
  	  for(i=0;i<G->uliRows[kk];i++){
  		  for(j=0;j<G->uliColumns[kk];j++) fprintf(G->log_fPtr[kk],"%d\t",G->iVar[kk][i][j]);
  		  fprintf(G->log_fPtr[kk],"\n");
  	  }
  	  
  	  ///	  
  	  //	  printf("%d QUI %d %d %LG %LG\n",id,uliRows[kk],uliColumns[kk],buffer_lines,uliStep_max); fflush(NULL);
  	  //	  MPI_Barrier(MPI_COMM_WORLD);
  	  
  	  
  	  /*Allocate the buffer for the output*/
  	  G->ldM_matrix[kk] = (LD **) malloc(sizeof(void*) * (G->buffer_lines));
  	  for(i=0;i<G->buffer_lines;i++) G->ldM_matrix[kk][i] = (LD *) malloc (sizeof(LD) * (G->uliColumns[kk]+1));
	  
	  if(kk == 0) {
	  	G->ldM_send = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
		G->ldM_recv = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
	  }
	  
/*
  	  // Buffer matrix 
  	  ldM_send[kk] = (LD **) malloc(size * sizeof(LD *) );
  	  for(i=0;i<size;i++) ldM_send[kk][i] = (LD *) malloc(uliColumns[kk] * sizeof(LD) );
  	  
  	  //Buffer vectors
  	  ldM_sendV[kk] = malloc(uliColumns[kk] * sizeof(LD));
  	  ldDelta_M[kk] = malloc(uliColumns[kk] * sizeof(LD));
*/  	  
  	  /*Write on the output file*/
  	  write_M(kk,t,G);
  	  
  	  getFS(kk,G);
  	  //fprintf(G->fs_fPtr[kk], "%LG\n", G->ldFS[kk]);
  }

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  printf ("%d: Data acquisition %d-%d in %.2f\n",id,nMstart,nMstop,elapsed/1000000);
  
  /********************************************/
  /********************************************/
  /********************************************/
  /*Iterative loop of the stochastic algorithm*/
  /********************************************/
  /********************************************/
  /********************************************/
  while (t < uliStep_max) {
  	if(id==0) if(uliStep % (ULI) G->every == 0) printf("%LG\n",t);

  	mymintau = LONG_MAX;
	
	//EACH membrane
  	for(k=nMstart;k<=nMstop;k++) {
  	        kk = k-nMstart; 
  	        
  	        get_HOR(kk,G);						/* Compute the highest order of reaction */
  	        
  	        memcpy(G->Mbkp[kk], G->M[kk], G->uliColumns[kk] * sizeof(LD));	/* Backup M and FS */
  	        G->ldFSbkp[kk] = G->ldFS[kk];
  	        
		//printf("%d: %d\n",k,G->iFlag_SSA[kk]); fflush(NULL);
		
  	        stepkind[kk] = dpp_step1(size, r, kk, a0c, t2, t, uliStep, G); /* Dynamic step */
  	        if(G->tau[kk] < mymintau) mymintau = G->tau[kk];
		
		//printf("%d: %LG %LG %d\n",k,tau[kk],mymintau,iFlag_SSA[kk]); fflush(NULL);
		
  	}

	//ALL membranes
	globmintau = get_min_tau(mymintau);

	cleanbuffers(nMembranes,G);

	//printf("%d: step1\n",id); fflush(NULL);

	//EACH membrane	
  	for(k=nMstart;k<=nMstop;k++) {
		kk = k-nMstart;
  	        G->tau[kk] = globmintau;
  	        dpp_step2(size, r, kk, stepkind[kk], a0c, t2, G);
	}

	//printf("%d: step2\n",id); fflush(NULL);	

	//ALL membranes
	receive(nMstart,nMstop,nMembranes, G);

	do {
		for(k=nMstart;k<=nMstop;k++) getFS(k-nMstart, G);
		
		cn = checkNegFS(nMstart, nMstop, G);
		
		cleanbuffers(nMembranes, G);
				
		if(cn == 1) {
			for(k=nMstart;k<=nMstop;k++) dpp_step3(size,r,k-nMstart, G);

			//ALL membranes
			receive(nMstart,nMstop,nMembranes, G);
		}
	} while(cn == 1);

	//printf("%d: step3\n",id); fflush(NULL);	

	t3 = (LD)0.0;

	uliStep++;

	for(k=nMstart;k<=nMstop;k++) {
		kk = k-nMstart;
		
		t1 = dpp_step4(kk, t,uliStep_max, uliStep, G);
		
		//printf("%d %LG %LG\n",k,t1,t3);
		
		if(t3 < t1) t3 = t1;
		
  		//fprintf(G->fs_fPtr[kk], "%LG\n", G->ldFS[kk]);
	}
	t = t3;
  }

  //printf("%d: step4\n",id); fflush(NULL);


  for(k=nMstart;k<=nMstop;k++) {
  	kk = k-nMstart;

  	// Add a vector to the output buffer 
  	if (uliStep % (ULI) G->every != 0){
		add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
  		G->uliRow_counter[kk]++;
  	}

  	// Write in the output file what is currently in the buffer 
  	write_matrix(G->m_fPtr[kk],G->ldM_matrix[kk],G->uliIndexCol,kk,G);

	
  	//fclose(G->fs_fPtr[kk]);	
	fclose(G->m_fPtr[kk]); 
	fclose(G->log_fPtr[kk]); 
	//fclose(G->rprob_fPtr[kk]);
  }

}



int main( int argc, char *argv[]) {
  LD step_max=10, array3[3];
  int id, size;
  FILE *fPtr;

  /*Random Number Generator definition*/
  gsl_rng *r;
  const gsl_rng_type *T_gsl;
  //struct rlimit rll;
  struct globale G;
  
  G.uliMAXSEED = pow(2,32)-1;
  G.uliIndexCol=0;
  G.every=0.0;
  G.iTgtcheck=0;
  
  /*Random Number Generator initialization*/
  gsl_rng_env_setup();
  T_gsl = gsl_rng_default;
  r = gsl_rng_alloc(T_gsl);

  /*MPI Processes initialization*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

//  dpp_update = dpp_update_tgt_vector;

  //Potrebbe servire
/*	
  rll.rlim_cur = RLIM_INFINITY; 
  rll.rlim_max = RLIM_INFINITY;
  setrlimit(RLIMIT_STACK,&rll);
*/  
  /* Read the data common to al the processes */
  if(id==0) {
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
//  		  fscanf(fPtr,"%LG", &G.every);
//  		  fscanf(fPtr,"%LG", &G.buffer_lines);

  		  fscanf(fPtr,"%d", &G.every);
  		  fscanf(fPtr,"%d", &G.buffer_lines);
  		  printf ("Every: %d , buffer_lines: %d\n",G.every, G.buffer_lines);
  		  fclose (fPtr);
  	  }
  	  
  	  array3[0] = step_max; array3[1] = (LD)G.every; array3[2] = (LD)G.buffer_lines;
  	  
  }

  /*Broadcast the time value*/
  MPI_Bcast(array3,3,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);

  if(id != 0) {
  	  step_max = array3[0]; G.every = (int)array3[1];  G.buffer_lines = (int)array3[2];
  }

  /* Do the job */
  initialize(step_max, size, r, &G);

  printf("%d: fine\n",id); fflush(NULL);
  MPI_Barrier(MPI_COMM_WORLD);
  
  /*Finalizing the communicator*/
  MPI_Finalize();

  return (0);
}
