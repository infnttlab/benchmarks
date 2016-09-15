/*
** Caso lungo **
Data acquisition 0-80 in 0.02
Data processing  0-80 in 528.32 (2.00e-04 per step)
Output production in     0.17

End: 2647727 steps, 27 write operations with a buffer of 61.8 MB
End of the program in 528.89 sec.

** Caso breve **
Data acquisition 0-80 in 0.28
Data processing  0-80 in 19.07 (1.81e-04 per step)
Output production in     0.02

End: 105216 steps, 2 write operations with a buffer of 61.8 MB
End of the program in 19.76 sec.
*/

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
  char file[60],app[20];
  FILE* fPtr;
  int nMembranes, nMstart, nMstop, nMnum, k, kk;
  struct timeval tstart, tstop;
  float elapsed, ctime, rtime;
  int* stepkind;
  LD *a0c, *t2;
  int cn;
  LD t1, t3;
  LD t;
  int uliStep=0;
  int ee, indaux;
  int buffer_size =0;
  unsigned*** unLeft_side, ***unRight_side;
  int ***iTgt_matrix;
  

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
  
  nMnum = nMembranes;
  nMstart = 0;
  nMstop = nMembranes-1;

  unLeft_side	        = (unsigned***) malloc(nMnum * sizeof(void*)); //temporanei
  unRight_side          = (unsigned***) malloc(nMnum * sizeof(void*));
  iTgt_matrix		= (int ***) malloc(nMnum * sizeof(void*));
    
  ////restano costanti
  G->side_index		= (int*) malloc((nMnum+1) * sizeof(int)); //costanti
  G->side_index2	= (int*) malloc((nMnum+1) * sizeof(int)); //costanti
  
  G->ldMembSize		= (LD *)  malloc(nMnum * sizeof(LD));  
/*
  G->ldM_feed		= (LD **) malloc(nMnum * sizeof(void*));
  G->ldMolSize		= (LD **) malloc(nMnum * sizeof(void*));
  G->ldC_vector		= (LD **) malloc(nMnum * sizeof(void*));
*/
  G->iTgt_vector	= (int **)  malloc(nMnum * sizeof(void*));

  G->uliRows		= (ULI*) malloc(nMnum * sizeof(ULI));
  G->uliColumns		= (ULI*) malloc(nMnum * sizeof(ULI)); // colonne == specie quindi in realta' e' uno scalare

  ////ausiliari 
  G->m_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*)); // per i due file di output
  G->log_fPtr		= (FILE **) malloc(nMnum * sizeof(FILE*));

  G->ldM_matrix		= (LD **) malloc(nMnum * sizeof(void*));  //buffer e contatore
  G->uliRow_counter	= (ULI*) malloc(nMnum * sizeof(ULI));

  ////realmente private e cambiano
  G->M			= (LD **) malloc(nMnum * sizeof(void*)); 
  G->Mbkp  		= (LD **) malloc(nMnum * sizeof(void*)); 

  G->uliK_rule		= (ULI**) malloc(nMnum * sizeof(void*)); 
  G->uliCritical	= (ULI**) malloc(nMnum * sizeof(void*));
/*
  G->ldR_prob		= (LD **) malloc(nMnum * sizeof(void*));  
  G->ldR_prob_c		= (LD **) malloc(nMnum * sizeof(void*));  
*/
  G->ldFS  		= (LD *) malloc(nMnum * sizeof(LD));
  G->ldFSbkp		= (LD *) malloc(nMnum * sizeof(LD));
  
  G->tau		= (LD *) malloc(nMnum * sizeof(LD));
  
  stepkind		= (int*) malloc(nMnum * sizeof(int));//passaggio variabili
  a0c  			= (LD *) malloc(nMnum * sizeof(LD));
  t2  			= (LD *) malloc(nMnum * sizeof(LD));
  
  G->iFlag_SSA		= (int*) malloc(nMnum * sizeof(int)); 		
  G->iFlagStep		= (int*) malloc(nMnum * sizeof(int));
  G->ldTau_ssa		= (LD *) malloc(nMnum * sizeof(LD));

////////////////////////verificare
  
  //G->checkFS		= (int*) malloc(nMnum * sizeof(int));
  //G->uliR_app		= (ULI**) malloc(nMnum * sizeof(void*));
  //G->uliOrder		= (ULI**) malloc(nMnum * sizeof(void*)); //usata solo in get_HOR ?????
  //G->indexes		= (int **) malloc(nMnum * sizeof(void*));//inutile, o meglio letto tutte le volte solo per scriverlo nel log file

  //G->ldHOR 		= (LD **) malloc(nMnum * sizeof(void*));  //vedere se get_HOR puo' essere chiamata in dpp_step1 e quindi eliminata questa var

////////////////////////
  
  for(i=0;i<nMnum;i++) { G->tau[i] = 0.0; G->iFlag_SSA[i] = 0; G->ldTau_ssa[i] = 0.0; G->uliRow_counter[i] = 0;}

  gettimeofday(&tstart,NULL);

  t = 0.0;
  
  G->iotime = 0.0;
  G->side_index[0] = 0;
  G->side_index2[0] = 0;
    
  for(k=nMstart;k<=nMstop;k++) {  
  	kk = k-nMstart;

  	/************
  	 ************
  	 OUTPUT FILES
  	 ************
  	 ************/
  	
  	//Create a log file name labelled with the process id
  	strcpy(file,"output/log_");
	sprintf(app,"%d",k);
	
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
  	  
  	/***********
  	 ***********
  	 INPUT FILES
  	 ***********
  	 ***********/
  	/* Read the input files */
  	G->iTgtcheck = load_files1(k,nMstart,G, unLeft_side, unRight_side, iTgt_matrix);
	
  	//G->uliR_app[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Application vector*/
  	//G->uliOrder[kk]    = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Order of reactions vector*/
  	//G->ldHOR[kk]	   = (LD *) malloc(sizeof(LD) * G->uliColumns[kk]);	/* Highest order of reaction vector*/
		
  	G->uliK_rule[kk]   = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Tau leaping application vector*/
  	G->uliCritical[kk] = (ULI*) malloc(sizeof(ULI) * G->uliRows[kk]);	/* Critical reactions vector*/

  	//G->ldR_prob[kk]    = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	/* Propensity functions vector*/
  	//G->ldR_prob_c[kk]  = (LD *) malloc(sizeof(LD) * G->uliRows[kk]);	/* Critical reactions propensity functions vector*/
  	
  	
  	//Allocate the buffer for the output
	G->ldM_matrix[kk]   = (LD *) malloc(sizeof(LD) * (G->uliColumns[kk]+1)*(G->buffer_lines));

	buffer_size += sizeof(LD)*(G->uliColumns[kk]+1)*(G->buffer_lines);
		
	if(kk == 0) {
	      G->ldM_send = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
	      //G->ldM_recv = (LD *) malloc(nMembranes*G->uliColumns[kk]*sizeof(LD));
	}

  	//Write on the output file
	t1 = (LD) G->uliIndexCol+1;
	fwrite(&t1,sizeof(long double),1,G->m_fPtr[kk]);

	
  	//Write on the output file
  	//write_M(kk,t,G);

	add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
  	
	//Qui non dovrei avere problemi di buffer pieno, pero'...
	if(G->uliRow_counter[kk] == (ULI) G->buffer_lines) {
  		write_matrix(kk,G);
  	}
	//fprintf(G->fs_fPtr[kk], "%LG\n", G->ldFS[kk]);
  }

  G->unLeft_side  = (unsigned*) malloc(G->side_index[nMnum]*sizeof(unsigned));
  G->unRight_side = (unsigned*) malloc(G->side_index[nMnum]*sizeof(unsigned));
  G->iVar  	  = (int *)     malloc(G->side_index[nMnum] * sizeof(int));    
  G->iVarSend	  = (int *)     malloc(G->side_index[nMnum] * sizeof(int));
  G->iTgt_matrix  = (int *)     malloc(G->side_index[nMnum] * sizeof(int));
  
  G->ldM_feed	  = (LD *) malloc(nMnum*G->uliColumns[0]*sizeof(LD));
  G->ldMolSize	  = (LD *) malloc(nMnum*G->uliColumns[0]*sizeof(LD));
  
  G->ldC_vector	  = (LD *) malloc(G->side_index2[nMnum]*sizeof(LD));
  G->ldR_prob	  = (LD *) malloc(G->side_index2[nMnum]*sizeof(LD));
  G->ldR_prob_c	  = (LD *) malloc(G->side_index2[nMnum]*sizeof(LD));

  for(i=0;i<G->side_index[nMnum]; i++) { G->iVar[i] = 0; G->iVarSend[i] = 0; }
    
   //From matrixes to vectors
  for(k=nMstart;k<=nMstop;k++) { 
  	load_files2(k,nMstart,G);
   
  	kk = k-nMstart;
	//fprintf(G->log_fPtr[kk],"modifying the left_side\n");
  	
	//Print the variations matrix into the log file
  	fprintf(G->log_fPtr[kk], "\nVariations Matrix\n");
	
	if(G->iTgtcheck != 1) {
		for(i=0;i<G->uliRows[kk];i++) {
			for(j=0;j<G->uliColumns[kk];j++) 
				G->iTgt_matrix[G->side_index[kk]+(i*G->uliColumns[kk])+j] = iTgt_matrix[kk][i][j];
			free(iTgt_matrix[kk][i]);
		}
		free(iTgt_matrix[kk]);
	}
	
	for(i=0;i<G->uliRows[kk];i++) {
		for(j=0;j<G->uliColumns[kk];j++) {
			indaux = G->side_index[kk]+(i*G->uliColumns[kk])+j;

			G->unLeft_side [indaux] =  unLeft_side[kk][i][j];
			G->unRight_side[indaux] = unRight_side[kk][i][j];
						
			if(G->iTgtcheck == 1) {
				if(G->iTgt_vector[kk][i] == -1){
  	        			G->iVar[indaux] -= unLeft_side[kk][i][j];
  	        			G->iVar[indaux] += unRight_side[kk][i][j];
  	        		} //IF I have to send the output to another process...
  	        		else    G->iVar[indaux] -= unLeft_side[kk][i][j];	
			}
			else {
				if(G->iTgt_matrix[indaux] == k){ //era id, quando id == numero della membrana
  	        			G->iVar[indaux] -= unLeft_side [kk][i][j];
  	        			G->iVar[indaux] += unRight_side[kk][i][j];
  	        			G->iVarSend[indaux] = G->iVar[indaux];
  	        		}
  	        		else{
  	        			G->iVar[indaux]    -= unLeft_side [kk][i][j];
  	        			G->iVarSend[indaux] = unRight_side[kk][i][j];
  	        		}			
			}
			fprintf(G->log_fPtr[kk],"%d\t",G->iVar[indaux]);		
		}
		free(unLeft_side [kk][i]);
		free(unRight_side[kk][i]);
		fprintf(G->log_fPtr[kk],"\n");
	}
	free(unLeft_side[kk]);
	free(unRight_side[kk]);
  }
  free(unLeft_side);
  free(unRight_side);
  free(iTgt_matrix);
  
  for(k=nMstart;k<=nMstop;k++) getFS(k-nMstart,G);

  //E da qui iVar e iVarSend restano costanti, cosi' come unLeft e unRight side
  //iVarSend solo se tgt_matrix e usato in dpp_update_tgt_matrix

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  rtime = elapsed/1000000;
  
  gettimeofday(&tstart,NULL);
  
  ee = 0;
  /********************************************/
  /********************************************/
  /********************************************/
  /*Iterative loop of the stochastic algorithm*/
  /********************************************/
  /********************************************/
  /********************************************/
  while (t < uliStep_max) {									
  	if(uliStep % (ULI) G->every == 0) {	
		printf("%12LG ",t);
		if(ee == 8) { printf("\n"); ee=0; }
		else ee++;
		//fflush(NULL);									
	}

  	mymintau = LONG_MAX;
	
	//EACH membrane
  	for(k=nMstart;k<=nMstop;k++) {									//94%%
  	        kk = k-nMstart; 
  	        
  	        //get_HOR(kk,G); //Compute the highest order of reaction, direttamente in dpp_step1	//di cui 46.6%	
		 
  	        //Backup M and FS 
  	        memcpy(G->Mbkp[kk], G->M[kk], G->uliColumns[kk] * sizeof(LD));				//e 1% qui
  	        G->ldFSbkp[kk] = G->ldFS[kk];
  	        
		//printf("%d: %d\n",k,G->iFlag_SSA[kk]); fflush(NULL);
		
		//Dynamic step
  	        stepkind[kk] = dpp_step1(r, kk, a0c+kk, t2+kk, t, uliStep, G); 				//e qui l'altro 46.1%
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
  	        //dpp_step2(r, uliStep_max, kk, a0c+kk, t2+kk, G, nMembranes);
		dpp_step2(r, kk, stepkind[kk], a0c+kk, t2+kk, G, nMembranes);
	}

	//printf("%d: step2\n",id); fflush(NULL);	

	//ALL membranes
	receive(nMstart,nMstop,nMembranes, G);

	do {
		//EACH membrane
		for(k=nMstart;k<=nMstop;k++) getFS(k-nMstart, G);

		//ALL membranes
		cn = checkNegFS(nMstart, nMstop, G);
		
		cleanbuffers(nMembranes, G);
				
		if(cn == 1) {
			//EACH membrane
			for(k=nMstart;k<=nMstop;k++) dpp_step3(r, k-nMstart, G, nMembranes);

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
	//Write always the last status
	add_vector_in_matrix(G->ldM_matrix[kk],G->M[kk],G->uliColumns[kk],t,kk,G);
	
  	// Write in the output file what is currently in the buffer 
  	write_matrix(kk,G);
	
	fclose(G->m_fPtr[kk]); 
	fclose(G->log_fPtr[kk]); 
  } 
  
  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  ctime += (elapsed/1000000);
  
  printf("Output production in     %.2f\n",ctime);
  
  printf("\nEnd: %d steps, %d write operations with a buffer of %4.1f MB\n",\
  	uliStep,(uliStep/G->every)/G->buffer_lines+1,(float)buffer_size/1048576);
  
}



int main(int argc, char *argv[]) {
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
  
  //Random Number Generator initialization
  gsl_rng_env_setup();
  T_gsl = gsl_rng_default;
  r = gsl_rng_alloc(T_gsl);


  //Read the data common to al the processes
  system("mkdir -p output/");
  
  //Read the maximal time of the simulation
  if ((fPtr = fopen ("input/time_max","r")) == NULL)	  {
  	  printf ("\nCannot open time max\n");
  	  exit(0);
  }
  else {
  	  fscanf(fPtr,"%LG",&step_max);
  	  printf("Time max: %LG ",step_max);
  	  fclose(fPtr);
  }
  
  
  //Read the every input file
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

  //Do the job
  initialize(step_max, r, &G);

  gettimeofday(&tstop,NULL);

  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  elapsed /= 1000000;

  printf("End of the program in %.2f sec.\n", elapsed); fflush(NULL);
  return (0);
}
