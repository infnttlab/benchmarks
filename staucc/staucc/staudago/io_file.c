/*
Data acquisition 0-80 in 1.66
Data processing  0-80 in 63.48 (1.27e-02 per step)
==1059== Syscall param write(buf) points to uninitialised byte(s)
==1059==    at 0x3625CDAE60: __write_nocancel (in /lib64/libc-2.12.so)
==1059==    by 0x3625C71582: _IO_file_write@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)
==1059==    by 0x3625C72B34: _IO_do_write@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)  
==1059==    by 0x3625C711FC: _IO_file_xsputn@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)
==1059==    by 0x3625C6760C: fwrite (in /lib64/libc-2.12.so)                      
==1059==    by 0x4012F2: write_matrix (io_file.c:78)                              
==1059==    by 0x4048F5: initialize (main.c:436)                                  
==1059==    by 0x404C9F: main (main.c:512)                                        
==1059==  Address 0x4c0e00a is not stack'd, malloc'd or (recently) free'd         
==1059==  Uninitialised value was created by a stack allocation                   
==1059==    at 0x4029D9: initialize (main.c:37)                                   
==1059==                                                                          
==1059== Syscall param write(buf) points to uninitialised byte(s)                 
==1059==    at 0x3625CDAE60: __write_nocancel (in /lib64/libc-2.12.so)            
==1059==    by 0x3625C71582: _IO_file_write@@GLIBC_2.2.5 (in /lib64/libc-2.12.so) 
==1059==    by 0x3625C71449: _IO_file_xsputn@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)
==1059==    by 0x3625C6760C: fwrite (in /lib64/libc-2.12.so)                      
==1059==    by 0x4012F2: write_matrix (io_file.c:78)                              
==1059==    by 0x4048F5: initialize (main.c:436)                                  
==1059==    by 0x404C9F: main (main.c:512)                                        
==1059==  Address 0x4c3a13a is 4,090 bytes inside a block of size 800,000 alloc'd 
==1059==    at 0x4A05FDE: malloc (vg_replace_malloc.c:236)                        
==1059==    by 0x40326C: initialize (main.c:201)                                  
==1059==    by 0x404C9F: main (main.c:512)                                        
==1059==  Uninitialised value was created by a heap allocation                    
==1059==    at 0x4A05FDE: malloc (vg_replace_malloc.c:236)                        
==1059==    by 0x4036AF: initialize (main.c:238)                                  
==1059==    by 0x404C9F: main (main.c:512)                                        
==1059==                                                                          
==1059== Syscall param write(buf) points to uninitialised byte(s)                 
==1059==    at 0x3625CDAE60: __write_nocancel (in /lib64/libc-2.12.so)            
==1059==    by 0x3625C71582: _IO_file_write@@GLIBC_2.2.5 (in /lib64/libc-2.12.so) 
==1059==    by 0x3625C72B34: _IO_do_write@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)   
==1059==    by 0x3625C7230F: _IO_file_close_it@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)
==1059==    by 0x3625C660B7: fclose@@GLIBC_2.2.5 (in /lib64/libc-2.12.so)
==1059==    by 0x404940: initialize (main.c:439)
==1059==    by 0x404C9F: main (main.c:512)
==1059==  Address 0x4c0e00a is not stack'd, malloc'd or (recently) free'd
==1059==  Uninitialised value was created by a heap allocation
==1059==    at 0x4A05FDE: malloc (vg_replace_malloc.c:236)
==1059==    by 0x401EDB: load_files1 (io_file.c:230)
==1059==    by 0x403167: initialize (main.c:187)
==1059==    by 0x404C9F: main (main.c:512)
==1059==

gcc  -Wall -g -o staudago -I/usr/include/gls -L/usr/lib64 *.c  -lgsl -lgslcblas -lm

*/

#include "io_file.h"

//Count the number of rows and columns of a binary file
void uliRows_uliColumns(FILE *data, ULI *uliRows, ULI *uliColumns){
  long c;
  ULI t=0,n=0;

  while ((c=getc(data)) !=EOF){
  	  if (c == '\t')   t++;
  	  if (c == '\n')  break;
  }
  while ((c=getc(data)) !=EOF)
  {
  	  if (c == '\n')
  		  n++;
  }
  rewind(data);
  *uliRows = n + 1;
  *uliColumns = t + 1;
}



//Allocate and read an unsigned table from a file
unsigned **load_table(ULI *uliRows_Ptr, ULI *uliColumns_Ptr, FILE *fPtr) {
  unsigned **m, j=0, i=0;

  uliRows_uliColumns(fPtr, uliRows_Ptr, uliColumns_Ptr);

  m=malloc(*uliRows_Ptr*sizeof(void *) );
  
  for(i=0;i<*uliRows_Ptr;i++) m[i] = malloc(*uliColumns_Ptr*sizeof(unsigned) );
  
  for(j=0; j < *uliRows_Ptr; j++){
  	  for(i=0; i< *uliColumns_Ptr - 1; i++) fscanf(fPtr,"%i\t", &m[j][i]);
  	  i=*uliColumns_Ptr -1;
  	  fscanf(fPtr,"%i\n", &m[j][i]);
  }

  return m;
}



/*Allocate and read an integer table from a file*/
int **iLoad_table(FILE *fPtr, int id, struct globale *G) {
  int **m;
  unsigned  j=0, i=0;

  m= malloc(G->uliRows[id]*sizeof(void *) );
  for(i=0;i<G->uliRows[id];i++) m[i] = malloc(G->uliColumns[id]*sizeof(int) );

  for(j=0; j < G->uliRows[id]; j++){
	for(i=0; i< G->uliColumns[id] - 1; i++){
		fscanf(fPtr,"%d\t", &m[j][i]);
	}
	i=G->uliColumns[id]-1;
	fscanf(fPtr,"%d\n", &m[j][i]);
  }
  return m;
}


void add_vector_in_matrix(LD *matrix, LD *M_in, ULI col, LD time, int membid, struct globale *G) {
  matrix[(col+1)*G->uliRow_counter[membid]] = time;
  memcpy(&(matrix[( (col+1)*G->uliRow_counter[membid] )+1]), M_in, ( sizeof(LD)* col ));
  G->uliRow_counter[membid]++;
}



//Write the output matrix in the output file
void write_matrix(int membid, struct globale *G) { 
  struct timeval tstart, tstop;
  float elapsed;
	
  printf("%d: Scrivo %ld %ld bytes\n",membid,G->uliRow_counter[membid],sizeof(LD)*(G->uliColumns[membid]+1)*G->uliRow_counter[membid]);
	
  gettimeofday(&tstart,NULL);
 
  fwrite(G->ldM_matrix[membid], sizeof(LD)*(G->uliColumns[membid]+1),G->uliRow_counter[membid],G->m_fPtr[membid]);
  
  gettimeofday(&tstop,NULL);

  G->uliRow_counter[membid] = 0;
  
  elapsed = (tstop.tv_sec-tstart.tv_sec)*1000000;
  elapsed += (tstop.tv_usec-tstart.tv_usec);
  G->iotime += elapsed/1000000;	  
}



//Read the input files
int load_files1(int idf, int base, struct globale *G, unsigned*** LS, unsigned*** RS, int*** iTgt_matrix) {
  FILE *fPtr;
  char app[20], stringa[60];
  int i,j,id;
  ULI uliTmp, uliRowss, uliColumnss;
  int iTgtcheck = 0;
  int indexes;

  sprintf(app,"%d",idf);
  
  id = idf - base;
	
  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/left_side_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/left_side", "r")) == NULL){
		fprintf (G->log_fPtr[id], "\nCannot open the left_side\n");
		exit(0);
	}
		
  LS[id] = load_table(&uliRowss, &uliColumnss, fPtr);
  
  G->side_index[id+1]  = G->side_index[id] + (uliRowss * uliColumnss);
  G->side_index2[id+1] = G->side_index2[id] + uliRowss;
  
  //printf("%d: %d %d\n",id+1,G->side_index[id+1],G->side_index2[id+1]);
  ///printf("%3d: row %3d col %3d index %10d (%5d and %5d)\n",id,uliRowss , uliColumnss,G->side_index[id],uliRowss * uliColumnss,G->side_index[id-1]);
  
  fclose (fPtr);
  
  G->uliRows[id] = uliRowss;
  G->uliColumns[id] = uliColumnss;
  
  //Print the left hand side into the log file
  fprintf(G->log_fPtr[id],"loading left_side\n");
  for(i=0;i<uliRowss;i++){
	for(j=0;j<uliColumnss;j++) fprintf(G->log_fPtr[id],"%u\t",LS[id][i][j]);
	fprintf(G->log_fPtr[id],"\n");
  }

  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/right_side_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/right_side", "r")) == NULL){
		fprintf (G->log_fPtr[id],"\nCannot open the right_side\n");
		exit(0);
	}
	
  RS[id] = load_table(&uliRowss, &uliColumnss, fPtr);
  fclose (fPtr);
  
  if(uliRowss != G->uliRows[id] || G->uliColumns[id] != uliColumnss) 
  	printf("DIVERSI %lu %lu %lu %lu\n", G->uliRows[id],uliRowss, G->uliColumns[id],uliColumnss);

  //Print the right hand side into the log file
  fprintf(G->log_fPtr[id],"loading right_side\n");
  for(i=0;i<uliRowss;i++){
	for(j=0;j<uliColumnss;j++) fprintf(G->log_fPtr[id],"%u\t",RS[id][i][j]);
	fprintf(G->log_fPtr[id],"\n");
  }

/*
  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/c_vector_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/c_vector", "r")) == NULL){
		fprintf (G->log_fPtr[id],"\nCannot open the c_vector\n");
		exit(0);
	}
	
  G->ldC_vector[id] = malloc(sizeof(LD) * uliRowss);
  
  //Read and print the constants vector in the log file
  fprintf(G->log_fPtr[id],"c vector\n");
  for(j=0; j < uliRowss; j++) {
  	fscanf(fPtr,         "%LG\n", &G->ldC_vector[id][j]);
	fprintf(G->log_fPtr[id],"%LG\n",  G->ldC_vector[id][j]);	
  }
  fclose (fPtr);
  fprintf(G->log_fPtr[id],"\n");
*/
	
  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/tgt_vector_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	  if((fPtr = fopen ("input/tgt_vector", "r")) == NULL)
  		  fprintf (G->log_fPtr[id],"\nCannot open the tgt_vector\n");
  if(fPtr != NULL) {
	/*Read and print the targets vector in the log file*/
  	G->iTgt_vector[id] = malloc(uliRowss * sizeof(int));
  	iTgtcheck = 1;
  	fprintf(G->log_fPtr[id],"loading tgt vector\n");
  	for(j=0; j < uliRowss; j++) { 
		fscanf(fPtr,"%d\n", &G->iTgt_vector[id][j]); 
		fprintf(G->log_fPtr[id],"%d\n",G->iTgt_vector[id][j]); 
	}
  	fclose (fPtr);
  }
  else {
	strcpy(stringa,"input/tgt_matrix_");
  	strcat(stringa,app);
  	if ((fPtr = fopen (stringa, "r")) == NULL)
		if((fPtr = fopen ("input/tgt_matrix", "r")) == NULL){
  	        	fprintf (G->log_fPtr[id],"\nCannot open the tgt_matrix\n");
  	        	exit(0);
  	        }
  	if(fPtr != NULL) {
		/*Read and print the targets matrix in the log file*/
		iTgt_matrix[id] = iLoad_table(fPtr,id,G);
  	    	fclose (fPtr);
		
  	    	fprintf(G->log_fPtr[id],"loading tgt_matrix\n");
  	    	for(i=0;i<uliRowss;i++){
			for(j=0;j<uliColumnss;j++)
				fprintf(G->log_fPtr[id],"%d\t",iTgt_matrix[id][i][j]);
	        		fprintf(G->log_fPtr[id],"\n");
  	        	}
	}
  	fprintf(G->log_fPtr[id],"\n");
  }
  fprintf(G->log_fPtr[id],"\n");

  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/M_0_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	  if((fPtr = fopen ("input/M_0", "r")) == NULL){
		  fprintf (G->log_fPtr[id],"\nCannot open the M_0\n");
		  exit(0);
	  }
	
  if(fPtr != NULL) {
	  //Read and print the multiset vector in the log file
	uliRows_uliColumns(fPtr, &uliTmp, &uliColumnss);
	  
	if(G->uliColumns[id] != uliColumnss) printf("DIVERSI %lu %lu\n", G->uliColumns[id],uliColumnss);

	  G->M[id]    = malloc( sizeof(LD) * uliColumnss);
	  G->Mbkp[id] = malloc( sizeof(LD) * uliColumnss);
	  fprintf(G->log_fPtr[id],"M_0\t");
	  
	  for(i=0; i < uliColumnss - 1; i++) {
		  fscanf(fPtr,     "%LG\t", &G->M[id][i]);
		  fprintf(G->log_fPtr[id],"%LG\t", G->M[id][i]);
	  }
	  fscanf(fPtr,      "%LG", &G->M[id][i]);
	  fprintf(G->log_fPtr[id],"%LG\n", G->M[id][i]);

	  fclose (fPtr);
  }

/*
  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/M_feed_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	  if((fPtr = fopen ("input/M_feed", "r")) == NULL){
  		  fprintf (G->log_fPtr[id],"\nCannot open the M_feed\n");
  		  exit(0);
  	  }
  if(fPtr != NULL) {
  	  //Read and print the feeding multiset vector in the log file
  	  G->ldM_feed[id] = malloc(sizeof(LD) * uliColumnss);
  	  fprintf(G->log_fPtr[id],"M_feed\t");
  	  for(i=0; i < uliColumnss - 1; i++) {
  		  fscanf(fPtr,    "%LG\t", &G->ldM_feed[id][i]);
  		  fprintf(G->log_fPtr[id],"%LG\t", G->ldM_feed[id][i]);
  	  }
  	  fscanf(fPtr,      "%LG", &G->ldM_feed[id][i]);
  	  fprintf(G->log_fPtr[id],"%LG\n", G->ldM_feed[id][i]);
  	  fclose (fPtr);
  }
*/
  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/indexes_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	  if((fPtr = fopen ("input/indexes", "r")) == NULL){
  		  fprintf (G->log_fPtr[id],"\nCannot open indexes\n");
  		  exit(0);
  	  }
  if(fPtr != NULL) {
  	  /*Read and print the indexes in the log file*/
  	  fprintf(G->log_fPtr[id],"\n indexes \t");
  	  uliRows_uliColumns(fPtr, &uliTmp, &G->uliIndexCol);
/*
  	  if(uliTmp != uliRows[id] || uliColumns[id] != uliIndexCol) 
  		  printf("DIVERSI %d %d %d %d\n", uliRows[id],uliTmp, uliColumns[id],uliIndexCol);
*/
  	  
  	  fprintf(G->log_fPtr[id],"%lu\n",G->uliIndexCol);
  	  for (i=0; i<G->uliIndexCol; i++){
  		  fscanf(fPtr,     "%d\t", &indexes);
  		  fprintf(G->log_fPtr[id],"%d\t", indexes);
  	  }
  	  fprintf(G->log_fPtr[id],"\n");
  	  fclose (fPtr);
  }

  // The every input file is read only once so moved outside the main loop 

/*
  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/sizes_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	  if((fPtr = fopen ("input/sizes", "r")) == NULL){
  		  fprintf (G->log_fPtr[id],"\nCannot open sizes\n");
  		  exit(0);
  	  }
  if(fPtr != NULL) {
  	  G->ldMolSize[id] = malloc (uliColumnss * sizeof(LD));
  	  fprintf(G->log_fPtr[id],"\n MolSize \n");
  	  for(i=0; i < uliColumnss - 1; i++) {
  		  fscanf(fPtr,     "%LG\t", &G->ldMolSize[id][i]);
  		  fprintf(G->log_fPtr[id],"%LG\t",  G->ldMolSize[id][i]);
  	  }
  	  fscanf(fPtr,     "%LG\n", &G->ldMolSize[id][i]);
  	  fprintf(G->log_fPtr[id],"%LG\n",  G->ldMolSize[id][i]);

  	  fprintf(G->log_fPtr[id],"\n MembSize \n");	//Size of the membrane
  	  fscanf(fPtr,     "%LG", &G->ldMembSize[id]);
  	  fprintf(G->log_fPtr[id],"%LG\n",  G->ldMembSize[id]);

  	  fclose (fPtr);
  }
*/  
  return(iTgtcheck);
}

void load_files2(int idf, int base, struct globale *G) {
  FILE *fPtr;
  char app[20], stringa[60];
  int i,j,id;
  ULI uliTmp, uliRowss, uliColumnss;
  int indexes;

  sprintf(app,"%d",idf);
  
  id = idf - base;
  uliRowss = G->uliRows[id];
  uliColumnss = G->uliColumns[id];


  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/c_vector_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/c_vector", "r")) == NULL){
		fprintf (G->log_fPtr[id],"\nCannot open the c_vector\n");
		exit(0);
	}
	
  
  //Read and print the constants vector in the log file
  fprintf(G->log_fPtr[id],"c vector\n");
  indexes = G->side_index2[id];
  
  for(j=0; j < uliRowss; j++) {
  	fscanf(fPtr,            "%LG\n", &G->ldC_vector[indexes+j]);
	fprintf(G->log_fPtr[id],"%LG\n",  G->ldC_vector[indexes+j]);	
  }
  fprintf(G->log_fPtr[id],"\n");
  fclose (fPtr);
	
/*
  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/M_0_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	  if((fPtr = fopen ("input/M_0", "r")) == NULL){
		  fprintf (G->log_fPtr[id],"\nCannot open the M_0\n");
		  exit(0);
	  }
	
  if(fPtr != NULL) {
	  //Read and print the multiset vector in the log file
	uliRows_uliColumns(fPtr, &uliTmp, &uliColumnss);
	  
	if(G->uliColumns[id] != uliColumnss) 
  		printf("DIVERSI %lu %lu\n", G->uliColumns[id],uliColumnss);

	  G->M[id]    = malloc( sizeof(LD) * uliColumnss);
	  G->Mbkp[id] = malloc( sizeof(LD) * uliColumnss);
	  fprintf(G->log_fPtr[id],"M_0\t");
	  
	  for(i=0; i < uliColumnss - 1; i++) {
		  fscanf(fPtr,     "%LG\t", &G->M[id][i]);
		  fprintf(G->log_fPtr[id],"%LG\t", G->M[id][i]);
	  }
	  fscanf(fPtr,      "%LG", &G->M[id][i]);
	  fprintf(G->log_fPtr[id],"%LG\n", G->M[id][i]);

	  fclose (fPtr);
  }

*/



  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/M_feed_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	if((fPtr = fopen ("input/M_feed", "r")) == NULL){
  		fprintf (G->log_fPtr[id],"\nCannot open the M_feed\n");
  		exit(0);
  	}
  if(fPtr != NULL) {
  	/* Read and print the feeding multiset vector in the log file*/
  	fprintf(G->log_fPtr[id],"M_feed\t");
	indexes = id*uliColumnss;
	
  	for(i=0; i < uliColumnss; i++) {
  	        fscanf(fPtr,	"%LG\t", &G->ldM_feed[indexes+i]);
  	        fprintf(G->log_fPtr[id],"%LG\t", G->ldM_feed[indexes+i]);
  	}
  	fprintf(G->log_fPtr[id],"\n");
	fclose (fPtr);
  }



  //Read the input file with the process id or the common input file
  strcpy(stringa,"input/sizes_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	if((fPtr = fopen ("input/sizes", "r")) == NULL){
  		fprintf (G->log_fPtr[id],"\nCannot open sizes\n");
  		exit(0);
  	}
  if(fPtr != NULL) {
  	fprintf(G->log_fPtr[id],"\n MolSize \n");
	indexes = id*uliColumnss;
		
  	for(i=0; i < uliColumnss; i++) {
		fscanf(fPtr,	 "%LG\t", &G->ldMolSize[indexes+i]);
  	        fprintf(G->log_fPtr[id],"%LG\t",  G->ldMolSize[indexes+i]);
  	}

  	fprintf(G->log_fPtr[id],"\n\n MembSize \n");    //Size of the membrane
  	fscanf(fPtr,	 "%LG", &G->ldMembSize[id]);
  	fprintf(G->log_fPtr[id],"%LG\n",  G->ldMembSize[id]);

  	fclose (fPtr);
  }
}

