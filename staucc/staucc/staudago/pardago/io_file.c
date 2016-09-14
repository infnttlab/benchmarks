#include "io_file.h"

/*Count the number of rows and columns of a binary file*/
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



/*Write the state into the output file*/
void write_M(int membid, LD t, struct globale *G) {
  unsigned c=0;

  fprintf(G->m_fPtr[membid],"%.10LG\t",t);
  for(c=0;c<G->uliIndexCol - 1;c++) fprintf(G->m_fPtr[membid],"%LG\t",G->M[membid][G->indexes[membid][c]-1]);
  fprintf(G->m_fPtr[membid],"%LG\n",G->M[membid][G->indexes[membid][c]-1]);
}



/*Allocate and read an unsigned table from a file*/
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



/*Deallocate a table*/
unsigned **unload_table(ULI *uliRows_Ptr, ULI *uliColumns_Ptr, unsigned **table_Ptr) {
  unsigned i;

  for(i=0;i<*uliRows_Ptr;i++) free(table_Ptr[i]);
  free(table_Ptr);

  return 0;
}



/* Convert an integer into a char */
void itoa(int n, char s[]){
  int i,sign,j,c;

  if((sign = n)< 0) n = -n;
  i=0;
  do {
  	  s[i++] = n % 10 + '0';
  } while((n/=10)>0);

  if(sign < 0) s[i++] = '-';
  s[i]= '\0';

  for(i=0, j=strlen(s)-1;i<j; i++,j--){
  	  c=s[i];
  	  s[i]= s[j];
  	  s[j]=c;
  }
}


//Add a row to the output buffer
void add_vector_in_matrix(LD **matrix, LD *M_in, ULI col, LD time, int membid, struct globale *G) {
  matrix[G->uliRow_counter[membid]][0] = time;
  memcpy(&(matrix[G->uliRow_counter[membid]][1]), M_in, sizeof(LD)* col);
}



/*Write the output matrix in the output file*/
void write_matrix(FILE * fPtr, LD **matrix, ULI cols, int membid, struct globale *G) { 
   ULI c1=0,c2=0;

   for(c1=0; c1<G->uliRow_counter[membid];c1++){
   	   fprintf(fPtr,"%.10LG\t",matrix[c1][0]);
   	   for(c2=0;c2<cols-1;c2++) fprintf(fPtr,"%.10LG\t",matrix[c1][G->indexes[membid][c2]]);
   	   fprintf(fPtr,"%.10LG\n",matrix[c1][G->indexes[membid][c2]]);
   }
}



/*Read the input files*/
int load_files(int idf, int base, struct globale *G) {
  FILE *fPtr;
  char app[20], stringa[40];
  int i,j,id;
  ULI uliTmp, uliRowss, uliColumnss;
  int iTgtcheck;

  itoa(idf,app);

  id = idf - base;
	
  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/left_side_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/left_side", "r")) == NULL){
		fprintf (G->log_fPtr[id], "\nCannot open the left_side\n");
		exit(0);
	}
		
  G->unLeft_side[id] = load_table(&uliRowss, &uliColumnss, fPtr);
  fclose (fPtr);
  
  G->uliRows[id] = uliRowss;
  G->uliColumns[id] = uliColumnss;
  
  /*Print the left hand side into the log file*/
  fprintf(G->log_fPtr[id],"loading left_side\n");
  for(i=0;i<uliRowss;i++){
	for(j=0;j<uliColumnss;j++) fprintf(G->log_fPtr[id],"%u\t",G->unLeft_side[id][i][j]);
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
	
  G->unRight_side[id] = load_table(&uliRowss, &uliColumnss, fPtr);
  fclose (fPtr);
  
  if(uliRowss != G->uliRows[id] || G->uliColumns[id] != uliColumnss) 
  	printf("DIVERSI %d %d %d %d\n", G->uliRows[id],uliRowss, G->uliColumns[id],uliColumnss);

  /*Print the right hand side into the log file*/
  fprintf(G->log_fPtr[id],"loading right_side\n");
  for(i=0;i<uliRowss;i++){
	for(j=0;j<uliColumnss;j++) fprintf(G->log_fPtr[id],"%u\t",G->unRight_side[id][i][j]);
	fprintf(G->log_fPtr[id],"\n");
  }

  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/c_vector_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	if((fPtr = fopen ("input/c_vector", "r")) == NULL){
		fprintf (G->log_fPtr[id],"\nCannot open the c_vector\n");
		exit(0);
	}
	
  G->ldC_vector[id] = malloc(sizeof(LD) * uliRowss);
  
  /*Read and print the constants vector in the log file*/
  fprintf(G->log_fPtr[id],"c vector\n");
  for(j=0; j < uliRowss; j++) {
  	fscanf(fPtr,         "%LG\n", &G->ldC_vector[id][j]);
	fprintf(G->log_fPtr[id],"%LG\n",  G->ldC_vector[id][j]);	
  }
  fclose (fPtr);
  fprintf(G->log_fPtr[id],"\n");

	
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
		G->iTgt_matrix[id] = iLoad_table(fPtr,id,G);
  	    	fclose (fPtr);
		
  	    	fprintf(G->log_fPtr[id],"loading tgt_matrix\n");
  	    	for(i=0;i<uliRowss;i++){
			for(j=0;j<uliColumnss;j++)
				fprintf(G->log_fPtr[id],"%d\t",G->iTgt_matrix[id][i][j]);
	        		fprintf(G->log_fPtr[id],"\n");
  	        	}
	}
  	fprintf(G->log_fPtr[id],"\n");
  }
  fprintf(G->log_fPtr[id],"\n");

  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/M_0_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
	  if((fPtr = fopen ("input/M_0", "r")) == NULL){
		  fprintf (G->log_fPtr[id],"\nCannot open the M_0\n");
		  exit(0);
	  }
	
  if(fPtr != NULL) {
	  /*Read and print the multiset vector in the log file*/
	uliRows_uliColumns(fPtr, &uliTmp, &uliColumnss);
	  
	if(G->uliColumns[id] != uliColumnss) 
  		printf("DIVERSI %d %d\n", G->uliColumns[id],uliColumnss);

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

  /*Read the input file with the process id or the common input file*/
  strcpy(stringa,"input/M_feed_");
  strcat(stringa,app);
  if ((fPtr = fopen (stringa, "r")) == NULL)
  	  if((fPtr = fopen ("input/M_feed", "r")) == NULL){
  		  fprintf (G->log_fPtr[id],"\nCannot open the M_feed\n");
  		  exit(0);
  	  }
  if(fPtr != NULL) {
  	  /* Read and print the feeding multiset vector in the log file*/
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
  	  G->indexes[id] = malloc (G->uliIndexCol * sizeof(int));
  	  for (i=0; i<G->uliIndexCol; i++){
  		  fscanf(fPtr,     "%d\t", &G->indexes[id][i]);
  		  fprintf(G->log_fPtr[id],"%d\t",  G->indexes[id][i]);
  	  }
  	  fprintf(G->log_fPtr[id],"\n");
  	  fclose (fPtr);
  }

  // The every input file is read only once so moved outside the main loop 

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

  	  fprintf(G->log_fPtr[id],"\n MembSize \n");
  	  fscanf(fPtr,     "%LG", &G->ldMembSize[id]);
  	  fprintf(G->log_fPtr[id],"%LG",  G->ldMembSize[id]);

  	  fclose (fPtr);
  }
  
  return(iTgtcheck);
}
