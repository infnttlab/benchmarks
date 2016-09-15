//g++ -o scrivioutput scrivioutput.cc 

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char *argv[]) {
  FILE* fptr1, * fptr2;
  int col, i,j;
  long double* riga, lcol;
  char file[100];

/*
 *   if ((fptr1 = fopen (argv[1],"w")) == NULL)	  {
 *   	printf ("\nCannot open %s\n",argv[1]);
 *   		exit(0);
 *   		  } else printf("Creo %s\n",argv[1]);
 *   		    
 *   		      lrow = 10.0;
 *   		        
 *   		          fwrite(&lrow,sizeof(long double),1,fptr1);
 *   		            fclose(fptr1);
 *   		            */
  if ((fptr1 = fopen (argv[1],"r")) == NULL)	  {
	printf ("\nCannot open %s\n",argv[1]);
	exit(0);
  } else printf("Apro %s\n",argv[1]);
  
  sprintf(file,"%s.txt",argv[1]);
  
  if ((fptr2 = fopen (file,"w")) == NULL)	  {
	printf ("\nCannot open %s\n",file);
	exit(0);
  }  else printf("Creo %s\n",file);
  
  j = 0;
  fread(&lcol,sizeof(long double),1,fptr1);
  
  col = (int) lcol;
  
  //printf("%d %LG\n",col,lcol);

  riga = (long double*) malloc(col*sizeof(long double));
  
  while(fread(riga,sizeof(long double),col,fptr1) == col) {
  	fprintf(fptr2,"%.15LG\t",riga[0]);
	for(i=1;i<col-1;i++) fprintf(fptr2,"%.10LG\t",riga[i]);
	fprintf(fptr2,"%.10LG\n",riga[col-1]);
	j++;
  }
  printf("%d colonne e %d righe\n",col,j);

  fclose(fptr1); fclose(fptr2);

}

