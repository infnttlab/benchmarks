#include <stdio.h>
#include <stdlib.h>

int main( int argc, char *argv[]) {
  FILE* fptr1, * fptr2;
  int row, i,j;
  long double* riga, lrow;
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
  fread(&lrow,sizeof(long double),1,fptr1);
  
  row = (int) lrow;
  
  printf("%d %LG\n",row,lrow);

  riga = (long double*) malloc(row*sizeof(long double));
  
  while(fread(riga,sizeof(long double),row,fptr1) == row) {
	for(i=0;i<row-1;i++) fprintf(fptr2,"%.10LG\t",riga[i]);
	fprintf(fptr2,"%.10LG\n",riga[row-1]);
	j++;
  }
  printf("%d colonne e %d righe\n",row,j);

  fclose(fptr1); fclose(fptr2);

}

