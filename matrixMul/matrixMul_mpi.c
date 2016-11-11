#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>


int help_func(){
        printf("\nUsage: mpirun -np <P> ./a.out <ROW_A> <COL_A> <COL_B> <DEBUG>\n");
        printf("Where: MATRIX(ROW, COL)  and  <COL_A> == <ROW_B>; <P> : parallel process\n");
	printf("       <COL_A> must be divisible by <P>\n\n");
        printf("Default: A = (512,512) B = (512,512); DEBUG = 0\n\n");

        return 0;
}


void fill_matrix(int row, int col, float *matrix){
	int i, j;
	for (i=0; i<row; i++){
		for (j=0; j<col; j++)
			matrix[i*col+j] = 0.1f;
	}
}


void printMatrix(int row, int col, float *matrix){
	int i, j = 0;
	for (i=0; i<row; i++) {
		for (j=0; j<col; j++)
			printf("%f ", matrix[i*row+j]);
		printf("\n");
  	}
}


int main(int argc, char *argv[]){

	int myrank, nProc, i, j, k;
	MPI_Status status;
  
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	/* who am i */
	MPI_Comm_size(MPI_COMM_WORLD, &nProc); /* number of processors */

	int val_returned = 0;
        if(argc == 2 && ( (strcmp(argv[1], "--help")==0) || (strcmp(argv[1], "-h")==0) )){
                if(myrank==0)
			val_returned = help_func();
		MPI_Finalize();
		return 0;
        }
	else{
		int row_a = 512, col_a = 512;
                int row_b = 512, col_b = 512;
                int debug = 0;
                if(argc >= 2){
                        // change ROW_A
                        row_a = atoi(argv[1]);

                        if(argc >= 3){
                                // change ROW_A COL_A and ROW_B, where COL_A = ROW_B
                                col_a = row_b = atoi(argv[2]);
                                if(argc >= 4){
                                        col_b = atoi(argv[3]);
                                        if(argc == 5)
                                                debug = atoi(argv[4]);
                                }
                        }
                }

		if (row_a%nProc != 0) {
    			if (myrank == 0) printf("Width A not divisible by number of processors\n");
    			MPI_Finalize();
    			return 0;
  		}
		
		if(myrank == 0)
			 printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n",
				row_a, col_a, row_b, col_b, row_a, col_b);

		struct timeval tstart, tstop;
		double elapsed;
		if(myrank == 0){
			elapsed = 0.0;
			gettimeofday(&tstart,NULL);
		}

  		float *matrix_a;
  		float *matrix_b;
  		float *matrix_c;

  		matrix_b = (float*)malloc(col_a*col_b*sizeof(float));

  		if (myrank==0) {
    			matrix_a = (float*)malloc(row_a*col_a*sizeof(float));
    			fill_matrix(row_a, col_a, matrix_a);
    			fill_matrix(col_a, col_b, matrix_b);
  		}

  		float *rowsAxProc;
  		float *rowsCxProc;

  		rowsAxProc = (float*)malloc(nProc*col_a*sizeof(float));
  		rowsCxProc = (float*)malloc(nProc*col_b*sizeof(float));

  		MPI_Bcast (matrix_b, col_a*col_b, MPI_FLOAT, 0, MPI_COMM_WORLD);
  		MPI_Scatter (matrix_a, nProc*col_a, MPI_FLOAT, rowsAxProc, nProc*col_a, MPI_FLOAT, 0, MPI_COMM_WORLD);

		for(i=0; i<nProc; i++){
			for(j=0; j<col_b; j++){
				rowsCxProc[i*col_b+j] = 0.f;
				for(k=0; k<col_a; k++){
					if(debug)
						printf("** rank %d ** rowsCxProc[%d]: %f; rowsAxProc[%d]*B[%d]: %f*%f\n",
							myrank, i*col_b+j, rowsCxProc[i*col_b+j], i*col_a+k, k*col_b+j,
							rowsAxProc[i*col_a+k], matrix_b[k*col_b+j]);

					rowsCxProc[i*col_b+j] += rowsAxProc[i*col_a+k]*matrix_b[k*col_b+j];
				}
			}
  		}

 		free(rowsAxProc);

  		if(myrank == 0){
			matrix_c = (float*)malloc(row_a*col_b*sizeof(float));
  		}

  		MPI_Gather (rowsCxProc, nProc*col_b, MPI_FLOAT, matrix_c, nProc*col_b, MPI_FLOAT, 0, MPI_COMM_WORLD);

  		free(rowsCxProc); 

  		if (myrank==0) {
			if(debug){
    				printf("\n## Matrix A:\n");
    				printMatrix(row_a, col_a, matrix_a);
				printf("\n## Matrix B:\n");
    				printMatrix(col_a, col_b, matrix_b);
				printf("\n## Matrix C:\n");
    				printMatrix(row_a, col_b, matrix_c);
			}
    			free(matrix_a); free(matrix_c);
  		}

  		free(matrix_b);
		
		if(myrank == 0){
			gettimeofday(&tstop,NULL);
                	printf("\nTerminated.\n");
                	elapsed = (tstop.tv_sec - tstart.tv_sec) + ((tstop.tv_usec - tstart.tv_usec)/1000000.0);
                	printf("Data processing in %f s.\n\n", elapsed);
		}

	}

  	MPI_Finalize();
  	return 0;
}


