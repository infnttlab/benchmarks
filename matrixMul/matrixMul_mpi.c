#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>
#include "helper_string.h"


int help_func(){
        printf("\nUsage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB  | matrix(row,col), ColumnsA = RowsB divisibile per num.processi\n");
        printf("         -w=WarmUpData\n");
        printf("         -v=Verbose\n\n");
        printf("Default: A = (512,512) B = (512,512); WARMUP = 0; VERBOSE = 0\n\n");

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
                int debug = 0, perf = 0;

		if (checkCmdLineFlag(argc, (const char **)argv, "rA")){
                        row_a = getCmdLineArgumentInt(argc, (const char **)argv, "rA");
                }
                if (checkCmdLineFlag(argc, (const char **)argv, "cA")){
                        col_a = getCmdLineArgumentInt(argc, (const char **)argv, "cA");
                }
                if (checkCmdLineFlag(argc, (const char **)argv, "cB")){
                        col_b = getCmdLineArgumentInt(argc, (const char **)argv, "cB");
                }
                if (checkCmdLineFlag(argc, (const char **)argv, "w")){
                        perf = getCmdLineArgumentInt(argc, (const char **)argv, "w");
                }
                if (checkCmdLineFlag(argc, (const char **)argv, "v")){
                        debug = getCmdLineArgumentInt(argc, (const char **)argv, "v");
                }

		if (row_a%nProc != 0) {
    			if (myrank == 0) printf("Width A not divisible by number of processors\n");
    			MPI_Finalize();
    			return 0;
  		}
		
		if(myrank == 0)
			 printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n",
				row_a, col_a, row_b, col_b, row_a, col_b);

		//struct timeval tstart, tstop;
		double start_tot, time_tot, s_mm, time_wrm, time_mm,  delta_wrm, delta_mm, elapsed;
		double sum_tot = 0.0, sum_wrm = 0.0, sum_mm = 0.0;
		if(myrank == 0){
			elapsed = delta_mm = delta_wrm = 0.0;
		}

		start_tot = MPI_Wtime();
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
  		float *rowsCxProc;\

  		rowsAxProc = (float*)malloc(nProc*col_a*sizeof(float));
  		rowsCxProc = (float*)malloc(nProc*col_b*sizeof(float));

  		MPI_Bcast (matrix_b, col_a*col_b, MPI_FLOAT, 0, MPI_COMM_WORLD);
  		MPI_Scatter (matrix_a, nProc*col_a, MPI_FLOAT, rowsAxProc, nProc*col_a, MPI_FLOAT, 0, MPI_COMM_WORLD);

		if(perf){
			s_mm = MPI_Wtime();
			//Performs warmup operations:
			for(i=0; i<nProc; i++){
                        	for(j=0; j<col_b; j++){
                                	rowsCxProc[i*col_b+j] = 0.f;
                                	for(k=0; k<col_a; k++){
                                        	rowsCxProc[i*col_b+j] += rowsAxProc[i*col_a+k]*matrix_b[k*col_b+j];
                                	}
                        	}
                	}
			time_wrm = MPI_Wtime() - s_mm;
		}

		s_mm = MPI_Wtime();
		//computing matrix multipliaction:
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
		time_mm = MPI_Wtime() - s_mm;

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
		time_tot = MPI_Wtime() - start_tot;

		MPI_Reduce(&time_tot, &sum_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
		MPI_Reduce(&time_mm, &sum_mm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if(perf)
			MPI_Reduce(&time_wrm, &sum_wrm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if(myrank == 0){
                	printf("\nTerminated.\n");

			elapsed = sum_tot/(double)nProc;
			delta_mm = sum_mm/(double)nProc;
			if(perf){
				delta_wrm = sum_wrm/(double)nProc;
                		printf("Data processing in %f s (time warmup: %f s).\n\n", elapsed-delta_wrm, delta_wrm);
			}
			else
				printf("Data processing in %f s.\n\n", elapsed);

			double flops = 2.0*(double)row_a*(double)col_a*(double)col_b;
			double gigaf = (flops * 1.0e-9f) / delta_mm;
			printf("Performance: %f GFlop/s, Time: %f ms, Flop: %.0f\n\n", gigaf, delta_mm, flops);
		}

	}

  	MPI_Finalize();
  	return 0;
}


