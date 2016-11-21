#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>
#include "helper_string.h"


void help_func(){
        printf("\nUsage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB  | matrix(row,col), ColumnsA = RowsB divisibile per num.processi\n");
        printf("         -w=WarmUpData\n");
        printf("         -v=Verbose\n\n");
        printf("Default: A = (512,512) B = (512,512); WARMUP = 0; VERBOSE = 0\n\n");
}


void fill_matrix(int row, int col, float *matrix){
	int i, j;
	for (i=0; i<row; i++){
		for (j=0; j<col; j++){
			matrix[i*col+j] = 0.1f;
		}
	}
}

void printMatrix(int row, int col, float *matrix){
	int i, j = 0;
	for (i=0; i<row; i++) {
		for (j=0; j<col; j++)
			printf("%f ", matrix[i*col+j]);
		printf("\n");
  	}
}


int main(int argc, char *argv[]){

	int myrank, nProc, i, j, k;
	//MPI_Status status;
  
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	/* who am i */
	MPI_Comm_size(MPI_COMM_WORLD, &nProc); /* number of processors */

	//int val_returned = 0;
        if(argc == 2 && ( (strcmp(argv[1], "--help")==0) || (strcmp(argv[1], "-h")==0) )){
                if(myrank==0)
			help_func();
		MPI_Finalize();
		return 0;
        }
	else{
		int row_a = 512, col_a = 512, col_b = 512;
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

		if (row_a%nProc != 0){
                       	if (myrank == 0) printf("\nERROR: N.Row A(%d) is not divisible by n.Proc(%d).\n\n",row_a,nProc);
                }
		else{
			if(myrank == 0)
                         printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n\n",
                                row_a, col_a, col_a, col_b, row_a, col_b);

			double start_tot, time_tot, s_mm, time_wrm, time_mm, delta_wrm, delta_mm, elapsed;
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
  			float *rowsCxProc;  
			int rows4proc = row_a/nProc;

	  		rowsAxProc = (float*)malloc(rows4proc*col_a*sizeof(float));
			rowsCxProc = (float*)malloc(rows4proc*col_b*sizeof(float));
	
  			MPI_Bcast (matrix_b, col_a*col_b, MPI_FLOAT, 0, MPI_COMM_WORLD);
  			MPI_Scatter (matrix_a, rows4proc*col_a, MPI_FLOAT, rowsAxProc, rows4proc*col_a, MPI_FLOAT, 0, MPI_COMM_WORLD);

			if(perf){
				s_mm = MPI_Wtime();
				//Performs warmup operations:
				for(i=0; i<rows4proc; i++){
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
			for(i=0; i<rows4proc; i++){
				for(j=0; j<col_b; j++){
					rowsCxProc[i*col_b+j] = 0.f;
					for(k=0; k<col_a; k++){
						rowsCxProc[i*col_b+j] += rowsAxProc[i*col_a+k]*matrix_b[k*col_b+j];
					}
				}
			}
			time_mm = MPI_Wtime() - s_mm;
	
			free(rowsAxProc);

  			if(myrank == 0){
				matrix_c = (float*)malloc(row_a*col_b*sizeof(float));
  			}
			MPI_Gather (rowsCxProc, rows4proc*col_b, MPI_FLOAT, matrix_c, rows4proc*col_b, MPI_FLOAT, 0, MPI_COMM_WORLD);

  			free(rowsCxProc); 

  			if (myrank==0) {
				if(debug == 2){
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
				elapsed = sum_tot/(double)nProc;
				delta_mm = sum_mm/(double)nProc;
	
				double tTot;

				if(perf){
					delta_wrm = sum_wrm/(double)nProc;
                			tTot = elapsed-delta_wrm;
				}
				else
					tTot = elapsed;

				double flops = 2.0*(double)row_a*(double)col_a*(double)col_b;
				double gigaf = (flops * 1.0e-9f) / delta_mm;
			
				if(debug)
               		        	printf("\nThreads: %d,  Flop: %.0f,  GFlop: %f GFlop/s,  Time_mtxMul: %f s,  Time_tot: %f s\n\n",
                       		        	nProc, flops, gigaf, delta_mm, tTot);
	               		else
                        		printf("\n%d %.0f %f %f %f\n\n",
               	                		nProc, flops, gigaf, delta_mm, tTot);
			}
		}
	}

  	MPI_Finalize();
  	return 0;
}


