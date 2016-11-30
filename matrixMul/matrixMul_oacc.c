#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "helper_string.h"

void printMatrix(int row, int col, float *matrix){
        int i, j;
        for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                      printf("%f ", matrix[i*col+j]);
                }
                printf("\n");
        }
}

int help_func(){
        printf("\nUsage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB  | matrix(row,col), ColumnsA = RowsB\n");
        printf("         -w=WarmUpData\n");
        printf("         -v=Verbose\n\n");
        printf("Default: A = (512,512) B = (512,512); WARMUP = 0; VERBOSE = 0\n\n");

        return 0;
}

void fill_matrix(int row, int col, float *matrix){
        int i, j;
         #pragma acc kernels
        for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                        matrix[i*col+j] =  0.1f;
                }
        }
}

int main(int argc, char **argv){
        int val_returned = 0;
        if(argc == 2 && ( (strcmp(argv[1], "--help")==0) || (strcmp(argv[1], "-h")==0) )){
                val_returned = help_func();
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

                printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n",
			row_a, col_a, col_a, col_b, row_a, col_b);

                struct timeval tstart, tstop, mm_s, mm_e;
                double elapsed = 0.0;
		double t_mm = 0.0;
		double t_w = 0.0;

                gettimeofday(&tstart,NULL);

                float *matrix_a;
                float *matrix_b;
                float *matrix_c;

		matrix_a = (float*)malloc(row_a*col_a*sizeof(float));
		matrix_b = (float*)malloc(col_a*col_b*sizeof(float));
		matrix_c = (float*)malloc(row_a*col_b*sizeof(float));

                //fill matrix A and B with random float, range 0-1

                fill_matrix(row_a, col_a, matrix_a);
                fill_matrix(col_a, col_b, matrix_b);
                
		//matrix multiplication
		int i, j, k;

		if(perf){
			//Performs warmup operations
			gettimeofday(&mm_s,NULL);
			#pragma acc kernels copyin(matrix_a,matrix_b) copy(matrix_c)
                	for(i=0; i<row_a; i++){
                        	for(j=0; j<col_b; j++){
                                	matrix_c[i*col_b+j] = 0.f;
                                	for(k=0; k<col_a; k++){
                                        	matrix_c[i*col_b+j] += matrix_a[i*col_a+k]*matrix_b[k*col_b+j];
                                	}
                        	}
                	}
			gettimeofday(&mm_e,NULL);
			t_w = (mm_e.tv_sec - mm_s.tv_sec) + ((mm_e.tv_usec - mm_s.tv_usec)/1000000.0);
		}

		gettimeofday(&mm_s,NULL);
                #pragma acc kernels copyin(matrix_a,matrix_b) copy(matrix_c)
                for(i=0; i<row_a; i++){
                        for(j=0; j<col_b; j++){
				matrix_c[i*col_b+j] = 0.f;
                                for(k=0; k<col_a; k++){
                                        matrix_c[i*col_b+j] += matrix_a[i*col_a+k]*matrix_b[k*col_b+j];
                                }
                        }
                }
		gettimeofday(&mm_e,NULL);
		t_mm = (mm_e.tv_sec - mm_s.tv_sec) + ((mm_e.tv_usec - mm_s.tv_usec)/1000000.0);

		if(debug == 2){
                        //print all matrix:
                        printf("\n## Matrix A:\n");
                        printMatrix(row_a, col_a, matrix_a);
                        printf("\n## Matrix B:\n");
                        printMatrix(col_a, col_b, matrix_b);
                        printf("\n## Matrix C:\n");
                        printMatrix(row_a, col_b, matrix_c);
                }

		free(matrix_a); free(matrix_b); free(matrix_c);

                gettimeofday(&tstop,NULL);

                elapsed = (tstop.tv_sec - tstart.tv_sec) + ((tstop.tv_usec - tstart.tv_usec)/1000000.0);

		double t_tot;
		if(perf)
                	t_tot = elapsed-t_w;
		else
			t_tot = elapsed;

		double flops = 2.0*(double)row_a*(double)col_a*(double)col_b;
		double giga = (flops*1.0e-9f)/t_mm;
		
	//	printf("\nProcessing: %f GFlop/s, Time: %f s, Flop: %.0f\n\n", giga, t_m, flops);

		if(debug)
                        printf("\nFlop: %.0f,  GFlop: %f GFlop/s,  Time_mtxMul: %f s,  Time_tot: %f s\n\n",
                                flops, giga, t_mm, t_tot);
                else
                        printf("\n%.0f %f %f %f\n\n",
                                flops, giga, t_mm, t_tot);
        }
        return val_returned;
}
