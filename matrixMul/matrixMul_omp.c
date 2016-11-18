#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "helper_string.h"

int help_func(){
        printf("\nUsage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB  | matrix(row,col), ColumnsA = RowsB\n");
	printf("         -p=Threads\n");
        printf("         -w=WarmUpData\n");
        printf("         -v=Verbose\n\n");
        printf("Default: A = (512,512) B = (512,512); THREADS = 2; WARMUP = 0; VERBOSE = 0\n\n");

        return 0;
}

void printMatrix(int row, int col, float *matrix){
        int i, j;
        for(i=0; i<row; i++){
                for(j=0; j<col; j++){
                      printf("%f ", matrix[i*col+j]);
                }
		printf("\n");
        }
}

int main(int argc, char **argv){
        int val_returned = 0;
        if(argc == 2 && ( (strcmp(argv[1], "--help")==0) || (strcmp(argv[1], "-h")==0) )){
                val_returned = help_func();
        }
        else{
		int debug = 0, perf = 0;
                int num_thread = 2;
                int row_a = 512, col_a = 512, col_b = 512;
		
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
		if (checkCmdLineFlag(argc, (const char **)argv, "p")){
                        num_thread = getCmdLineArgumentInt(argc, (const char **)argv, "p");
                }

                printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n\n",
			row_a, col_a, col_a, col_b, row_a, col_b);

                double start_time, run_time, start_mm, t_mm, t_warm;

		int th; int right = 0;
		for(th=2;th<=num_thread; th++){
			omp_set_num_threads(th);
                	start_time = omp_get_wtime();

	                float *matrix_a;
        	        float *matrix_b;
                	float *matrix_c;

			matrix_a = (float*)malloc(row_a*col_a*sizeof(float));
			matrix_b = (float*)malloc(col_a*col_b*sizeof(float));
			matrix_c = (float*)malloc(row_a*col_b*sizeof(float));

			int i,j,k;
	
			#pragma omp parallel shared(matrix_a, matrix_b, matrix_c) private(i,j,k)
			{
		   		#pragma omp for collapse(2) schedule(static)
		   		for(i=0; i<row_a; i++){
               				for(j=0; j<col_a; j++){
	                	 	       matrix_a[i*col_a+j] =  0.1f;
        		        	}
	        	   	}

			   	#pragma omp for collapse(2) schedule(static)
			   	for(i=0; i<col_a; i++){
                			for(j=0; j<col_b; j++){
                		        	matrix_b[i*col_b+j] =  0.1f;
		        	        }
			           }

				   if(perf){
		   			#pragma omp master
                        			 start_mm = omp_get_wtime();
	                   		#pragma omp for schedule(static)
        	           		for(i=0; i<row_a; i++){
                	        		for(j=0; j<col_b; j++){
                        	        		matrix_c[i*col_b+j] = 0.f;
	                        	        	for(k=0; k<col_a; k++){
        	                        	        	matrix_c[i*col_b+j] += matrix_a[i*col_a+k]*matrix_b[k*col_b+j];
	                	                	}
        	                		}
	        	           	}	
        	        	   	#pragma omp master
                	        		t_warm = omp_get_wtime() - start_mm;
		  	  	}

			   	#pragma omp master
					 start_mm = omp_get_wtime();
			   	#pragma omp for schedule(static) 
	                	for(i=0; i<row_a; i++){
        	                	for(j=0; j<col_b; j++){
						matrix_c[i*col_b+j] = 0.f;
                        	        	for(k=0; k<col_a; k++){
                                	        	matrix_c[i*col_b+j] += matrix_a[i*col_a+k]*matrix_b[k*col_b+j];
	                                	}
	        	                }
        	        	}
				#pragma omp master
					t_mm = omp_get_wtime() - start_mm;
			}

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
		
                	run_time = omp_get_wtime() - start_time;

			double time_tot;
			if(perf)
                		time_tot = run_time-t_warm;
			else
				time_tot = run_time;
		
			double flops = 2.0*(double)row_a*(double)col_a*(double)col_b;
			double gigaflop = (flops * 1.0e-9f) / t_mm;

//			printf("Performance: %f GFlop/s, Time: %f s, Flops: %.0f\n\n", gigaflop, t_mm, flops);
		
			if(debug)
                        	printf("Threads: %d,  Flop: %.0f,  GFlop: %f GFlop/s,  Time_mtxMul: %f s,  Time_tot: %f s\n",
                                	th, flops, gigaflop, t_mm, time_tot);
                	else
                        	printf("%d %.0f %f %f %f\n",
                                	th, flops, gigaflop, t_mm, time_tot);

		}
        }
	printf("\n");
        return val_returned;
}
