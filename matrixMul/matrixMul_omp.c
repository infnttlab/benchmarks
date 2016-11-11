#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

int help_func(){
        printf("\nUsage: ./a.out <ROW_A> <COL_A> <COL_B> <THREADS> <DEBUG>\n");
        printf("Where: MATRIX(ROW, COL)  and  <COL_A> == <ROW_B>\n\n");
        printf("Default: A = (512,512) B = (512,512); THREADS = 2; DEBUG = 0\n\n");

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
		int debug = 0;
                int num_thread = 2;
                int row_a = 512, col_a = 512;
                int row_b = 512, col_b = 512;
                if(argc >= 2){
                        // change ROW_A
                        row_a = atoi(argv[1]);

                        if(argc >= 3){
                                // change ROW_A COL_A and ROW_B, where COL_A = ROW_B
                                col_a = row_b = atoi(argv[2]);
                                if(argc >= 4){
                                        col_b = atoi(argv[3]);
                                        if(argc >= 5){
                                                num_thread = atoi(argv[4]);
						if(argc == 6)
							debug = atoi(argv[5]);
                                        }
                                }
                        }
                }
                printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n",
			row_a, col_a, row_b, col_b, row_a, col_b);

                struct timeval tstart, tstop;
                double elapsed = 0.0;
                double start_time, run_time;

                omp_set_num_threads(num_thread);
                start_time = omp_get_wtime();

                gettimeofday(&tstart,NULL);

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

		   #pragma omp for schedule(static) 
                   for(i=0; i<row_a; i++){
                        for(j=0; j<col_b; j++){
				matrix_c[i*col_b+j] = 0.f;
                                for(k=0; k<col_a; k++){
                                        matrix_c[i*col_b+j] += matrix_a[i*col_a+k]*matrix_b[k*col_b+j];
                                }
                        }
                   }
		}

		if(debug){
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
                run_time = omp_get_wtime() - start_time;

                printf("\nTerminated.\n");
                elapsed = (tstop.tv_sec - tstart.tv_sec) + ((tstop.tv_usec - tstart.tv_usec)/1000000.0);
                printf("Data processing with %d threads in %.10f s using \"gettimeofday\" and  %.10f s using \"OMP timer\".\n\n", num_thread,elapsed,run_time);
        }
        return val_returned;
}
