#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

int help_func(){
        printf("\nUsage: ./a.out <ROW_A> <COL_A> <COL_B> <THREADS>\n");
        printf("Where: MATRIX(ROW, COL)  and  <COL_A> == <ROW_B>\n\n");
        printf("Default: A = (512,512) B = (512,512); THREADS = 2\n\n");

        return 0;
}

void fill_matrix(int row, int col, float matrix[row][col]){
 //       float range_max = 1.f;
        int i, j;
        #pragma omp parallel for  shared(row,col,matrix) private (i,j)
        for(i=0; i<row; i++){
                for(j=0; j<col; j++){
   //                     matrix[i][j] = ((float)rand()/(float)(RAND_MAX))*range_max;
                        matrix[i][j] =  0.1f;
//                      printf("%f ", matrix[i][j]);
                }
//              printf("\n");
        }
}

int main(int argc, char **argv){
        int val_returned = 0;
        if(argc == 2 && ( (strcmp(argv[1], "--help")==0) || (strcmp(argv[1], "-h")==0) )){
                val_returned = help_func();
        }
        else{
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
                                        if(argc == 5){
                                                num_thread = atoi(argv[4]);
                                        }
                                }
                        }
                }
                printf("\nMatrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d)\n", row_a, col_a, row_b, col_b, row_a, col_b);

                struct timeval tstart, tstop;
                float elapsed = 0.f;

                float start_time, run_time;

                omp_set_num_threads(num_thread);
                start_time = omp_get_wtime();

                gettimeofday(&tstart,NULL);

                float matrix_a[row_a][col_a];
                float matrix_b[row_b][col_b];
                float matrix_c[row_a][col_b];

                //fill matrix A and B with random float, range 0-1
//                srand(time(NULL));

//              printf("\n## Matrix A:\n");
                fill_matrix(row_a, col_a, matrix_a);
//              printf("\n## Matrix B:\n");
                fill_matrix(row_b, col_b, matrix_b);
                //matrix multiplication
//              printf("\n## Matrix C:\n");
                int i, j, k;
                #pragma omp parallel for  shared(matrix_a,matrix_b,matrix_c) private (i,j,k)
                for(i=0; i<row_a; i++){
                        for(j=0; j<col_b; j++){
                                //float sum_el = 0.f;
                                for(k=0; k<col_a; k++){
                                       // sum_el += matrix_a[i][k]*matrix_b[k][j];
                                        matrix_c[i][j] += matrix_a[i][k]*matrix_b[k][j];
                                }
                                //matrix_c[i][j] = sum_el;
//                              printf("%f ", matrix_c[i][j]);
                        }
//                      printf("\n");
                }
                gettimeofday(&tstop,NULL);
                run_time = omp_get_wtime() - start_time;

                printf("\nTerminated.\n");
                elapsed = (tstop.tv_sec - tstart.tv_sec) + ((tstop.tv_usec - tstart.tv_usec)/1000000.0);
                printf("Data processing with %d threads in %.6f s using \"gettimeofday\" and  %.6f s using \"OMP timer\".\n\n", num_thread,elapsed,run_time);
        }
        return val_returned;
}
