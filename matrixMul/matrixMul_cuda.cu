#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

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
        printf("\nUsage: ./a.out <ROW_A> <COL_A> <COL_B> <DIM_BLOCK>\n");
        printf("Where:   MATRIX(ROW, COL)  and  <COL_A> == <ROW_B>\n");
        printf("         DIM_BLOCK: [1-32]; BLOCK(dimBlock, dimBlock)\n\n");
        printf("Default: A = (512,512); B = (512,512); DIM_BLOCK = 16;\n\n");

        return 0;
}

__global__
void matrixFillKernel(int rowMax, int colMax,  float *d_matrix){
        int col = blockIdx.x * blockDim.x + threadIdx.x;
        int row = blockIdx.y * blockDim.y + threadIdx.y;
        int i = row*colMax+col;
        if(row<rowMax && col<colMax){
                d_matrix[i] = 0.1f;
        }
}

__global__
void matrixMulKernel(int row_a, int col_a, int col_b, float* d_matrix_a, float* d_matrix_b, float* d_matrix_c){
        int col = blockIdx.x * blockDim.x + threadIdx.x;
        int row = blockIdx.y * blockDim.y + threadIdx.y;

        int k;

	d_matrix_c[row*col_b + col] = 0.f;
        if(col<col_b && row<row_a){
                for(k=0; k<col_a; k++){
                        d_matrix_c[row*col_b + col] += d_matrix_a[row*col_a + k] * d_matrix_b[k*col_b + col];
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
		int debug = 0;
                if(argc >= 2){
                        // change ROW_A
                        row_a = atoi(argv[1]);

                        if(argc >= 3){
                                // change ROW_A COL_A and ROW_B, where COL_A = ROW_B
                                col_a = atoi(argv[2]);
                                if(argc >= 4){
                                        col_b = atoi(argv[3]);
					if(argc == 5)
						debug = atoi(argv[4]);
                                }
                        }
                }

                int dimBlock = 16;
                if(argc >= 5){
                        dimBlock = atoi(argv[4]);
                }

                if(dimBlock>32){
                        val_returned =  help_func();
                }
                else{
                        struct timeval tstart, tstop;
                        double elapsed = 0.f;

                        cudaEvent_t startCUDA, stopCUDA;
                        double timeCUDA;

                        cudaEventCreate(&startCUDA);
                        cudaEventCreate(&stopCUDA);

                        gettimeofday(&tstart,NULL);
                        cudaEventRecord(startCUDA, 0);

                        float *matrix_a = (float*)malloc(row_a*col_a * sizeof(float));
                        float *matrix_b = (float*)malloc(col_a*col_b * sizeof(float));
                        float *matrix_c = (float*)malloc(row_a*col_b * sizeof(float));

                        float *d_matrix_a;
                        float *d_matrix_b;
                        float *d_matrix_c;

                        cudaMalloc((void**)&d_matrix_a, (row_a*col_a * sizeof(float)));
                        cudaMalloc((void**)&d_matrix_b, (col_a*col_b * sizeof(float)));
                        cudaMalloc((void**)&d_matrix_c, (row_a*col_b * sizeof(float)));

                        dim3 block(dimBlock,dimBlock);

                        //matrix with theadsPerRow = row_a threadsPerCol = col_a
                        //              dimBlock.x = col                        dimBlock.y = row
                        dim3 gridA( (int)ceil(col_a/(float)dimBlock) , (int)ceil(row_a/(float)dimBlock) );
                        dim3 gridB( (int)ceil(col_b/(float)dimBlock) , (int)ceil(col_a/(float)dimBlock)  );
                        dim3 gridC( (int)ceil(col_b/(float)dimBlock) , (int)ceil(row_a/(float)dimBlock)  );

			if(debug){
                        	printf("\n### Matrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d);\n",
					row_a, col_a, col_a, col_b, row_a, col_b);
                        	printf("### dimBlock = %d; gridA(%d,%d); gridB(%d,%d); gridC(%d,%d);\n",
					dimBlock, (int)ceil(col_a/(float)dimBlock),(int)ceil(row_a/(float)dimBlock),
					(int)ceil(col_b/(float)dimBlock),(int)ceil(col_a/(float)dimBlock),
					(int)ceil(col_b/(float)dimBlock) , (int)ceil(row_a/(float)dimBlock));
                        
				int col_gA = (int)ceil(col_a/(float)dimBlock);
        	                int row_gA = (int)ceil(row_a/(float)dimBlock);
                	        int col_gB = (int)ceil(col_b/(float)dimBlock);
                        	int row_gB = (int)ceil(col_a/(float)dimBlock);
	                        int col_gC = (int)ceil(col_b/(float)dimBlock);
        	                int row_gC = (int)ceil(row_a/(float)dimBlock);

	                        int totThA = col_gA*row_gA*dimBlock*dimBlock;
        	                int totThB = col_gB*row_gB*dimBlock*dimBlock;
                	        int totThC = col_gC*row_gC*dimBlock*dimBlock;

	                        printf("\n******************** TEST ***********************\n");
        	                printf("- totThA = %d VS totElA = %d\n", totThA,col_a*row_a);
                	        printf("- totThB = %d VS totElB = %d\n", totThB,col_a*col_b);
	                        printf("- totThC = %d VS totElC = %d\n", totThC,col_b*row_a);
        	                printf("*************************************************\n");
			}

                        matrixFillKernel<<<gridA,block>>>(row_a,col_a,d_matrix_a);
                        matrixFillKernel<<<gridB,block>>>(col_a,col_b,d_matrix_b);

                        matrixMulKernel<<<gridC,block>>>(row_a,col_a,col_b,d_matrix_a,d_matrix_b,d_matrix_c);

			if(debug){
                        	cudaMemcpy(matrix_a, d_matrix_a, (row_a*col_a)*sizeof(float), cudaMemcpyDeviceToHost);
                        	cudaMemcpy(matrix_b, d_matrix_b, (col_a*col_b)*sizeof(float), cudaMemcpyDeviceToHost);
                        	cudaMemcpy(matrix_c, d_matrix_c, (row_a*col_b)*sizeof(float), cudaMemcpyDeviceToHost);
			}

                        cudaFree(d_matrix_c); cudaFree(d_matrix_a); cudaFree(d_matrix_b);

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

                        cudaEventRecord(stopCUDA, 0);
                        cudaEventSynchronize(stopCUDA);
                        cudaEventElapsedTime(&timeCUDA, startCUDA, stopCUDA);

                        gettimeofday(&tstop,NULL);
                        printf("\nTerminated.\n");
                        elapsed = (tstop.tv_sec - tstart.tv_sec) + ((tstop.tv_usec - tstart.tv_usec)/1000000.0);

                        printf("Data processing in %f s using \"gettimeofday\" and %f ms using \"CUDA Events\".\n\n",
				elapsed, timeCUDA);
                }
        }
        return val_returned;
}
