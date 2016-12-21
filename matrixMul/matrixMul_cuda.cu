#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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
        printf("\nUsage:   -rA=RowsA(d:512)     -cA=ColumnsA(d:512)  -cB=ColumnsB(d:512) | matrix(row,col), ColumnsA = RowsB\n");
	printf("         -db=DimBlock(d:16)                             | DimBlock(in threads): [1-32], block(DimBlock, DimBlock)\n");
        printf("         -w=WarmUpData(d:0)\n");
        printf("         -v=Verbose(d:0)\n\n");
       // printf("Default: A = (512,512) B = (512,512); DIM_BLOCK = 16; WARMUP = 0; VERBOSE = 0\n\n");

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
		int debug = 0, perf = 0;
		int dimBlock = 16;

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
		if (checkCmdLineFlag(argc, (const char **)argv, "db")){
                        dimBlock = getCmdLineArgumentInt(argc, (const char **)argv, "db");
                }


                if(dimBlock>32){
                        val_returned =  help_func();
                }
                else{
                        cudaEvent_t start_all, end_all,
					start_fill, end_fill,
					start_mm, end_mm,
					start_free, end_free;

			cudaEventCreate(&start_all);
                        cudaEventCreate(&end_all);
			cudaEventCreate(&start_fill);
                        cudaEventCreate(&end_fill);
			cudaEventCreate(&start_mm);
                        cudaEventCreate(&end_mm);
			cudaEventCreate(&start_free);
                        cudaEventCreate(&end_free);

                        cudaEventRecord(start_all, 0);

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

			cudaEventRecord(end_all, 0);
                        cudaEventSynchronize(end_all);
			float timeAlloc;
                        cudaEventElapsedTime(&timeAlloc, start_all, end_all);

                        printf("\n### Matrix A = (%d,%d); Matrix B = (%d,%d); AxB = (%d,%d);\n",
					row_a, col_a, col_a, col_b, row_a, col_b);
                        printf("### dimBlock = %d; gridA(%d,%d); gridB(%d,%d); gridC(%d,%d);\n",
					dimBlock,
					(int)ceil(col_a/(float)dimBlock),(int)ceil(row_a/(float)dimBlock),
					(int)ceil(col_b/(float)dimBlock),(int)ceil(col_a/(float)dimBlock),
					(int)ceil(col_b/(float)dimBlock),(int)ceil(row_a/(float)dimBlock)
				);

			if(debug){                        
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


			cudaEventRecord(start_fill, 0);

                        matrixFillKernel<<<gridA,block>>>(row_a,col_a,d_matrix_a);
                        matrixFillKernel<<<gridB,block>>>(col_a,col_b,d_matrix_b);

			cudaEventRecord(end_fill, 0);
                        cudaEventSynchronize(end_fill);
			float timeFill;
                        cudaEventElapsedTime(&timeFill, start_fill, end_fill);

			if(perf){
				//Performs warmup operation
				printf("\nPreforming warmup...\n");
				matrixMulKernel<<<gridC,block>>>(row_a,col_a,col_b,d_matrix_a,d_matrix_b,d_matrix_c);
			}

			printf("\nComputing matrix multimplication...\n");
			cudaEventRecord(start_mm, 0);
                        matrixMulKernel<<<gridC,block>>>(row_a,col_a,col_b,d_matrix_a,d_matrix_b,d_matrix_c);

			cudaEventRecord(end_mm, 0);
                        cudaEventSynchronize(end_mm);
			float timeMtxMul;
                        cudaEventElapsedTime(&timeMtxMul, start_mm, end_mm);

			if(debug == 2){
                        	cudaMemcpy(matrix_a, d_matrix_a, (row_a*col_a)*sizeof(float), cudaMemcpyDeviceToHost);
                        	cudaMemcpy(matrix_b, d_matrix_b, (col_a*col_b)*sizeof(float), cudaMemcpyDeviceToHost);
                        	cudaMemcpy(matrix_c, d_matrix_c, (row_a*col_b)*sizeof(float), cudaMemcpyDeviceToHost);
			}


			cudaEventRecord(start_free, 0);

                        cudaFree(d_matrix_c); cudaFree(d_matrix_a); cudaFree(d_matrix_b);

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

			cudaEventRecord(end_free, 0);
                        cudaEventSynchronize(end_free);
			float timeFree;
                        cudaEventElapsedTime(&timeFree, start_free, end_free);

                        printf("\nTerminated.\n");

			double flops4mtxmul = 2.0*(double)row_a*(double)col_a*(double)col_b;
			double gigaFlops = (flops4mtxmul * 1.0e-9f) / (timeMtxMul / 1000.0f);

		//	printf("\nPerformance: %f GFlop/s; Time: %f ms; Flop: %f\n\n", gigaFlops, timeMtxMul, flops4mtxmul);
			if(debug){
                        	printf("\nDimBlock: %d,  Flop: %.0f,  GFlop: %f GFlop/s,  Time_mtxMul: %f s\n",
                                	dimBlock, flops4mtxmul, gigaFlops, timeMtxMul/1000.0f);
				printf("Time_tot:");
				printf("\nA. timeAllocation: %f s;\nB. timeComputation: %f s (fill: %f s, matrixMul: %f s);\nC. timeFree: %f s;\nD. TOTAL: %f s\n\n",
					timeAlloc/1000.0f,
					(timeFill+timeMtxMul)/1000.0f,
					timeFill/1000.0f,
					timeMtxMul/1000.0f,
					timeFree/1000.0f,
					(timeAlloc+timeFill+timeMtxMul+timeFree)/1000.0f
				);
			}
                	else
                        	 printf("\n%d %.0f %f %f %f\n\n",
                                        dimBlock, flops4mtxmul, gigaFlops, timeMtxMul/1000.0f, (timeAlloc+timeFill+timeMtxMul+timeFree)/1000.0f);
                }
        }
        return val_returned;
}
