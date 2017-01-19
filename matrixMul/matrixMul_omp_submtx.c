// System includes
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <assert.h>

// Helper functions and utilities to work with CUDA
//#include <helper_functions.h>
#include "helper_string.h"

// omp
#include <omp.h>

void matrixMul_seq(float *C, float *A, float *B, int sizeMatrix, int BLOCK_SIZE)
{
 int cur_row, cur_col;

 #pragma omp parallel shared(A,B,C) private(cur_row,cur_col)
 {
 // #pragma omp master
 //  printf("\n\n xxxxx get_omp_threads: %d xxxxx\n\n", omp_get_num_threads());
 
  #pragma omp for collapse(2) schedule(static)
  for(cur_row=0; cur_row<sizeMatrix; cur_row+=BLOCK_SIZE){
      for(cur_col=0; cur_col<sizeMatrix; cur_col+=BLOCK_SIZE){

    int bx = cur_col/BLOCK_SIZE;
    int by = cur_row/BLOCK_SIZE;

    // Index of the first sub-matrix of A processed by the block
    int aBegin = sizeMatrix * BLOCK_SIZE * by;

    // Index of the last sub-matrix of A processed by the block
    int aEnd   = aBegin + sizeMatrix - 1;

    // Step size used to iterate through the sub-matrices of A
    int aStep  = BLOCK_SIZE;

    // Index of the first sub-matrix of B processed by the block
    int bBegin = BLOCK_SIZE * bx;

    // Step size used to iterate through the sub-matrices of B
    int bStep  = BLOCK_SIZE * sizeMatrix;

    float Cs[BLOCK_SIZE][BLOCK_SIZE];

     int sr1,sc1;
     #pragma unroll
     for(sr1=0; sr1<BLOCK_SIZE; sr1++){
            for(sc1=0; sc1<BLOCK_SIZE; sc1++){
                    Cs[sr1][sc1] = 0.0;
            }
     }

    // Loop over all the sub-matrices of A and B
    // required to compute the block sub-matrix
    for (int a = aBegin, b = bBegin;
         a <= aEnd;
         a += aStep, b += bStep)
    {
        // Declaration of the shared memory array As used to
        // store the sub-matrix of A
        /*__shared__ */float As[BLOCK_SIZE][BLOCK_SIZE];

        // Declaration of the shared memory array Bs used to
        // store the sub-matrix of B
        /*__shared__ */float Bs[BLOCK_SIZE][BLOCK_SIZE];

        // Load the matrices from device memory
        // to shared memory; each thread loads
        // one element of each matrix
	int sr,sc;
	for(sr=0; sr<BLOCK_SIZE; sr++){
		for(sc=0; sc<BLOCK_SIZE; sc++){
        		As[sr][sc] = A[a + sizeMatrix * sr + sc];
        		Bs[sr][sc] = B[b + sizeMatrix * sr + sc];
		}
	}

        // Multiply the two matrices together;
        // each thread computes one element
        // of the block sub-matrix
	
	#pragma unroll
	for(sr=0; sr<BLOCK_SIZE; sr++){
		for(sc=0; sc<BLOCK_SIZE; sc++){
			float Csub = 0;
        		for (int k = 0; k < BLOCK_SIZE; ++k)
        		{
            			Csub += As[sr][k] * Bs[k][sc];
        		}
			Cs[sr][sc] += Csub;
		}
	}
    }

    // Write the block sub-matrix to device memory;
    // each thread writes one element
    //int sr,sc;
    int c = sizeMatrix * BLOCK_SIZE * by + BLOCK_SIZE * bx;
    for(sr1=0; sr1<BLOCK_SIZE; sr1++){
	for(sc1=0; sc1<BLOCK_SIZE; sc1++){
    		C[c + sizeMatrix * sr1 + sc1] = Cs[sr1][sc1];
	}
    }

  }// fine cilco su cur_col
 }// fine ciclo su cur_row
}// fine regione parallela

}// fine funzione

void constantInit(float *data, int size, float val)
{
    for (int i = 0; i < size; ++i)
    {
        data[i] = val*i;
    }
}

void printMatrix(float *data, int size){
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			printf("%.0f ", data[i*size+j]);
		}
		printf("\n");
	}
}


/**
 * Run a simple test of matrix multiplication using CUDA
 */
int matrixMultiply(int argc, char **argv, int block_size, int &sizeMatrix, int &perf, int &verb, int &max_th)
{
  //threads:
  int nth;
  for(nth=2; nth<=max_th; nth++){

    omp_set_num_threads(nth);

    // Allocate host memory for matrices A and B
    unsigned int size_A = sizeMatrix * sizeMatrix;
    unsigned int mem_size_A = sizeof(float) * size_A;
    float *h_A = (float *)malloc(mem_size_A);
    unsigned int size_B = sizeMatrix * sizeMatrix;
    unsigned int mem_size_B = sizeof(float) * size_B;
    float *h_B = (float *)malloc(mem_size_B);

    // Initialize host memory
    const float valB = 10.0f;
    constantInit(h_A, size_A, 1.0f);
    constantInit(h_B, size_B, valB);

    // Allocate host matrix C
    unsigned int mem_size_C = sizeMatrix * sizeMatrix * sizeof(float);
    float *h_C = (float *) malloc(mem_size_C);

    int ii, jj;

    if(perf){
//    	printf("\n----- Performing warmup operation... (perf = 1)\n");
        matrixMul_seq(h_C, h_A, h_B, sizeMatrix,block_size);
  //  	printf("done\n");
    }
    
    // Execute the kernel
    int nIter = 1;

    struct timeval start, stop;

    gettimeofday(&start,NULL);
    matrixMul_seq(h_C, h_A, h_B, sizeMatrix,block_size);
    gettimeofday(&stop,NULL);

    double secTotal = 0.0;
    secTotal = (stop.tv_sec - start.tv_sec) + ((stop.tv_usec - start.tv_usec)/1000000.0);
    
    // Compute and print the performance
    double secPerMatrixMul = secTotal / nIter;
    double flopsPerMatrixMul = 2.0 * (double)sizeMatrix * (double)sizeMatrix * (double)sizeMatrix;
    double flopsPerMatrixMul1 = (2.0*(double)block_size + 1.0)*(((double)sizeMatrix * (double)sizeMatrix * (double)sizeMatrix)/(double)block_size);
   // printf("flopsPerMatrixMul = %f\n", flopsPerMatrixMul);
    double gigaFlops = (flopsPerMatrixMul * 1.0e-9f) / (secPerMatrixMul );
    double gigaFlops1 = (flopsPerMatrixMul1 * 1.0e-9f) / secPerMatrixMul;
   // printf("A. Performance= %.2f GFlop/s, Time= %.8f sec, Size= %.0f Ops\n",gigaFlops,secPerMatrixMul,flopsPerMatrixMul);
    printf("%d %.2f %.8f\n",nth, gigaFlops1,secPerMatrixMul);

    if(verb){
        //print matrix
        printf("\n##### A #####\n");
         printMatrix(h_A, sizeMatrix);
        printf("\n##### B #####\n");
        printMatrix(h_B, sizeMatrix);
    	printf("\n##### C #####\n");
    	printMatrix(h_C, sizeMatrix);
    }

    // Clean up memory
    free(h_A);
    free(h_B);
    free(h_C);

  }// fine cilco sui threads omp
  return EXIT_SUCCESS;
}


/**
 * Program main
 */
int main(int argc, char **argv)
{
//    printf("[Matrix Multiply Using CUDA] - Starting...\n");

    if (checkCmdLineFlag(argc, (const char **)argv, "help") ||
        checkCmdLineFlag(argc, (const char **)argv, "?"))
    {
        printf("Usage: -size=SizeMatrix(d:512)  -subSize=SizeSubMatrix(d:64)\n");
	printf("       -maxTh=MaxNumberOfThreads(d:4)\n");
        printf("       -perf=Performance(d:0)   -v=Verbose(d:0)\n\n");
        printf("Note:  Size matrix must be a multiple of the size sub-matrix.\n\n");

        exit(EXIT_SUCCESS);
    }

    // Use a larger block size for Fermi and above
   // int block_size = (deviceProp.major < 2) ? 16 : 32;
    int block_size = 64;
    int max_th = 4;

    // width of Matrix A
    if (checkCmdLineFlag(argc, (const char **)argv, "subSize"))
    {
        block_size = getCmdLineArgumentInt(argc, (const char **)argv, "subSize");
    }

    if (checkCmdLineFlag(argc, (const char **)argv, "maxTh"))
    {
        max_th = getCmdLineArgumentInt(argc, (const char **)argv, "maxTh");
    }

    int sizeMatrix = 512;

    int perf = 0;
    int verb = 0;

     if(checkCmdLineFlag(argc, (const char **)argv, "v")){
        verb = getCmdLineArgumentInt(argc, (const char **)argv, "v");
    }

     if(checkCmdLineFlag(argc, (const char **)argv, "perf")){
	perf = getCmdLineArgumentInt(argc, (const char **)argv, "perf");
    }

    if (checkCmdLineFlag(argc, (const char **)argv, "size"))
    {
        sizeMatrix = getCmdLineArgumentInt(argc, (const char **)argv, "size");
    }

    if (sizeMatrix % block_size != 0)
    {
        printf("Error: The size matrix must be a multiple of the size sub-matrix!\n");
        exit(EXIT_FAILURE);
    }

    printf("# ###################################### #\n# MatrixA(%d,%d), MatrixB(%d,%d)\n", sizeMatrix, sizeMatrix, sizeMatrix, sizeMatrix);
    printf("# From 2 To %d threads\n", max_th);
    printf("# Block_size: %d\n",block_size);
    if(perf)
        printf("# PERF ENABLED\n");
    else
	printf("# PERF DISABLED\n");
    printf("# ###################################### #\n\n");


    int matrix_result = matrixMultiply(argc, argv, block_size, sizeMatrix, perf, verb, max_th);

    exit(matrix_result);
}
