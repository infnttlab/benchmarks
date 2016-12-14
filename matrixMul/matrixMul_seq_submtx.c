// System includes
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <assert.h>

// Helper functions and utilities to work with CUDA
//#include <helper_functions.h>
#include "helper_string.h"

/**
 * Matrix multiplication (CUDA Kernel) on the device: C = A * B
 * wA is A's width and wB is B's width
 */
//template<int BLOCK_SIZE> void
void matrixMul_seq(float *C, float *A, float *B, int sizeMatrix, int cur_row, int cur_col, int BLOCK_SIZE)
//matrixMulCUDA(float *C, float *A, float *B, int wA, int wB)
{
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
}

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
int matrixMultiply(int argc, char **argv, int block_size, int &sizeMatrix, int &perf, int &verb)
{
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
    	printf("\n----- Performing warmup operation... (perf = 1)\n");
	for(ii=0; ii<sizeMatrix; ii+=block_size){
                for(jj=0; jj<sizeMatrix; jj+=block_size){
                        matrixMul_seq(h_C, h_A, h_B, sizeMatrix, ii, jj,block_size);
                }
        }
    	printf("done\n");
    }
    
    // Execute the kernel
    int nIter = 1;

    struct timeval start, stop;

   // for (int j = 0; j < nIter; j++)
   // {
//	int ii, jj;
	gettimeofday(&start,NULL);
	for(ii=0; ii<sizeMatrix; ii+=block_size){
		for(jj=0; jj<sizeMatrix; jj+=block_size){
			matrixMul_seq(h_C, h_A, h_B, sizeMatrix, ii, jj,block_size);
		}
	}
	gettimeofday(&stop,NULL);
   // }

    double secTotal = 0.0;
    secTotal = (stop.tv_sec - start.tv_sec) + ((stop.tv_usec - start.tv_usec)/1000000.0);
    
    // Compute and print the performance
    double secPerMatrixMul = secTotal / nIter;
    double flopsPerMatrixMul = 2.0 * (double)sizeMatrix * (double)sizeMatrix * (double)sizeMatrix;
    printf("flopsPerMatrixMul = %f\n", flopsPerMatrixMul);
    double gigaFlops = (flopsPerMatrixMul * 1.0e-9f) / (secPerMatrixMul );
    printf(
        "Performance= %.2f GFlop/s, Time= %.8f sec, Size= %.0f Ops\n",
        gigaFlops,
        secPerMatrixMul,
        flopsPerMatrixMul);

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
        printf("Usage: -device=n (n >= 0 for deviceID)\n");
        printf("       -size=SizeMatrix  -subSize=SizeSubMatrix\n");
        printf("       -perf=Performance -v=Verbose\n\n");
        printf("Note:  Size matrix must be a multiple of the size sub-matrix.\n\n");

        exit(EXIT_SUCCESS);
    }

    // Use a larger block size for Fermi and above
   // int block_size = (deviceProp.major < 2) ? 16 : 32;
    int block_size = 2;

    // width of Matrix A
    if (checkCmdLineFlag(argc, (const char **)argv, "subSize"))
    {
        block_size = getCmdLineArgumentInt(argc, (const char **)argv, "subSize");
    }

    printf("\n++++ block_size: %d ++++\n", block_size);

    int sizeMatrix = 5*2*block_size;

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

    printf("#### MatrixA(%d,%d), MatrixB(%d,%d) ####\n", sizeMatrix, sizeMatrix, sizeMatrix, sizeMatrix);

    int matrix_result = matrixMultiply(argc, argv, block_size, sizeMatrix, perf, verb);

    exit(matrix_result);
}
