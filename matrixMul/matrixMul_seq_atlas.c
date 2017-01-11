// System includes
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <assert.h>

#include "helper_string.h"

//for blas:
#include <stdlib.h>

extern "C"
{
   int sgemm_(char *, char *, int *, int *, int *, float *, float *, int *, 
              float *, int *, float *, float *, int *);
}


//le blas lavorano sulle COLONNE, quindi riempio le matrici PER COLONNE:
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


int matrixMultiply(int argc, char **argv, int &sizeMatrix, int &perf, int &verb)
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
    char transA = 'N', transB = 'N';
    float one = 1.f, zero = 0.f;
    if(perf){
    	printf("\n----- Performing warmup operation... (perf = 1)\n");
	sgemm_(&transA, &transB, &sizeMatrix, &sizeMatrix, &sizeMatrix, &one, h_A, &sizeMatrix, h_B, &sizeMatrix, &zero, h_C, &sizeMatrix);
    	printf("done\n");
    }
    
    // Execute the kernel
    int nIter = 1;

    struct timeval start, stop;

   // for (int j = 0; j < nIter; j++)
   // {
//	int ii, jj;

    gettimeofday(&start,NULL);

    //sgemm_(&transA, &transB, &rowsA, &colsB, &common, &one, A, &rowsA, B, &common, &zero, C, &rowsA);

    sgemm_(&transA, &transB, &sizeMatrix, &sizeMatrix, &sizeMatrix, &one, h_A, &sizeMatrix, h_B, &sizeMatrix, &zero, h_C, &sizeMatrix);

/*	for(ii=0; ii<sizeMatrix; ii+=block_size){
		for(jj=0; jj<sizeMatrix; jj+=block_size){
			matrixMul_seq(h_C, h_A, h_B, sizeMatrix, ii, jj,block_size);
		}
	}*/
     gettimeofday(&stop,NULL);
   // }

    double secTotal = 0.0;
    secTotal = (stop.tv_sec - start.tv_sec) + ((stop.tv_usec - start.tv_usec)/1000000.0);
    
    // Compute and print the performance
    double secPerMatrixMul = secTotal / nIter;
    double flopsPerMatrixMul = 2.0 * (double)sizeMatrix * (double)sizeMatrix * (double)sizeMatrix;
   // double flopsPerMatrixMul1 = (2.0*(double)block_size + 1.0)*(((double)sizeMatrix * (double)sizeMatrix * (double)sizeMatrix)/(double)block_size);
    printf("flopsPerMatrixMul = %f\n", flopsPerMatrixMul);
    double gigaFlops = (flopsPerMatrixMul * 1.0e-9f) / (secPerMatrixMul );
   // double gigaFlops1 = (flopsPerMatrixMul1 * 1.0e-9f) / secPerMatrixMul;
    printf("A. Performance= %.2f GFlop/s, Time= %.8f sec, Size= %.0f Ops\n",gigaFlops,secPerMatrixMul,flopsPerMatrixMul);
   // printf("B. Performance= %.2f GFlop/s, Time= %.8f sec, Size= %.0f Ops\n",gigaFlops1,secPerMatrixMul,flopsPerMatrixMul1);

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

//    int block_size = 2;
    int sizeMatrix = 512;

    int perf = 0;
    int verb = 0;


    if (checkCmdLineFlag(argc, (const char **)argv, "help") ||
        checkCmdLineFlag(argc, (const char **)argv, "?"))
    {
        printf("\nUsage: -size=SizeMatrix(d:%d)\n",sizeMatrix);
        printf("       -perf=Performance(d:%d)   -v=Verbose(d:%d)\n\n",perf,verb);

        exit(EXIT_SUCCESS);
    }

    // width of Matrix A
 /*   if (checkCmdLineFlag(argc, (const char **)argv, "subSize"))
    {
        block_size = getCmdLineArgumentInt(argc, (const char **)argv, "subSize");
    }
*/
    //printf("\n++++ block_size: %d ++++\n", block_size);

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

/*    if (sizeMatrix % block_size != 0)
    {
        printf("Error: The size matrix must be a multiple of the size sub-matrix!\n");
        exit(EXIT_FAILURE);
    }
*/
    printf("#### MatrixA(%d,%d), MatrixB(%d,%d) ####\n", sizeMatrix, sizeMatrix, sizeMatrix, sizeMatrix);

    int matrix_result = matrixMultiply(argc, argv, sizeMatrix, perf, verb);

    exit(matrix_result);
}
