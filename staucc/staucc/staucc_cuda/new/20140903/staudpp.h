#ifndef STAUDPP_H
#define STAUDPP_H

#define numSimsMin              1
#define numSimsMax              1024
#define numMemsMin    	        1

#define	numTrials		1
/////////////////////////
#define numSims		        32  // numero di blocchi lanciati = numero di simulazioni
#define threadsPerBlock         256 // 27 numero di thread per blocco (che si sparticono le membrane di una simulazione)
#define numMems			1000 // 27 ogni thread computa (numMems/threadsPerBlock) membrane, che quindi deve essere intero
#define numSpecs 		14 // il numero di colonne di left_side
#define numReacts 		40 // il numero di righe di left_side
/////////////////////////
/*
#define numSims		        1  // numero di blocchi lanciati = numero di simulazioni
#define threadsPerBlock         27 // numero di thread per blocco (che si sparticono le membrane di una simulazione)
#define numMems			27 // ogni thread computa (numMems/threadsPerBlock) membrane, che quindi deve essere intero
#define numSpecs 		10 // il numero di colonne di left_side
#define numReacts 		69 // il numero di righe di left_side
*/
/*
#define numSims		            32 // con time_max = 10
#define	numTrials				1
#define numMems			        256
#define threadsPerBlock         64
#define numSpecs 		        16
#define numReacts 		        32
*/

/*
#define numSims		            32 // con time_max = 1
#define	numTrials				1
#define numMems			        4
#define threadsPerBlock         4
#define numSpecs 		        3
#define numReacts 		        3
*/
/*
#define numSims		            1 // con time_max = 10
#define	numTrials				1
#define numMems			        36
#define threadsPerBlock         32
#define numSpecs 		        14
*/
/*
#define numSims		            1 // con time_max = 10
#define	numTrials				1
#define numMems			        144
#define threadsPerBlock         16
#define numSpecs 		        14
#define numReacts				32
*/
/*
#define numSims		            1 // con time_max = 10
#define numTrials 				1
#define numMems			        900
#define threadsPerBlock         192
#define numSpecs 		        14
*/
/*
#define numSims		        	1 // time_max = 10
#define numTrials 				1
#define threadsPerBlock			640
#define numMems			   	8100
#define numSpecs 		    14
*/

#include <stdio.h>
#include <float.h>
#include <string>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>
#include <curand_kernel.h>

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

#endif
