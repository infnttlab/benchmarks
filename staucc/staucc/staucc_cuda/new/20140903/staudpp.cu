#include "staudpp.h"
//nvcc -Wall -I/home/cuda-5.0/samples/common/inc -gencode arch=compute_20,code=sm_20 staudpp.cu -o staudpp
//nvcc  -I/home/cuda-5.0/samples/common/inc -gencode arch=compute_20,code=sm_20 staudpp.cu -o staudpp

// Scalar parameters

__constant__ float stepMax;
float stepMax_h;
__constant__ unsigned tgtKind;
unsigned tgtKind_h;

// Vector parameters

float *Mfeed_h, *Mfeed;
float *molSize_h, *molSize;
float *membSize_h, *membSize;

unsigned *leftSide_h, *leftSide;
unsigned *rightSide_h, *rightSide;
	
int *var_h, *var;
int *varSend_h, *varSend;	// varSend used only if tgtMatrix used (tgtKind=0)

int *tgtVector_h, *tgtVector;
int *tgtMatrix_h, *tgtMatrix;

float *Cvector_h, *Cvector;

// Output variables

FILE *log_fPtr;
FILE **buffer_fPtr;

__constant__ bool buffBool;
bool buffBool_h;
__constant__ unsigned buffEvery;
unsigned buffEvery_h;
__constant__ unsigned buffRows;
unsigned buffRows_h;
__constant__ unsigned buffCols;
unsigned buffCols_h;
__constant__ float eps;
float eps_h;

unsigned *buffRowCounterAllSims_h, *buffRowCounterAllSims;
float *buffMAllSims_h, *buffMAllSims;

// Simulation variables

unsigned long *seed_h, *seed;
curandState *rngStatesAllSims;

unsigned *step_h, *step;
float *t_h, *t;

float *MAllSims_h, *MAllSims;
float *MbkpAllSims;
unsigned *k_ruleAllSims;
unsigned *criticalAllSims;
unsigned *orderAllSims;
float *R_probAllSims;
float *R_probCAllSims;
float *M_sendAllSims; // ex SHARED memory
float *t1AllSims; // ex SHARED memory
float *tauAllSims; // ex SHARED memory
float *FSAllSims; // ex SHARED memory
float *FSbkpAllSims;
float *tauSSAAllSims;
unsigned *flagSSAAllSims;
unsigned *flagStepAllSims;

float *tauPreMinAllSims;
float *tau1AllSims;
float *tau2AllSims;
float *a0cAllSims;
float *a0AllSims;

unsigned *HORAllSims; // ex float HOR
float *auxVec1AllSims; // ex float mu, unsigned (anche float va bene) current
float *auxVec2AllSims; // ex float sigma





///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// DEVICE code /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////          

__global__ void initRNG(curandState *rngStatesAllSims, unsigned long *seed) {
    
	// membrane partition among threads
	int numMemsPerThread = numMems / blockDim.x;

	curandState *rngStates = rngStatesAllSims + blockIdx.x*numMems;

	for (int m=0; m<numMemsPerThread; m++)
		curand_init(seed[blockIdx.x], m*blockDim.x + threadIdx.x, 0, &rngStates[m*blockDim.x + threadIdx.x]);
}

__global__ void dynamicStep (curandState *rngStatesAllSims, 
  float *buffMAllSims, unsigned *buffRowCounterAllSims,
  unsigned *step, float *MAllSims, float *t, 
  unsigned *leftSide, unsigned *rightSide,
  int *var, int *varSend,
  int *tgtVector, int *tgtMatrix,
  float *Cvector,
  float *Mfeed, float *molSize, float *membSize,
  unsigned *k_ruleAllSims, unsigned *criticalAllSims, unsigned *orderAllSims, float *R_probAllSims, float *R_probCAllSims,
  float *M_sendAllSims, float *t1AllSims, float *tauAllSims, float *FSAllSims,
  float *FSbkpAllSims, 
  float *tauSSAAllSims, 
  float *MbkpAllSims,
  unsigned *flagSSAAllSims, unsigned *flagStepAllSims,
  float *tauPreMinAllSims, float *tau1AllSims, float *tau2AllSims, float *a0cAllSims, float *a0AllSims,
  unsigned *HORAllSims, float *auxVec1AllSims, float *auxVec2AllSims) {
  	  
  unsigned s = blockIdx.x;
  unsigned bDim = blockDim.x;
  unsigned thIdx = threadIdx.x;

  // membrane partition among threads
  unsigned numMemsPerThread = numMems / bDim;

  // per-thread scalar variables

  unsigned localStep;
  float localT;

  int m;
  int r, c, e;

  float aux1; // ex float a0cTmp, float numR, float maxTmp, unsigned (anche float va bene) currentTmp	       
  float aux2; // ex float minNumR, float tau1Tmp, unsigned HORTmp
  // anziche' avere unsigned HORTmp e float HOR ho fatto il contrario per avere meno variabili
  // tanto la conversione float -> unsigned dovrebbe farla lo stesso...

  //float rnd;
  float rnd_t, alpha;

  float splitShared; // ex t3, cn, globmintau, minFS

  curandState *rngStates = rngStatesAllSims + s*numMems;

  float *buffM = buffMAllSims + s*numMems*buffRows*buffCols;
  unsigned *buffRowCounter = buffRowCounterAllSims + s*numMems;

  unsigned *k_rule = k_ruleAllSims + s*numMems*numReacts;
  unsigned *critical = criticalAllSims + s*numMems*numReacts;
  unsigned *order = orderAllSims + s*numMems*numReacts;
  float *R_prob = R_probAllSims + s*numMems*numReacts;
  float *R_probC = R_probCAllSims + s*numMems*numReacts;

  float *M_send = M_sendAllSims + s*numMems*numSpecs;
  float *t1 = t1AllSims + s*numMems;
  float *tau = tauAllSims + s*numMems;
  float *FS = FSAllSims + s*numMems;

  float *FSbkp = FSbkpAllSims + s*numMems;
  float *tauSSA = tauSSAAllSims + s*numMems;
  unsigned *flagSSA = flagSSAAllSims + s*numMems;
  unsigned *flagStep = flagStepAllSims + s*numMems;
  float *tauPreMin = tauPreMinAllSims + s*numMems;
  float *tau1 = tau1AllSims + s*numMems;
  float *tau2 = tau2AllSims + s* numMems;
  float *a0c = a0cAllSims + s*numMems;
  float *a0 = a0AllSims + s*numMems;

  float *Mbkp = MbkpAllSims + s*numMems*numSpecs;
  float *M = MAllSims + s*numMems*numSpecs;
  unsigned *HOR = HORAllSims + s*numMems*numSpecs;   

  float *auxVec1 = auxVec1AllSims + s*blockDim.x*numSpecs;	   
  float *auxVec2 = auxVec2AllSims + s*blockDim.x*numSpecs;	

  ////////////////////////////////// Elaboration:

  localStep = 0;
  localT = 0.0f;

  // ex getFS
  for (m=0; m<numMemsPerThread; m++) {
  	  //if(m*bDim + thIdx > numMems) break;
  	  FS[m*bDim + thIdx] =  membSize[m*bDim + thIdx];
  	  for (c=0; c<numSpecs; c++) 
  		  FS[m*bDim + thIdx] -= (M[c*numMems+m*bDim+thIdx] * molSize[c*numMems+m*bDim+thIdx]);
  }

  while (localT<stepMax) {


  	  /////////////////////////////////////////////////////////
  	  // ex dpp_step1

  	  for (m=0; m<numMemsPerThread; m++) {
  		  FSbkp[m*bDim + thIdx] = FS[m*bDim + thIdx];
  		  for (c=0; c<numSpecs; c++) {
  			  Mbkp[c*numMems+m*bDim+thIdx] = M[c*numMems+m*bDim+thIdx];
  			  if (Mfeed[c*numMems+m*bDim+thIdx]>0)
  				  M[c*numMems+m*bDim+thIdx] = Mfeed[c*numMems+m*bDim+thIdx];
  		  }
  	  }

  	  ////////////////////////////////////////////////////////
  	  // ex rule_max_app()

  	  for (m=0; m<numMemsPerThread; m++) {
  		  aux1 = 0.0f;
  		  for (r=0; r<numReacts; r++) {
  			  aux2 = FLT_MAX;
  			  for (c=0; c<numSpecs; c++) {
  				  if (leftSide[r*numMems*numSpecs + c*numMems+m*bDim+thIdx] != 0)
  					  aux1 = M[c*numMems+m*bDim+thIdx] / (float)leftSide[r*numMems*numSpecs + c*numMems+m*bDim+thIdx];
  				  if (aux1<aux2)
  					  aux2 = aux1;
  			  }
  			  critical[r*numMems+m*bDim+thIdx] = (aux2>0.0f && aux2<10.0f);
  		  }
  	  }

  	  ////////////////////////////////////////////////////////
  	  // ex rule_prob()

  	  // ex getFS
  	  for (m=0; m<numMemsPerThread; m++) {
  		  FS[m*bDim + thIdx] = membSize[m*bDim + thIdx];
  		  for (c=0; c<numSpecs; c++) 
  			  FS[m*bDim + thIdx] -= (M[c*numMems+m*bDim+thIdx] * molSize[c*numMems+m*bDim+thIdx]);
  	  }

  	  for (m=0; m<numMemsPerThread; m++) {
  		  for (r=0; r<numReacts; r++) {
  			  R_prob[r*numMems+m*bDim+thIdx] = Cvector[r*numMems+m*bDim+thIdx];
  			  if (tgtVector[r*numMems+m*bDim+thIdx] == -1) {
  			     R_prob[r*numMems+m*bDim+thIdx] /= FS[m*bDim + thIdx];
  			     R_prob[r*numMems+m*bDim+thIdx] /= membSize[m*bDim + thIdx];
  			  }
  			  for (c=0; c<numSpecs; c++) {
  				  switch (leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]) {
  					  case 0:
  						  break;
  					  case 1:
  						  R_prob[r*numMems+m*bDim+thIdx] *= M[c*numMems+m*bDim+thIdx];
  						  break;
  					  case 2:
  						  R_prob[r*numMems+m*bDim+thIdx] *= (M[c*numMems+m*bDim+thIdx] * (M[c*numMems+m*bDim+thIdx] - 1)) * 0.5f;
  						  break;
  					  default:
  						  for (e=1; e<=leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]; e++)
  							  R_prob[r*numMems+m*bDim+thIdx] *= (M[c*numMems+m*bDim+thIdx] - e + 1.0f) / e;
  						  break;
  				  }
  			  }
  		  }
  	  }

  	  ////////////////////////////////////////////////////////

  	  for (m=0; m<numMemsPerThread; m++) {
  		  a0[m*bDim + thIdx] = 0.0f;
  		  for (r=0; r<numReacts; r++)
  			   a0[m*bDim + thIdx] += R_prob[r*numMems+m*bDim+thIdx];
  	  }

  	  for (m=0; m<numMemsPerThread; m++) {

  		  if (flagSSA[m*bDim + thIdx] == 0) {
  			  if (a0[m*bDim + thIdx]>0.0f) {

  				  ////////////////////////////////////////////////////////
  				  // ex get_HOR()

  				  for (r=0; r<numReacts; r++) {
  					  order[r*numMems+m*bDim+thIdx] = 0;
  					  for (c=0; c<numSpecs; c++)
  						  order[r*numMems+m*bDim+thIdx] += leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  				  }

  				  for (c=0; c<numSpecs; c++) {
  					  HOR[c*numMems+m*bDim+thIdx] = 0;
  					  auxVec1[c*bDim+thIdx] = 0.0f;
  					  for (r=0; r<numReacts; r++) {
  					  
  						  aux2 = 0.0f;
  						  aux1 = 0.0f;
  						  if (leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]!=0) {
  							  switch (order[r*numMems+m*bDim+thIdx]) {
  							     case 1:
  								     aux2  = 1.0f;
  								     aux1 = 1.0f;
  								     break;
  							     case 2:
  								     if (leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] == 1) {
  									     aux2  = 2.0f;
  									     aux1 = 2.0f;
  								     } else {
  									     aux2  = 2.0f + (1.0f / (M[c*numMems+m*bDim+thIdx] - 1.0f));
  									     aux1 = 3.0f;
  								     }
  								     break;
  							     case 3:
  								     if (leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] == 1) {
  									     aux2  = 3.0f;
  									     aux1 = 4.0f;
  								     } else if (leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] == 2) {
  									     aux2  = 1.5f * (2.0f + (1.0f / (M[c*numMems+m*bDim+thIdx] - 1.0f)));
  									     aux1 = 5.0f;
  								     } else {
  									     aux2  = 3.0f + (1.0f / (M[c*numMems+m*bDim+thIdx] - 1.0f)) + (2.0f / (M[c*numMems+m*bDim+thIdx] - 2.0f));
  									     aux1 = 6.0f;
  								     }
  								     break;
  							     default: // shouldn't be reached
  								     break;
  							  }
  						  }
  						  if (aux1 > auxVec1[c*bDim+thIdx]) {
  							  auxVec1[c*bDim+thIdx] = aux1;
  							  HOR[c*numMems+m*bDim+thIdx] = (unsigned) aux2;
  						  }
  					  }
  				  }

  			  ////////////////////////////////////////////////////////
  			  // ex get_tau1()

  			  for (c=0; c<numSpecs; c++) {
  				  auxVec1[c*bDim+thIdx] = 0.0f;
  				  auxVec2[c*bDim+thIdx] = 0.0f;
  				  for (r=0; r<numReacts; r++) 
  					  if (critical[r*numMems+m*bDim+thIdx] == 0) {
  						  auxVec1[c*bDim+thIdx] +=	(float)var[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]    * R_prob[r*numMems+m*bDim+thIdx];
  						  auxVec2[c*bDim+thIdx] += powf((float)var[r*numMems*numSpecs + c*numMems + m*bDim + thIdx],2) * R_prob[r*numMems+m*bDim+thIdx];
  					  }
  			  }

  			  aux2 = FLT_MAX;
  			  tau1[m*bDim + thIdx] = FLT_MAX;
  			  for (c=0; c<numSpecs; c++) {
  				  //aux1 = fmaxf((0.03f*M[c*numMems+m*bDim+thIdx]) / HOR[c*numMems+m*bDim+thIdx], 1);
				  aux1 = fmaxf((eps*M[c*numMems+m*bDim+thIdx]) / HOR[c*numMems+m*bDim+thIdx], 1);
  				  if (auxVec1[c*bDim+thIdx]!=0 && auxVec2[c*bDim+thIdx]!=0.0f)
  					  aux2 = fminf(aux1/fabs(auxVec1[c*bDim+thIdx]), powf(aux1,2)/auxVec2[c*bDim+thIdx]);
  				  tau1[m*bDim + thIdx] = fminf(aux2, tau1[m*bDim + thIdx]);
  			  }

  			  for (r=0; r<numReacts; r++) 
  				  k_rule[r*numMems+m*bDim+thIdx] = 0;

  				  for (r=0; r<numReacts; r++)
  					  if (critical[r*numMems+m*bDim+thIdx]==1)
  						  aux1 += R_prob[r*numMems+m*bDim+thIdx];
  				  do rnd_t = curand_uniform(&rngStates[m*bDim + thIdx]);
  				  while (rnd_t==1.0f);
  				  tau2[m*bDim + thIdx] = logf(1.0f/rnd_t)/aux1;

  				  if (tau1[m*bDim + thIdx]<tau2[m*bDim + thIdx]) {
  					  flagStep[m*bDim + thIdx] = 2;
  					  tau[m*bDim + thIdx] = tau1[m*bDim + thIdx];
  				  } else {
  					  flagStep[m*bDim + thIdx] = 3;
  					  tau[m*bDim + thIdx] = tau2[m*bDim + thIdx];
  					  a0c[m*bDim + thIdx] = aux1;
  					  tauPreMin[m*bDim + thIdx] = tau2[m*bDim + thIdx];
  				  }

  			  ////////////////////////////////////////////////////////

  		   } else {
  		      flagStep[m*bDim + thIdx] = 4;
  		      tau[m*bDim + thIdx] = FLT_MAX;
  		   }
  	    } else {
  		  flagStep[m*bDim + thIdx] = 5;
  		  tau[m*bDim + thIdx] = tauSSA[m*bDim + thIdx];
  	    }

  	  }

  	  ////////////////////////////////////////////////////////////////////////////////////

  	  __syncthreads();
  	  
  	  splitShared = FLT_MAX;	  
  	  for (m=0; m<numMems; m++)
  		  if (tau[m]<splitShared)
  			  splitShared = tau[m];

  	  for (m=0; m<numMemsPerThread; m++)
  		  tau[m*bDim + thIdx] = splitShared;

  	  ////////////////////////////////////////////////////////////////////////////////////
  	  // dpp_step2

  	  ////////////////////////////////////////////////////////
  	  // ex get_tau2()

  	  for (m=0; m<numMemsPerThread; m++) {
  		  
  		  for (c=0; c<numSpecs; c++)
  			  M_send[c*numMems+m*bDim+thIdx] = 0.0f;

  		  if (flagStep[m*bDim+thIdx]==2)
  			  for (r=0; r<numReacts; r++)
  				  if (critical[r*numMems+m*bDim+thIdx] == 0.0f)
  					  k_rule[r*numMems+m*bDim+thIdx] = curand_poisson(&rngStates[m*bDim+thIdx], tau[m*bDim+thIdx] * R_prob[r*numMems+m*bDim+thIdx]);
  		    
  		  if (flagStep[m*bDim+thIdx]==3) {
  			  for (r=0; r<numReacts; r++) {
  				  if (critical[r*numMems+m*bDim+thIdx] == 0.0f) {
  					  k_rule[r*numMems+m*bDim+thIdx] = curand_poisson(&rngStates[m*bDim+thIdx], tau[m*bDim+thIdx] * R_prob[r*numMems+m*bDim+thIdx]);
  					  R_probC[r*numMems+m*bDim+thIdx] = 0.0f;
  				  } else
  					  R_probC[r*numMems+m*bDim+thIdx] = R_prob[r*numMems+m*bDim+thIdx];
  			  }
  		  
  			  r=0;
  			  if (tauPreMin[m*bDim+thIdx]==tau[m*bDim+thIdx]) {
  				  do rnd_t = curand_uniform(&rngStates[m*bDim+thIdx]);
  				  while (rnd_t==1.0f);
  				  rnd_t *= a0c[m*bDim+thIdx];
  				  alpha = R_probC[r*numMems+m*bDim+thIdx];
  				  while (alpha<rnd_t) {
  					  r++;
  					  alpha += R_probC[r*numMems+m*bDim+thIdx];
  				  }
  				  k_rule[r*numMems+m*bDim+thIdx] = 1;
  			  }
  		  }

  		  ////////////////////////////////////////////////////////

  		  if (flagStep[m*bDim+thIdx]==2 || flagStep[m*bDim+thIdx]==3 || (tauSSA[m*bDim+thIdx] == tau[m*bDim+thIdx] && (flagStep[m*bDim+thIdx] == 1 || flagStep[m*bDim+thIdx] == 5)) ) {
  			  flagSSA[m*bDim+thIdx] = 0;
  			  if (tgtKind == 1) {
  				  for (r=0; r<numReacts; r++) {
  					  if (k_rule[r*numMems+m*bDim+thIdx]>0) {
  						  if (tgtVector[r*numMems+m*bDim+thIdx] == -1) // execute an internal rule
  							  for (c=0; c<numSpecs; c++)
  								  M[c*numMems+m*bDim+thIdx] += (float)var[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] * (float)k_rule[r*numMems+m*bDim+thIdx];
  						  else if (tgtVector[r*numMems+m*bDim+thIdx] == -2) { // execute a nondeterministic communication rule
  							  for (c=0; c<numSpecs; c++)
  								  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  							  for (e=0; e<k_rule[r*numMems+m*bDim+thIdx]; e++) {
  								  for (c=0; c<numSpecs; c++)
  									  atomicAdd(M_send + c*numMems + (int)ceilf((curand_uniform(&rngStates[m*bDim+thIdx])*numMems)), (float)rightSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]);
  							  }
  						  } else // execute a deterministic communication rule
  							  for (c=0; c<numSpecs; c++) {
  								  atomicAdd(M_send + c*numMems + tgtVector[r*numMems+m*bDim+thIdx], (float)k_rule[r*numMems+m*bDim+thIdx] * (float)rightSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]);
  								  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  							  }
  					  }
  				  }
  			  } else {
  				  for (r=0; r<numReacts; r++)
  					  if (k_rule[r*numMems+m*bDim+thIdx]>0)
  						  for (c=0; c<numSpecs; c++)
  							  atomicAdd(M_send + c*numMems + tgtMatrix[r*numMems*numSpecs + c*numMems + m*bDim + thIdx], (float)k_rule[r*numMems+m*bDim+thIdx] * (float)varSend[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] );
  				  // update the internal state of the process
  				  for (r=0; r<numReacts; r++)
  					  for (c=0; c<numSpecs; c++)
  						  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  			  }
  		  }

  		  if ( tauSSA[m*bDim+thIdx] != tau[m*bDim+thIdx] && (flagStep[m*bDim+thIdx] == 1 || flagStep[m*bDim+thIdx] == 5) )
  			  tauSSA[m*bDim+thIdx] -= tau[m*bDim+thIdx];

  	  }

  	  ////////////////////////////////////////////////////////////////////////////////////

  	  __syncthreads();

  	  // ex receive
  	  for (m=0; m<numMemsPerThread; m++)
  		  for (c=0; c<numSpecs; c++)
  			  M[c*numMems+m*bDim+thIdx] += M_send[c*numMems+m*bDim+thIdx];

  	  do {
  			  
  		  // ex getFS
  		  for (m=0; m<numMemsPerThread; m++) {
  			  FS[m*bDim + thIdx] = membSize[m*bDim + thIdx];
  			  for (c=0; c<numSpecs; c++) 
  				  FS[m*bDim + thIdx] -= (M[c*numMems+m*bDim+thIdx] * molSize[c*numMems+m*bDim+thIdx]);
  		  }

  		  __syncthreads();
  			  
  		  ////////////////////////////////////////////////////////
  		  // ex checkNegFS()

  		  splitShared = FLT_MAX;
  		  for (m=0; m<numMems; m++)
  			  if (FS[m]<splitShared)
  				  splitShared = FS[m];

  		  splitShared = (splitShared<0);
  	  
  		  ////////////////////////////////////////////////////////
  		  
  		  if (splitShared) {

  			  ////////////////////////////////////////////////////////////////////////////////////
  			  // dpp_step3

  			  for (m=0; m<numMemsPerThread; m++) {
  				  
  				  FS[m*bDim + thIdx] = FSbkp[m*bDim + thIdx];
  				  for (c=0; c<numSpecs; c++) {
  					  M[c*numMems+m*bDim+thIdx] = Mbkp[c*numMems+m*bDim+thIdx];
  					  M_send[c*numMems+m*bDim+thIdx] = 0.0f;
  				  }
  				  for (r=0; r<numReacts; r++) 
  					  k_rule[r*numMems+m*bDim+thIdx] = 0;

  				  tau[m*bDim + thIdx] /= 2.0f;
  				  tauSSA[m*bDim + thIdx] += tau[m*bDim + thIdx];

  				  if (flagStep[m*bDim + thIdx]==2 || flagStep[m*bDim + thIdx]==3)
  					  for(r=0; r<numReacts; r++)
  						  if (critical[r*numMems+m*bDim+thIdx] == 0.0f)
  							  k_rule[r*numMems+m*bDim+thIdx] = curand_poisson(&rngStates[m*bDim + thIdx], tau[m*bDim + thIdx] * R_prob[r*numMems+m*bDim+thIdx]);

  				  if (tgtKind == 1) {
  					  for (r=0; r<numReacts; r++) {
  						  if (k_rule[r*numMems+m*bDim+thIdx]>0) {
  							  if (tgtVector[r*numMems+m*bDim+thIdx] == -1) // execute an internal rule
  								  for (c=0; c<numSpecs; c++)
  									  M[c*numMems+m*bDim+thIdx] += (float)var[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] * (float)k_rule[r*numMems+m*bDim+thIdx];
  							  else if (tgtVector[r*numMems+m*bDim+thIdx] == -2) { // execute a nondeterministic communication rule
  								  for (c=0; c<numSpecs; c++)
  									  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  								  for (e=0; e<k_rule[r*numMems+m*bDim+thIdx]; e++) {
  									  for (c=0; c<numSpecs; c++)
  										  atomicAdd(M_send + c*numMems + (int)ceilf((curand_uniform(&rngStates[m*bDim+thIdx])*numMems)), (float)rightSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]);
  								  }
  							  } else // execute a deterministic communication rule
  								  for (c=0; c<numSpecs; c++) {
  									  atomicAdd(M_send + c*numMems + tgtVector[r*numMems+m*bDim+thIdx], (float)k_rule[r*numMems+m*bDim+thIdx] * (float)rightSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx]);
  									  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  								  }
  						  }
  					  }
  				  } else {
  					  for (r=0; r<numReacts; r++)
  						  if (k_rule[r*numMems+m*bDim+thIdx]>0)
  							  for (c=0; c<numSpecs; c++)
  								  atomicAdd(M_send + c*numMems + tgtMatrix[r*numMems*numSpecs + c*numMems + m*bDim + thIdx], (float)k_rule[r*numMems+m*bDim+thIdx] * (float)varSend[r*numMems*numSpecs + c*numMems + m*bDim + thIdx] );
  					  // update the internal state of the process
  					  for (r=0; r<numReacts; r++)
  						  for (c=0; c<numSpecs; c++)
  							  M[c*numMems+m*bDim+thIdx] -= (float)k_rule[r*numMems+m*bDim+thIdx] * (float)leftSide[r*numMems*numSpecs + c*numMems + m*bDim + thIdx];
  				  }
  			  }

  			  ////////////////////////////////////////////////////////////////////////////////////

  			  __syncthreads();

  			  
  			  // ex receive
  			  for (m=0; m<numMemsPerThread; m++)
  				  for (c=0; c<numSpecs; c++)
  					  M[c*numMems+m*bDim+thIdx] += M_send[c*numMems+m*bDim+thIdx];
  		  }

  	  } while (splitShared);

  	  localStep++;

  	  ////////////////////////////////////////////////////////////////////////////////////
  	  // dpp_step4

  	  for (m=0; m<numMemsPerThread; m++) {
  		  
  		  if (tau[m*bDim + thIdx] == FLT_MAX)
  			  tau[m*bDim + thIdx] = stepMax - localT + 1.0f;
  		  t1[m*bDim + thIdx] = localT + tau[m*bDim + thIdx];

  		  ////////////////////////////////////////////////////////////////////////////////////
  		  if (buffBool) {
  			  if (remainderf(localStep,buffEvery)==0.0f && buffRowCounter[m*bDim+thIdx]<buffRows) {
  				  
  				  buffM[buffRowCounter[m*bDim+thIdx]*numMems*buffCols + 0*numMems + m*bDim+thIdx] = t1[m*bDim + thIdx];
  				  for (c=0; c<numSpecs; c++)
  					  memcpy(buffM + buffRowCounter[m*bDim+thIdx]*numMems*buffCols + (c+1)*numMems + m*bDim, M + c*numMems + m*bDim, bDim * sizeof(float));
  				  buffRowCounter[m*bDim+thIdx]++;

  			  } else
  				  printf("buffer full!\n");
  		  }

  	  }

  	  __syncthreads();

  	  splitShared = 0;
  	  for (m=0; m<numMems; m++)
  			  if (t1[m]>splitShared)
  				  splitShared = t1[m];
      localT = splitShared;

  }

  step[s] = localStep;
  t[s] = localT;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// HOST code /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void CUDAcall(cudaError_t cudaResult, std::string variable) {

if (cudaResult != cudaSuccess) {
	variable += ": ";
	variable += cudaGetErrorString(cudaResult);
	throw std::runtime_error(variable);
}
}

/* Count the number of rows and columns of a binary file */
void uliRows_uliColumns_old(FILE *fPtr, unsigned *uliRows_h_Ptr, unsigned *uliCols_h_Ptr){
  long e;
  unsigned c=0, r=0;

  while ((e=getc(fPtr))!=EOF) {
  	if (e == '\t') c++;
  	if (e == '\n') break;
  }
  while ((e=getc(fPtr))!=EOF)
	if (e == '\n') r++;

  rewind(fPtr);
  *uliRows_h_Ptr = r + 1;
  *uliCols_h_Ptr = c; // + 1; vedi sotto
}

/* Count the number of rows and columns of a binary file */
void uliRows_uliColumns(FILE *fPtr, unsigned *uliRows_h_Ptr, unsigned *uliCols_h_Ptr){
  long e, pos;
  unsigned c=0, r=0;

  while ((e=getc(fPtr)) != '\n') ; //printf("%c",(char)c) ; //read #num\n

  pos = ftell(fPtr);
    
  while ((e=getc(fPtr))!=EOF) {
	if (e == '\t') c++;
  	if (e == '\n') break;
  }

  while (( (e=getc(fPtr))!=EOF) && (e != '#') )
  	if (e == '\n') r++;

  fseek(fPtr, pos, SEEK_SET);
  *uliRows_h_Ptr = r + 1;
  *uliCols_h_Ptr = c;// + 1; occhio che dopo l'ultimo elemento ci sia \t, altrimenti va messo +1 + 1;
}



/* Allocate and read an unsigned table from a file as a 1D vector */
unsigned* uLoad_table(FILE *fPtr, unsigned *uliRows_Ptr, unsigned *uliCols_Ptr) {
  unsigned j=0, i=0;
  unsigned *m;

  uliRows_uliColumns(fPtr, uliRows_Ptr, uliCols_Ptr);

  m = (unsigned *)malloc((*uliRows_Ptr)*(*uliCols_Ptr)*sizeof(unsigned) );
  for(j=0; j<*uliRows_Ptr; j++){
  	  for(i=0; i<*uliCols_Ptr-1; i++)
	  	fscanf(fPtr,"%i\t", &m[j*(*uliCols_Ptr)+i]);
  	  i = *uliCols_Ptr-1;
  	  fscanf(fPtr,"%i\n", &m[j*(*uliCols_Ptr)+i]);
  }
  return m;
}


/* Allocate and read an int table from a file as a 1D vector 
int* iLoad_table(FILE *fPtr, int mId) {
  int *m;
  unsigned  j=0, i=0;

  m = (int *) malloc( numReacts*numSpecs*sizeof(int) );
  for (j=0; j<numReacts; j++) {
	  for (i=0; i<numSpecs-1; i++)
		  fscanf(fPtr,"%d\t", &m[j*numSpecs+i]);
	  i=numSpecs-1;
	  fscanf(fPtr,"%d\n", &m[j*numSpecs+i]);
  }
  return m;
}
*/

/* Write the output matrix in the binary MULTI output file */
void write_matrix_bin(int m, int s) {
  
  if (buffBool_h) {
	  printf("sim %d mem %d: scrivo %ld righe (%ld bytes)\n", s, m, buffRowCounterAllSims_h[s*numMems + m], buffCols_h*buffRowCounterAllSims_h[s*numMems + m]*sizeof(float));
	  fwrite(&m, sizeof(unsigned),1,buffer_fPtr[s]);
	  for (unsigned c=0; c<buffCols_h; c++)
		fwrite(buffMAllSims_h + s*buffRows_h*buffCols_h*numMems + buffRowCounterAllSims_h[s*numMems + m]*buffCols_h*numMems + c*numMems + m, sizeof(float), 1, buffer_fPtr[s]);
	  buffRowCounterAllSims_h[s*numMems + m] = 0;
  } else {
	  fwrite(&m, sizeof(unsigned),1,buffer_fPtr[s]);
	  for (unsigned c=0; c<numSpecs; c++)
		  fwrite(MAllSims_h + s*numSpecs*numMems + c*numMems + m, sizeof(float), 1, buffer_fPtr[s]);
  }

}
/* Write the output matrix in the txt MULTI output file */
void write_matrix_txt(int m, int s) {
  
	if (buffBool_h) {
	  printf("sim %d mem %d: writing %d lines\n", s, m, buffRowCounterAllSims_h[s*numMems + m]);
	  for (unsigned r=0; r<buffRowCounterAllSims_h[s*numMems + m]; r++) {
		  fprintf(buffer_fPtr[s],"mem%4d: %4G\t", m, buffMAllSims_h[s*buffRows_h*buffCols_h*numMems + buffRowCounterAllSims_h[s*numMems + m]*buffCols_h*numMems + 0*numMems + m]);
		  for (unsigned c=0; c<numSpecs; c++)
			  fprintf(buffer_fPtr[s],"%4G\t", buffMAllSims_h[s*buffRows_h*buffCols_h*numMems + buffRowCounterAllSims_h[s*numMems + m]*buffCols_h*numMems + (c+1)*numMems + m]);
		  fprintf(buffer_fPtr[s],"\n");
	  }
	  buffRowCounterAllSims_h[s*numMems + m] = 0;
  } else {
	  fprintf(buffer_fPtr[s],"mem%4d: %4G\t", m, MAllSims_h[s*numSpecs*numMems + 0*numMems + m]);
	  for (unsigned c=1; c<numSpecs; c++)
		  fprintf(buffer_fPtr[s],"%4G\t", MAllSims_h[s*numSpecs*numMems + c*numMems + m]);
  }

  fprintf(buffer_fPtr[s],"\n");
}

// Read the first part of input files: LEFT/RIGHT_SIDE, TGT_VECTOR/MATRIX, INDEXES
void load_files1(int mId, unsigned** LS, unsigned** RS, int** TM, int** TV, FILE *flmatrix, FILE *frmatrix, FILE *ftgtmatrix, FILE *findexes) {
//  char input1[100], input2[100], app[20], file[60];
  unsigned int i, j;

  //sprintf(app,"%d", mId);

  unsigned rowsL, rowsR, colsL, colsR;

  fprintf(log_fPtr, "\n---> MEMBRANE: %d\n\n", mId);

  fprintf(log_fPtr,"Loaded stepMax:\n%g\n",stepMax_h);

  // Read the LEFT_SIDE input file with the mId or the common input file
/*
  sprintf(input1, "input/");
  sprintf(input2, "input/");
  strcpy(file,"left_side_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"left_side");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
	if ((fPtr = fopen (input2, "r")) == NULL)
		fprintf (log_fPtr, "\nCannot open left_side!\n");
*/

  LS[mId] = uLoad_table(flmatrix, &rowsL, &colsL);
//  fclose (fPtr);

  if (numSpecs != colsL)  fprintf(log_fPtr, "left cols %u != numSpecs %u --> using numSpecs\n", colsL, numSpecs);
  if (numReacts != rowsL) fprintf(log_fPtr, "left rows %u != numReacts %u --> using numReacts\n", rowsL, numReacts);

  // Print the LEFT HAND SIDE into the LOG file
  fprintf(log_fPtr,"\nLoaded left_side:\n");
  for (i=0; i<numReacts; i++) {
	for (j=0; j<numSpecs; j++)
		fprintf(log_fPtr,"%u\t",LS[mId][i*numSpecs+j]);
	fprintf(log_fPtr,"\n");
  }

  // Read the RIGHT_SIDE input file with the mId or the common input file
/*
  sprintf(input1, "input/");
  sprintf(input2, "input/");
  strcpy(file,"right_side_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"right_side");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
	if ((fPtr = fopen (input2, "r")) == NULL)
		fprintf (log_fPtr, "\nCannot open right_side!\n");
*/	
  RS[mId] = uLoad_table(frmatrix, &rowsR, &colsR);
//  fclose (fPtr);

  if (numSpecs != colsR)  fprintf(log_fPtr, "right cols %u != numSpecs %u --> using numSpecs\n", colsR, numSpecs);
  if (numReacts != rowsR) fprintf(log_fPtr, "right rows %u != numReacts %u --> using numReacts\n", rowsR, numReacts);

  // Print the RIGHT HAND SIDE into the LOG file
  fprintf(log_fPtr,"\nLoaded right_side:\n");
  for (i=0; i<numReacts; i++){
	  for (j=0; j<numSpecs; j++)
		  fprintf(log_fPtr,"%u\t",RS[mId][i*numSpecs+j]);
	  fprintf(log_fPtr,"\n");
  }

  // Read the TGT_VECTOR/TGT_MATRIX input file with the mId or the common input file
/*
  sprintf(input1, "input", numMems);
  sprintf(input2, "input", numMems);
  strcpy(file,"tgt_vector_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"tgt_vector");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
  	  if ((fPtr = fopen (input2, "r")) == NULL)
  		  fprintf (log_fPtr,"\nCannot open tgt_vector!\n");
*/


 	 
//  if (fPtr != NULL) {
  TV[mId] = (int *) malloc(numReacts * sizeof(int));
  tgtKind_h = 1;
  // Print the TGT_VECTOR into the LOG file

  fscanf(ftgtmatrix,"%*s\n"); //Read #num\n
  fprintf(log_fPtr,"\nLoaded tgt_vector:\n");

  for (j=0; j<numReacts-1; j++) {
  	fscanf(ftgtmatrix,"%d\t", &TV[mId][j]);
	fprintf(log_fPtr,"%d\t",TV[mId][j]);
  }

  fscanf(ftgtmatrix,"%d\n", &TV[mId][j]);
  fprintf(log_fPtr,"%d\n\n",TV[mId][j]);
  
/*
  	fclose (fPtr);
  } else {
	sprintf(input1, "input%d/", numMems);
    sprintf(input2, "input%d/", numMems);
    strcpy(file,"tgt_matrix_");
    strcat(file,app);
    strcat(input1,file);
    strcpy(file,"tgt_matrix");
    strcat(input2,file);
	if ((fPtr = fopen (input1, "r")) == NULL)
		if((fPtr = fopen (input2, "r")) == NULL)
  	        	fprintf (log_fPtr,"\nCannot open tgt_matrix\n");
  	        
  	if (fPtr != NULL) {
		TM[mId] = iLoad_table(fPtr,mId);
  	    fclose (fPtr);
  	    // Print the TGT_MATRIX in the LOG file
  	    fprintf(log_fPtr,"\nLoaded tgt_matrix:\n");
  	    for (i=0; i<numReacts; i++) {
  	    	for (j=0; j<numSpecs; j++)
  	    		fprintf(log_fPtr,"%d\t",TM[mId][i*numSpecs+j]);
	        fprintf(log_fPtr,"\n");
  	   }
	}
  	fprintf(log_fPtr,"\n");
  }
*/
  // Read the INDEXES input file with the mId or the common input file
  unsigned tmp;
  int indexes;
/*
  sprintf(input1, "input%d/", numMems);
  sprintf(input2, "input%d/", numMems);
  strcpy(file,"indexes_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"indexes");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
  	  if ((fPtr = fopen (input2, "r")) == NULL)
  		  fprintf (log_fPtr,"\nCannot open indexes!\n");
  		
  if(fPtr != NULL) {
*/

  // Print the indexes in the LOG file
  fprintf(log_fPtr,"\nLoaded indexes:\t");
  uliRows_uliColumns_old(findexes, &tmp, &buffCols_h);
  fprintf(log_fPtr,"\n-- %d --\n",buffCols_h);
  
  for (i=0; i<buffCols_h; i++) {
  	  fscanf(findexes, "%d\t", &indexes);
  	  fprintf(log_fPtr,"%d\t", indexes);
  }
  fscanf(findexes, "%d\n", &indexes);
  fprintf(log_fPtr,"%d\n\n", indexes);

  rewind(findexes);// THIS IS A SINGLE FILE, otherwise this part of the code has to be modified

  buffCols_h++; //perche'???

  fflush(log_fPtr);
//  }
}


// Read the second part of input files: C_VECTOR, M_0, M_FEED, SIZES
void load_files2(int mId, FILE * fm0, FILE *fcmatrix, FILE * fmfeed, FILE * fsize) {
//  FILE *fPtr;
  long double temp;
//  char input1[100], input2[100], app[20], file[60];
  unsigned int i,j;
/*
  sprintf(app,"%d", mId);

  // Read the C_VECTOR input file with the mId or the common input file
  sprintf(input1, "input%d/", numMems);
  sprintf(input2, "input%d/", numMems);
  strcpy(file,"c_vector_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"c_vector");
  strcat(input2,file);
*/
  fprintf(log_fPtr, "\n---> MEMBRANE: %d\n\n", mId);
/*
  if ((fPtr = fopen (input1, "r")) == NULL)
  	  if((fPtr = fopen (input2, "r")) == NULL)
  		  fprintf (log_fPtr,"\nCannot open c_vector!\n");
  	  
  if(fPtr != NULL) {
*/  

  fprintf(log_fPtr,"\nLoaded c_vector:\n");
  fscanf(fcmatrix,"%*s\n"); //Read #num\n
  
  for(j=0; j<numReacts-1; j++) {
	fscanf(fcmatrix, "%LG\t", &temp);
  	Cvector_h[j*numMems+mId] = (float) temp;
  	// Print the CONSTANTS VECTOR in the LOG file
  	fprintf(log_fPtr,"%g\t", Cvector_h[j*numMems+mId]);
  }

  fscanf(fcmatrix, "%LG\n", &temp);
  Cvector_h[j*numMems+mId] = (float) temp;
  // Print the CONSTANTS VECTOR in the LOG file
  fprintf(log_fPtr,"%g\n", Cvector_h[j*numMems+mId]);


//  }

  // Read the M_0 input file with the mId or the common input file

  unsigned tmp, colsM0;
/*
  sprintf(input1, "input%d/", numMems);
  sprintf(input2, "input%d/", numMems);
  strcpy(file,"M_0_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"M_0");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
  	  if((fPtr = fopen (input2, "r")) == NULL)
  		  fprintf (log_fPtr,"\nCannot open M_0!\n");

  if (fPtr != NULL) {
*/
  uliRows_uliColumns(fm0, &tmp, &colsM0);
  if(numSpecs != colsM0) printf("numSpecs %d != cols in M_0 %d --> using numSpecs\n", numSpecs, colsM0);
  fprintf(log_fPtr,"\nLoaded M_0:\n");
  for(i=0; i<numSpecs-1; i++) {
	fscanf(fm0, "%LG\t", &temp);
	MAllSims_h[0*numMems*numSpecs + i*numMems + mId] = (float) temp;
	for (unsigned s=0; s<numSims; s++) MAllSims_h[s*numMems*numSpecs + i*numMems + mId] = MAllSims_h[0*numMems*numSpecs + i*numMems + mId];
	// Print the MULTISET vector in the LOG file
	fprintf(log_fPtr,"%g\t", MAllSims_h[0*numMems*numSpecs + i*numMems + mId]);
  }
  fscanf(fm0, "%LG\n", &temp);
  MAllSims_h[0*numMems*numSpecs + i*numMems + mId] = (float) temp;
  for (unsigned s=0; s<numSims; s++) MAllSims_h[s*numMems*numSpecs + i*numMems + mId] = MAllSims_h[0*numMems*numSpecs + i*numMems + mId];
  fprintf(log_fPtr,"%g\n", MAllSims_h[0*numMems*numSpecs + i*numMems + mId]);
	  
	  
/*  fclose (fPtr);
  }
*/
  // Read the M_FEED input file with the mId or the common input file

/*
  sprintf(input1, "input%d/", numMems);
  sprintf(input2, "input%d/", numMems);
  strcpy(file,"M_feed_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"Mfeed");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
  if((fPtr = fopen (input2, "r")) == NULL)
  	  fprintf (log_fPtr,"\nCannot open Mfeed!\n");

  if(fPtr != NULL) {
*/

  fprintf(log_fPtr,"\nLoaded Mfeed:\n");
  fscanf(fmfeed,"%*s\n"); //Read #num\n
    
  for (i=0; i<numSpecs-1; i++) {
	  fscanf(fmfeed, "%LG\t", &temp);
	  Mfeed_h[i*numMems+mId] = (float) temp;
	  // Print the FEEDING MULTISET VECTOR in the LOG file
	  fprintf(log_fPtr,"%g\t", Mfeed_h[i*numMems+mId]);
  }
  fscanf(fmfeed, "%LG", &temp);
  Mfeed_h[i*numMems+mId] = (float) temp;
  // Print the FEEDING MULTISET VECTOR in the LOG file
  fprintf(log_fPtr,"%g\n\n", Mfeed_h[i*numMems+mId]);

/*	  
  	  fclose (fPtr);
  }
*/
  // Read the SIZES input file with the mId or the common input file
/*
  sprintf(input1, "input%d/", numMems);
  sprintf(input2, "input%d/", numMems);
  strcpy(file,"sizes_");
  strcat(file,app);
  strcat(input1,file);
  strcpy(file,"sizes");
  strcat(input2,file);
  if ((fPtr = fopen (input1, "r")) == NULL)
  	  if((fPtr = fopen (input2, "r")) == NULL)
  		  fprintf (log_fPtr,"\nCannot open sizes.txt\n");

  if(fPtr != NULL) {
*/
  fprintf(log_fPtr,"\nLoaded molSize:\n");
  fscanf(fsize,"%*s\n"); //Read #num\n  
  
  for(i=0; i<numSpecs-1; i++) {
	fscanf(fsize, "%LG\t", &temp);
	molSize_h[i*numMems+mId] = (float) temp;
	// Print the MOLECULAR SIZES in the LOG file
	fprintf(log_fPtr,"%g\t", molSize_h[i*numMems+mId]);
  }
  fscanf(fsize, "%LG\n", &temp);
  molSize_h[i*numMems+mId] = (float) temp;
  // Print the MOLECULAR SIZES in the LOG file
  fprintf(log_fPtr,"%g\n\n", molSize_h[i*numMems+mId]);
  
//  fprintf(log_fPtr,"\n");

  fprintf(log_fPtr,"\nLoaded membSize:\n");
  fscanf(fsize, "%LG\n", &temp);
  membSize_h[mId] = (float) temp;
  fprintf(log_fPtr,"%g\n",  membSize_h[mId]);
/*
  fclose (fPtr);
  }
*/
  fflush(log_fPtr);
}


// Set HOST data reading from input files
void CPUinitializeData(bool bBufferOut) {
  FILE *flmatrix, *frmatrix, *ftgtmatrix, *findexes;
  FILE *fm0, *fcmatrix, *fmfeed, *fsize;

  char file[100], input[100];

  printf("Number of simulations: %d\n", numSims);
  if (numSims < numSimsMin || numSims > numSimsMax)
     printf("\nSpecified number of simulations (%d) is invalid, must be between %d and %d.\n", numSims, numSimsMin, numSimsMax);

  buffBool_h = bBufferOut;

  printf("\nSeeds:\n");
  seed_h = (unsigned long *) malloc(numSims*sizeof(unsigned long));
  for (unsigned s=0; s<numSims; s++) {
  	  seed_h[s] = s + (unsigned long)time(NULL);
  	  printf("seed[%4d] = %lu\n", s, seed_h[s]);
      // Check requested seed is valid
      if (seed_h[s] == 0)
  	  printf("\nSpecified seed[%d] is invalid, must be non-zero.\n", s);
  }

  /* Read data common to all membranes */
  FILE *fPtr;
  long double temp;

  // Read time_max
  sprintf(input, "input/");
  strcpy(file,"time_max");
  strcat(input,file);
  if ((fPtr = fopen (input,"r")) == NULL)
      printf ("\nCannot open time_max\n");
  else {
      fscanf(fPtr,"%LG",&temp);
      stepMax_h = (float) temp;
      printf("\nTime max: %g\n", stepMax_h);
      fclose(fPtr);
  }

  // Read every
  sprintf(input, "input/");
  strcpy(file,"every");
  strcat(input,file);
  if ((fPtr = fopen (input, "r")) == NULL) {
  	  buffEvery_h = 1;
          buffRows_h = 10000;
  	  printf ("\nCannot open every - loading default %d and %d\n", buffEvery_h, buffRows_h);
  } else {
	fscanf(fPtr,"%d", &buffEvery_h);
      	fscanf(fPtr,"%d", &buffRows_h);
      	if (buffBool_h) {
  		printf ("\nOutput saved every: %d iteration(s)", buffEvery_h);
		printf("\nNumber of buffer lines: %d\n", buffRows_h);
      	} else  printf("\nOutput saved at final iteration\n");
	fclose (fPtr);
  }

  sprintf(input, "input/");
  strcpy(file,"eps");
  strcat(input,file);
  // Read eps
  if ((fPtr = fopen (input,"r")) == NULL)	  {
  	printf ("\nCannot open eps\n");
	eps_h = 0.03f;
  }
  else {
  	fscanf(fPtr,"%LG",&temp);
	eps_h = (float) temp;
  	printf("EPS: %g ",eps_h);
	fclose(fPtr);
  }


  // Read numMembranes
  unsigned numMemsFromFile;
  sprintf(input, "input/");
  strcpy(file,"numMembranes.txt");
  strcat(input,file);
  if ((fPtr = fopen (input,"r")) == NULL)
  	  printf ("\nCannot obtain the numbers of membranes, the file numMembranes.txt is required\n");
  else {
  	  fscanf(fPtr,"%d",&numMemsFromFile);
  	  fclose(fPtr);
  	  if (numMems != numMemsFromFile) printf("\nDefined number of membranes differs from that read from numMebranes\n");
  }

  /* HOST allocation (1st part) */

  // Output files

  if (buffBool_h) {
	buffRowCounterAllSims_h = (unsigned*) malloc(numSims * numMems * sizeof(unsigned));
	for (unsigned s=0; s<numSims; s++)
        	for (unsigned m=0; m<numMems; m++)
        		buffRowCounterAllSims_h[s*numMems + m] = 0;
  }

  buffer_fPtr	  = (FILE **)  malloc(numSims * sizeof(FILE*));

  // Temporary constants

  unsigned **unLeft_side2D_h, **unRight_side2D_h;
  int **iTgt_matrix2D_h, **iTgt_vector2D_h;
  unLeft_side2D_h  = (unsigned**) malloc(numMems * sizeof(void*));
  unRight_side2D_h = (unsigned**) malloc(numMems * sizeof(void*));
  iTgt_matrix2D_h  = (int **)	  malloc(numMems * sizeof(void*));
  iTgt_vector2D_h  = (int **)	  malloc(numMems * sizeof(void*));

  // Costants

  Mfeed_h	  = (float *) malloc(numMems*numSpecs * sizeof(float));
  molSize_h		  = (float *) malloc(numMems*numSpecs * sizeof(float));

  membSize_h			  = (float *) malloc(numMems * sizeof(float));

  // Variables

  t_h		  = (float *)	 malloc(numSims * sizeof(float));
  step_h	  = (unsigned *) malloc(numSims * sizeof(unsigned));

  MAllSims_h  = (float *)    malloc(numSims * numMems * numSpecs * sizeof(float));

  /* HOST assignments (1st part) */

  tgtKind_h = 0;

  char mem_sim[20];

  // Create a LOG file
  sprintf(file,"output%u_%u_%u/log.txt", numMems, threadsPerBlock, numSims);
  if ((log_fPtr = fopen (file,"w")) == NULL) {
  	  printf ("\nCannot write in the current directory - Log file\n");
  	  exit(0);
  }



// Read the first part of input files: LEFT/RIGHT_SIDE, TGT_VECTOR/MATRIX, INDEXES

  if((flmatrix = fopen ("input//L_matrix", "r")) == NULL){
	printf("Cannot open the left side\n");
	exit(0);
  }

  if((frmatrix = fopen ("input/R_matrix", "r")) == NULL){
	printf("Cannot open the right side\n");
	exit(0);
  }
  

  if((ftgtmatrix = fopen ("input/tgt_matrix", "r")) == NULL){
	printf("Cannot open the tgt matrix\n");
	exit(0);
  }


  if((findexes = fopen ("input/indexes", "r")) == NULL){
	printf("Cannot open indexes\n");
	exit(0);
  }

  /* For each membrane */
  for (unsigned m=0; m<numMems; m++) {
  	  // Read the first part of input files
  	  load_files1(m,unLeft_side2D_h, unRight_side2D_h, iTgt_matrix2D_h, iTgt_vector2D_h, flmatrix, frmatrix, ftgtmatrix, findexes);
  }


  fclose(flmatrix); fclose(frmatrix); fclose(ftgtmatrix); fclose(findexes);

  printf("Data 1 Acquired\n"); fflush(NULL);

  /* For each simulation */
  for (unsigned s=0; s<numSims; s++) {
  	  sprintf(mem_sim, "s%d.txt", s);
  		  
  	  // create a MULTI file labelled with the sId
  	  sprintf(file,"output%u_%u_%u/multi_",numMems, threadsPerBlock, numSims);
  	  strcat(file,mem_sim);
  	  if ((buffer_fPtr[s] = fopen (file,"w")) == NULL)
  		  printf ("\nCannot write in the current directory - Output file for sim = %u \n", s);
  		  
  	  // write on the MULTI file
  	  //fwrite(&buffCols,sizeof(unsigned),1,buffer_fPtr[s*numMems + m]);
  	  fprintf(buffer_fPtr[s], "%u\n", buffCols_h);
  }

  /* HOST allocation (2nd part) */

  /* Constants */
  leftSide_h	  = (unsigned*) malloc(numMems*numSpecs*numReacts * sizeof(unsigned));
  rightSide_h	  = (unsigned*) malloc(numMems*numSpecs*numReacts * sizeof(unsigned));
  var_h 	  = (int *)	malloc(numMems*numSpecs*numReacts * sizeof(int));
  varSend_h		  = (int *)	malloc(numMems*numSpecs*numReacts * sizeof(int));
  tgtMatrix_h	  = (int *)	malloc(numMems*numSpecs*numReacts * sizeof(int));
  tgtVector_h	  = (int *)	    malloc(numMems*numSpecs*numReacts * sizeof(int));
  Cvector_h		  = (float *)	malloc(numMems*numSpecs*numReacts * sizeof(float));

  if (buffBool_h)
  	  buffMAllSims_h   = (float *)  malloc(numSims * numMems*buffRows_h*buffCols_h * sizeof(float));

  for (unsigned i=0; i<numMems*numSpecs*numReacts; i++) {
  	  var_h[i] = 0;
  	  varSend_h[i] = 0;
  }

  /* For each membrane */

  if((fm0 = fopen ("input/M_0", "r")) == NULL){
	printf("Cannot open M_0\n");
	exit(0);
  }

  if((fcmatrix = fopen ("input/C_matrix", "r")) == NULL){
	printf("Cannot open the C values\n");
	exit(0);
  }


  if((fmfeed = fopen ("input/M_feed", "r")) == NULL){
	printf("Cannot open M_feed\n");
	exit(0);
  }


  if((fsize = fopen ("input/Size_matrix", "r")) == NULL){
	printf("Cannot open the size matrix\n");
	exit(0);
  }

  for (unsigned m=0; m<numMems; m++) {
  	  
  	  load_files2(m, fm0, fcmatrix, fmfeed, fsize);

  	  /* Convert iTgt_matrix2D_h and iTgt_vector2D_h into tgtMatrix_h and tgtVector_h */
  	  if (tgtKind_h == 1) {
  		  for (unsigned r=0; r<numReacts; r++)
  			  tgtVector_h[r*numMems+m] = iTgt_vector2D_h[m][r];
  		  free(iTgt_vector2D_h[m]);
  	  } else {
  		  for(unsigned r=0; r<numReacts; r++)
  			  for(unsigned c=0; c<numSpecs; c++)
  					  tgtMatrix_h[r*numMems*numSpecs + c*numMems + m] = iTgt_matrix2D_h[m][r*numSpecs + c];
  		  free(iTgt_matrix2D_h[m]);
  		  }

  	  /* Convert unLeft_side2D_h/unRight_side2D_h into leftSide_h/rightSide_h */
  	  /* Compute the VARIANTIONS MATRIX and print it into the LOG file */
  	  
  	  fprintf(log_fPtr, "\nComputed variations matrix var:\n");

  	  for(unsigned r=0; r<numReacts; r++) {
  		  for(unsigned c=0; c<numSpecs; c++) {
  			  leftSide_h [r*numMems*numSpecs + c*numMems + m] = unLeft_side2D_h[m][r*numSpecs+c];
  			  rightSide_h[r*numMems*numSpecs + c*numMems + m] = unRight_side2D_h[m][r*numSpecs+c];

  		      if (tgtKind_h == 1) {
  			  if(tgtVector_h[r*numMems+m] == -1) {
  				  var_h[r*numMems*numSpecs + c*numMems + m] -= leftSide_h [r*numMems*numSpecs + c*numMems + m];
  				  var_h[r*numMems*numSpecs + c*numMems + m] += rightSide_h[r*numMems*numSpecs + c*numMems + m];
  			  } else
  				  var_h[r*numMems*numSpecs + c*numMems + m] -= leftSide_h [r*numMems*numSpecs + c*numMems + m];
  		      } else {
  			  if (tgtMatrix_h[r*numMems*numSpecs + c*numMems + m] == m) {
  				  var_h[r*numMems*numSpecs + c*numMems + m] -= leftSide_h [r*numMems*numSpecs + c*numMems + m];
  				  var_h[r*numMems*numSpecs + c*numMems + m] += rightSide_h[r*numMems*numSpecs + c*numMems + m];
  				  varSend_h[r*numMems*numSpecs + c*numMems + m] = var_h[r*numMems*numSpecs + c*numMems + m];
  			  } else {
  				  var_h[r*numMems*numSpecs + c*numMems + m]	  -= leftSide_h [r*numMems*numSpecs + c*numMems + m];
  				  varSend_h[r*numMems*numSpecs + c*numMems + m] = rightSide_h[r*numMems*numSpecs + c*numMems + m];
  			  }
  		      }
  		      fprintf(log_fPtr,"%d\t",var_h[r*numMems*numSpecs + c*numMems + m]);
  		  }
  		  fprintf(log_fPtr,"\n");
  	  }
  	  free(unLeft_side2D_h[m]);
  	  free(unRight_side2D_h[m]);

  }
  free(unLeft_side2D_h); free(unRight_side2D_h); free(iTgt_matrix2D_h); free(iTgt_vector2D_h);
  fclose(fm0); fclose(fmfeed); fclose(fsize); fclose(fcmatrix);

  printf("Data 2 Acquired\n"); fflush(NULL);
}



float run(StopWatchInterface **timerOUT) {

  // Determine max threads per block
  int device = 0;
  cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  cudaCheckErrors("get device properties");

  // Check requested size is valid
  if (threadsPerBlock < 1 || threadsPerBlock > static_cast<unsigned int>(deviceProperties.maxThreadsPerBlock))
     printf("specified block size (%d) is invalid, must be between %d and %d for device %d.\n", threadsPerBlock, 1, deviceProperties.maxThreadsPerBlock, device);

  cudaFuncAttributes funcAttributes;
  dim3 block;
  dim3 grid;

  grid.x  = numSims;
  block.x = threadsPerBlock;

  // Get initRNG function properties and check the maximum block size
  cudaFuncGetAttributes(&funcAttributes, initRNG);
  cudaCheckErrors("function get attributes");
  if (block.x > (unsigned)funcAttributes.maxThreadsPerBlock)
  	  printf("Block X dimension is too large for initRNG kernel");
  	  
  // Get dynamicStep function properties and check the maximum block size
  cudaFuncGetAttributes(&funcAttributes, dynamicStep);
  cudaCheckErrors("function get attributes");
  if (block.x > (unsigned)funcAttributes.maxThreadsPerBlock)
  	  printf("Block X dimension is too large for dynamicStep kernel");

  // Check the maximum block size for device
  if (block.x > (unsigned)deviceProperties.maxThreadsDim[0])
  	  printf("Block X dimension for 'initRNG' kernel is too large for device");
  // Check the maximum grid size for device
  if (grid.x > (unsigned)deviceProperties.maxGridSize[0])
  	  printf("Grid X dimension is too large for device");

  // Allocate memory for RNG states 
  curandState *d_rngStates = 0;
  cudaMalloc((void **)&d_rngStates, grid.x * numMems * sizeof(curandState));
  cudaCheckErrors("memory allocation on device for RNG states");
  
  // Initialize RNG
  initRNG<<<grid, block>>>(d_rngStates, seed);

  StopWatchInterface *timerIN = NULL;
  sdkCreateTimer(&timerIN);

  sdkStartTimer(timerOUT);
  sdkStartTimer(&timerIN);

 
  dynamicStep<<<grid, block>>>(d_rngStates,
  	  buffMAllSims, buffRowCounterAllSims,
  	  step, MAllSims, t,
  	  leftSide, rightSide,
  	  var, varSend,
  	  tgtVector, tgtMatrix,
  	  Cvector,
  	  Mfeed, molSize, membSize,
  	  k_ruleAllSims, criticalAllSims, orderAllSims, R_probAllSims, R_probCAllSims,
  	  M_sendAllSims, t1AllSims, tauAllSims, FSAllSims,
  	  FSbkpAllSims, 
  	  tauSSAAllSims, 
  	  MbkpAllSims,
  	  flagSSAAllSims, flagStepAllSims,
  	  tauPreMinAllSims, tau1AllSims, tau2AllSims, a0cAllSims, a0AllSims,
  	  HORAllSims, auxVec1AllSims, auxVec2AllSims);

  CUDAcall(cudaDeviceSynchronize(),"dynamicStep");

  sdkStopTimer(timerOUT);
  sdkStopTimer(&timerIN);

  float elapsed = sdkGetAverageTimerValue(&timerIN);

  sdkDeleteTimer(&timerIN);

  return elapsed;
}

void GPUmallocData(){

	CUDAcall(cudaMalloc((void **)&seed, numSims * sizeof(unsigned long)),"seed");

	if (buffBool_h) {
		CUDAcall(cudaMalloc((void **)&buffMAllSims,  numSims * numMems*buffRows_h*buffCols_h * sizeof(float)),"buffM");
		CUDAcall(cudaMalloc((void **)&buffRowCounterAllSims,  numSims * numMems * sizeof(unsigned)),"buffRowCounter");
	}
	
	// System costants

	CUDAcall(cudaMalloc((void **)&Mfeed, numMems*numSpecs * sizeof(float)),"Mfeed");
	CUDAcall(cudaMalloc((void **)&molSize, numMems*numSpecs * sizeof(float)),"molSize");
	CUDAcall(cudaMalloc((void **)&membSize, numMems * sizeof(float)),"membSize");
	
	CUDAcall(cudaMalloc((void **)&leftSide, numMems*numReacts*numSpecs * sizeof(unsigned)),"leftSide");
	CUDAcall(cudaMalloc((void **)&rightSide, numMems*numReacts*numSpecs * sizeof(unsigned)),"rightSide");

	CUDAcall(cudaMalloc((void **)&var, numMems*numReacts*numSpecs * sizeof(int)),"var");
	CUDAcall(cudaMalloc((void **)&varSend, numMems*numReacts*numSpecs * sizeof(int)),"varSend");

	CUDAcall(cudaMalloc((void **)&tgtVector, numMems*numReacts * sizeof(int)),"tgtVector");
	CUDAcall(cudaMalloc((void **)&tgtMatrix, numMems*numReacts*numSpecs * sizeof(int)),"tgtMatrix");

	CUDAcall(cudaMalloc((void **)&Cvector, numMems*numReacts * sizeof(float)),"Cvector");

	// System variables

	CUDAcall(cudaMalloc((void**)&t,		numSims * sizeof(float)), "t");
	CUDAcall(cudaMalloc((void**)&step,	numSims * sizeof(unsigned)), "step");
	CUDAcall(cudaMalloc((void **)&MAllSims,	numSims * numMems * numSpecs * sizeof(float)),"MAllSims");

	CUDAcall(cudaMalloc((void **)&k_ruleAllSims,   numSims * numMems*numReacts * sizeof(unsigned)),"uliK_rule");
	CUDAcall(cudaMalloc((void **)&criticalAllSims, numSims * numMems*numReacts * sizeof(unsigned)),"uliCritical");
	CUDAcall(cudaMalloc((void **)&orderAllSims,    numSims * numMems*numReacts * sizeof(unsigned)),"order");
	CUDAcall(cudaMalloc((void **)&R_probAllSims,   numSims * numMems*numReacts * sizeof(float)),"R_prob");
	CUDAcall(cudaMalloc((void **)&R_probCAllSims,  numSims * numMems*numReacts * sizeof(float)),"R_probC");
	
	// ex SHARED memory variables
	CUDAcall(cudaMalloc((void**)&M_sendAllSims, numSims * numMems * numSpecs * sizeof(float)),"M_send");
	CUDAcall(cudaMalloc((void**)&t1AllSims, numSims * numMems * sizeof(float)),"t1");
	CUDAcall(cudaMalloc((void**)&tauAllSims, numSims * numMems * sizeof(float)),"tau");
	CUDAcall(cudaMalloc((void**)&FSAllSims, numSims * numMems * sizeof(float)),"FS");
	//

	CUDAcall(cudaMalloc((void**)&FSbkpAllSims, numSims * numMems * sizeof(float)),"FSbkp");
	CUDAcall(cudaMalloc((void**)&tauSSAAllSims, numSims * numMems * sizeof(float)),"tauSSA");
	
	CUDAcall(cudaMalloc((void**)&MbkpAllSims, numSims * numSpecs * numMems * sizeof(float)),"Mbkp");

	CUDAcall(cudaMalloc((void**)&flagSSAAllSims, numSims * numMems * sizeof(unsigned)),"flagSSA");
	CUDAcall(cudaMalloc((void**)&flagStepAllSims, numSims * numMems * sizeof(unsigned)),"flagStep");

	CUDAcall(cudaMalloc((void**)&tauPreMinAllSims, numSims * numMems * sizeof(float)),"tauPreMin");
	CUDAcall(cudaMalloc((void**)&tau1AllSims, numSims * numMems * sizeof(float)),"tau1");
	CUDAcall(cudaMalloc((void**)&tau2AllSims, numSims * numMems * sizeof(float)),"tau2");
	CUDAcall(cudaMalloc((void**)&a0cAllSims, numSims * numMems * sizeof(float)),"a0c");
	CUDAcall(cudaMalloc((void**)&a0AllSims, numSims * numMems * sizeof(float)),"a0");

	CUDAcall(cudaMalloc((void**)&HORAllSims, numSims * numSpecs * numMems * sizeof(unsigned)),"HOR");       
	CUDAcall(cudaMalloc((void**)&auxVec1AllSims, numSims * threadsPerBlock * numSpecs * sizeof(float)),"auxVec1");         
	CUDAcall(cudaMalloc((void**)&auxVec2AllSims, numSims * threadsPerBlock * numSpecs * sizeof(float)),"auxVec2");     

}
void GPUfreeData() {

	CUDAcall(cudaFree(seed),"seed");
	
	if (buffBool_h) {
		CUDAcall(cudaFree(buffMAllSims),"buffM");
		CUDAcall(cudaFree(buffRowCounterAllSims),"buffRowCounter");
	}

	// System costants

	CUDAcall(cudaFree(Mfeed),"Mfeed");
	CUDAcall(cudaFree(molSize),"molSize");
	CUDAcall(cudaFree(membSize),"membSize");
			
	CUDAcall(cudaFree(leftSide),"leftSide");
	CUDAcall(cudaFree(rightSide),"rightSide");
			
	CUDAcall(cudaFree(var),"var");
	CUDAcall(cudaFree(varSend),"varSend");
		
	CUDAcall(cudaFree(tgtVector),"tgtVector");
	CUDAcall(cudaFree(tgtMatrix),"tgtMatrix");

	CUDAcall(cudaFree(Cvector),"Cvector");

	// System variables

	CUDAcall(cudaFree(t),"t");
	CUDAcall(cudaFree(step),"step");
				
	CUDAcall(cudaFree(MAllSims),"MAllSims");	

	CUDAcall(cudaFree(k_ruleAllSims),"k_rule");
	CUDAcall(cudaFree(criticalAllSims),"critical");
	CUDAcall(cudaFree(orderAllSims),"order");
	CUDAcall(cudaFree(R_probAllSims),"R_prob");
	CUDAcall(cudaFree(R_probCAllSims),"R_probC");
	CUDAcall(cudaFree(M_sendAllSims),"M_send"); // ex SHARED memory
	CUDAcall(cudaFree(t1AllSims),"t1"); // ex SHARED memory
	CUDAcall(cudaFree(tauAllSims),"tau"); // ex SHARED memory
	CUDAcall(cudaFree(FSAllSims),"FS"); // ex SHARED memory

	CUDAcall(cudaFree(FSbkpAllSims),"FSbkp");
	CUDAcall(cudaFree(tauSSAAllSims),"tauSSA");
	
	CUDAcall(cudaFree(MbkpAllSims),"Mbkp");

	CUDAcall(cudaFree(flagSSAAllSims),"flagSSA");
	CUDAcall(cudaFree(flagStepAllSims),"flagStep");

	CUDAcall(cudaFree(tauPreMinAllSims),"tauPreMin");
	CUDAcall(cudaFree(tau1AllSims),"tau1");
	CUDAcall(cudaFree(tau2AllSims),"tau2");
	CUDAcall(cudaFree(a0cAllSims),"a0c");
	CUDAcall(cudaFree(a0AllSims),"a0");

	CUDAcall(cudaFree(HORAllSims),"HOR");       
	CUDAcall(cudaFree(auxVec1AllSims),"auxVec1");          
	CUDAcall(cudaFree(auxVec2AllSims),"auxVec2");    

}

void CPUfreeData() {

  free(seed_h);

  if (buffBool_h) {
  	  free(buffMAllSims_h);
  	  free(buffRowCounterAllSims_h);
  }

  // System costants
  			  
  free(membSize_h);
  		  
  free(leftSide_h);
  free(rightSide_h);
  		  
  free(var_h);
  free(varSend_h);
  	  
  free(tgtVector_h);
  free(tgtMatrix_h);
  		  
  free(Mfeed_h);
  free(molSize_h);
  		  
  free(Cvector_h);
  		  
  // System variables

  free(t_h);
  free(step_h);

  free(MAllSims_h);	  

  fclose(log_fPtr);
  for (unsigned s=0; s<numSims; s++) fclose(buffer_fPtr[s]);
}


void uploadData() {

	CUDAcall(cudaMemcpy(seed, seed_h, numSims * sizeof(unsigned long), cudaMemcpyHostToDevice),"seed");

	if (buffBool_h) {
		CUDAcall(cudaMemcpy(buffMAllSims, buffMAllSims_h, numSims * numMems*buffRows_h*buffCols_h* sizeof(float), cudaMemcpyHostToDevice),"buffM");
		CUDAcall(cudaMemcpy(buffRowCounterAllSims,	buffRowCounterAllSims_h, numSims * numMems * sizeof(unsigned), cudaMemcpyHostToDevice),"buffRowCounter");
	}

	CUDAcall(cudaMemcpyToSymbol(buffEvery, 	&buffEvery_h, sizeof(unsigned), 0, cudaMemcpyHostToDevice), "buffEvery");
	CUDAcall(cudaMemcpyToSymbol(buffRows, 	&buffRows_h,  sizeof(unsigned), 0, cudaMemcpyHostToDevice), "buffRows");
	CUDAcall(cudaMemcpyToSymbol(buffCols, 	&buffCols_h,  sizeof(unsigned), 0, cudaMemcpyHostToDevice), "buffCols");
	CUDAcall(cudaMemcpyToSymbol(buffBool, 	&buffBool_h,  sizeof(bool),     0, cudaMemcpyHostToDevice), "buffBool");
	CUDAcall(cudaMemcpyToSymbol(eps, 	&eps_h,       sizeof(float),    0, cudaMemcpyHostToDevice), "eps");
	
	
	// System costants

	CUDAcall(cudaMemcpyToSymbol(stepMax, 		&stepMax_h,     sizeof(float),    0, cudaMemcpyHostToDevice), "stepMax");
	CUDAcall(cudaMemcpyToSymbol(tgtKind, 	&tgtKind_h, 	 sizeof(unsigned), 0, cudaMemcpyHostToDevice), "tgtKind");

	CUDAcall(cudaMemcpy(leftSide, leftSide_h, numMems*numReacts*numSpecs *sizeof(unsigned), cudaMemcpyHostToDevice),"leftSide");
	CUDAcall(cudaMemcpy(rightSide, rightSide_h, numMems*numReacts*numSpecs * sizeof(unsigned), cudaMemcpyHostToDevice),"rightSide");

	CUDAcall(cudaMemcpy(var, var_h, numMems*numReacts*numSpecs * sizeof(int), cudaMemcpyHostToDevice),"var");
	CUDAcall(cudaMemcpy(varSend, varSend_h, numMems*numReacts*numSpecs * sizeof(int), cudaMemcpyHostToDevice),"varSend");

	CUDAcall(cudaMemcpy(tgtMatrix, tgtMatrix_h, numMems*numReacts*numSpecs * sizeof(int), cudaMemcpyHostToDevice),"tgtMatrix");
	CUDAcall(cudaMemcpy(tgtVector, tgtVector_h, numMems*numReacts * sizeof(int), cudaMemcpyHostToDevice),"tgtVector");

	CUDAcall(cudaMemcpy(Cvector, Cvector_h, numMems*numReacts * sizeof(float), cudaMemcpyHostToDevice),"Cvector");

	CUDAcall(cudaMemcpy(membSize, membSize_h, numMems * sizeof(float), cudaMemcpyHostToDevice),"membSize");
	CUDAcall(cudaMemcpy(Mfeed, Mfeed_h, numMems*numSpecs * sizeof(float), cudaMemcpyHostToDevice),"Mfeed");
	CUDAcall(cudaMemcpy(molSize, molSize_h, numMems*numSpecs * sizeof(float), cudaMemcpyHostToDevice),"molSize");

	// System variables

	CUDAcall(cudaMemcpy(MAllSims, MAllSims_h,		numSims * numMems * numSpecs * sizeof(float), cudaMemcpyHostToDevice),"MAllSims");
	
	CUDAcall(cudaMemset(orderAllSims, 0, numSims * numMems*numReacts * sizeof(unsigned)),"orderAllSims");
	CUDAcall(cudaMemset(tauSSAAllSims, 0, numSims * numMems * sizeof(unsigned)),"tauSSAAllSims");
	CUDAcall(cudaMemset(flagSSAAllSims, 0, numSims * numMems * sizeof(unsigned)),"flagSSAAllSims");

	CUDAcall(cudaMemset(tauAllSims, 0, numSims * numMems * sizeof(float)),"tauAllSims");
}
void downloadData() {

	if (buffBool_h) {
		CUDAcall(cudaMemcpy(buffMAllSims_h, buffMAllSims, numSims * numMems*buffRows_h*buffCols_h* sizeof(float), cudaMemcpyDeviceToHost),"buffM");
		CUDAcall(cudaMemcpy(buffRowCounterAllSims_h, buffRowCounterAllSims,	numSims * numMems * sizeof(unsigned), cudaMemcpyDeviceToHost),"buffRowCounter");
	}

	// System variables

	CUDAcall(cudaMemcpy(t_h, t,	numSims * sizeof(float), cudaMemcpyDeviceToHost), "t");
	CUDAcall(cudaMemcpy(step_h, step, numSims * sizeof(unsigned), cudaMemcpyDeviceToHost), "step");

	CUDAcall(cudaMemcpy(MAllSims_h, MAllSims,	numSims * numMems*numSpecs*sizeof(float), cudaMemcpyDeviceToHost),"MAllSims");
}
void printResult() {

	for (unsigned m=0; m<numMems; m++)
		for (unsigned s=0; s<numSims; s++) {
			//write_matrix_bin(m, s);
			write_matrix_txt(m, s);
		}

}

int main(int argc, char **argv) {

  cudaDeviceProp deviceProperties;

  // Get number of available devices
  int deviceCount = 0;
  cudaError_t cudaResult = cudaGetDeviceCount(&deviceCount);
  cudaCheckErrors("get device count");

  // Get device properties
  int device = 0;
  cudaGetDeviceProperties(&deviceProperties, device);
  cudaCheckErrors("get device properties");

  // Check precision
  if (deviceProperties.major < 1 || (deviceProperties.major == 1 && deviceProperties.minor < 3))
  	  printf("device does not have float precision support");

  // Attach to GPU
  cudaSetDevice(device);
  cudaCheckErrors("set device");

  bool bBufferOut = 0;

  char folder[500];
  std::sprintf(folder, "mkdir output%u_%u_%u", numMems, threadsPerBlock, numSims);
  system(folder);

  char prompt[500];
  std::sprintf(prompt, "output%u_%u_%u/prompt_s%u_m%u_buff%d.txt", numMems, threadsPerBlock, numSims, numSims, numMems, bBufferOut);
  /*
  if (std::freopen(prompt, "w", stdout) == NULL) {
  	      printf ("cannot open %s", prompt);
  	      exit(0);
  }
  */
  StopWatchInterface *timer_ALLOC_UPLOAD = NULL, *timer_FREE_DOWNLOAD = NULL;
  StopWatchInterface *timer_CALC = NULL;
  sdkCreateTimer(&timer_ALLOC_UPLOAD);
  sdkCreateTimer(&timer_FREE_DOWNLOAD);
  sdkCreateTimer(&timer_CALC);

  float *elapsedTime_CALC = (float*) malloc(numTrials * sizeof (float));
  unsigned *numSteps = (unsigned*) malloc(numTrials * numSims * sizeof (unsigned));
  unsigned *numStepsAllSims = (unsigned*) malloc(numTrials * sizeof (unsigned));
  float timePerStepMedium_CALC = 0;

  size_t varSize=0, parSize=0;

  for (unsigned trial=0; trial<numTrials; trial++) {
	CPUinitializeData(bBufferOut);
  	
  	sdkStartTimer(&timer_ALLOC_UPLOAD);
  	GPUmallocData();
  	uploadData();
  	sdkStopTimer(&timer_ALLOC_UPLOAD);

	//Main Function
  	elapsedTime_CALC[trial] = run(&timer_CALC);

  	sdkStartTimer(&timer_FREE_DOWNLOAD);
  	downloadData();
  	GPUfreeData();
  	sdkStopTimer(&timer_FREE_DOWNLOAD);

  	numStepsAllSims[trial] = 0;
  	for (unsigned s=0; s<numSims; s++) {
  	        numSteps[trial*numSims + s] = step_h[s];
  	        numStepsAllSims[trial] += numSteps[trial*numSims + s];
  	}

  	timePerStepMedium_CALC += elapsedTime_CALC[trial]/numStepsAllSims[trial];

  	printResult();

  	if (trial==0) {
  	        varSize = numSims * (3 * numMems*numReacts * sizeof(unsigned) +
  	        			2 * numMems*numReacts * sizeof(float) +
  	        			4 * numMems * numSpecs * sizeof(float) +
  	        			10 * numMems * sizeof(float) +
  	        			2 * numSims * numMems * sizeof(unsigned) +
  	        			1 * numSims * numSpecs * numMems * sizeof(unsigned));

  	        parSize = numSims * numMems * sizeof(curandState) + // rngStates
  	        			numMems * sizeof(unsigned) + // numReacts
  	        			numSims * sizeof(unsigned) + // step
  	        			numSims * numMems * numSpecs * sizeof(float) + // M
  	        			numSims * sizeof(float) + // t
  	        			2 * numMems*numReacts*numSpecs * sizeof(unsigned) + // leftSide,rightSide
  	        			2 * numMems*numReacts*numSpecs * sizeof(int) + // var, varSend
  	        			numMems*numReacts * sizeof(int) + // tgtVector
  	        			numMems*numReacts*numSpecs * sizeof(int) + // tgtMatrix
  	        			numMems*numReacts * sizeof(float) + // Cvector
  	        			2 * numMems * numSpecs * sizeof(float) + // Mfeed, molSize
  	        			numMems * sizeof(float); // membSize

  	        if (bBufferOut)
  	        	parSize += numSims * numMems*buffRows_h*(numSpecs+1) * sizeof(float) + // buffM
  	        			numSims * numMems * sizeof(unsigned); // buffRowCounter
  	        /*
  	        k_rule, critical, order: numSims * numMems*numReacts * sizeof(unsigned)
  	        R_prob, R_probC: numSims * numMems*numReacts * sizeof(float)
  	        M_send, Mbkp, auxVec1, auxVec2: numSims * numMems * numSpecs * sizeof(float)
  	        t1, tau, FS, FSbkp, tauSSA, tauPreMin, tau1, tau2, a0c, a0: numSims * numMems * sizeof(float)
  	        flagSSA, flagStep: numSims*numMems*sizeof(unsigned)
  	        HOR: numSims*numSpecs*numMems*sizeof(unsigned)  
  	        rngStates: numSims * numMems * sizeof(curandState)
  	        buffM: numSims * numMems*buffRows_h*(numSpecs+1) * sizeof(float)
  	        buffRowCounter: numSims * numMems * sizeof(unsigned)
  	        numReacts: numMems * sizeof(unsigned)
  	        step: numSims * sizeof(unsigned)
  	        M: numSims * numMems * numSpecs * sizeof(float)
  	        t: numSims * sizeof(float)
  	        leftSide: numMems*numReacts*numSpecs * sizeof(unsigned)
  	        rightSide: numMems*numReacts*numSpecs * sizeof(unsigned)
  	        var: numMems*numReacts*numSpecs * sizeof(int)
  	        varSend: numMems*numReacts*numSpecs * sizeof(int)
  	        tgtVector: numMems*numReacts * sizeof(int)
  	        tgtMatrix: numMems*numReacts*numSpecs * sizeof(int)
  	        Cvector: numMems*numReacts * sizeof(float)
  	        Mfeed: numMems*numSpecs * sizeof(float)
  	        molSize: numMems*numSpecs * sizeof(float)
  	        membSize: numMems * sizeof(float)
  	        */
  	}

  	CPUfreeData();

  }

  float elapsedTimeMedium_ALLOC_UPLOAD = sdkGetAverageTimerValue(&timer_ALLOC_UPLOAD);
  float elapsedTimeMedium_FREE_DOWNLOAD = sdkGetAverageTimerValue(&timer_FREE_DOWNLOAD);
  float elapsedTimeMedium_CALC = sdkGetAverageTimerValue(&timer_CALC);
  timePerStepMedium_CALC /= numTrials;

  printf("\nnumMems = %u \nnumSims = %u \nnumTrials = %u", numMems, numSims, numTrials);

  printf("\nTotal GLOBAL memory occupied by simulations and system parameters: %u (vars) + %u (pars) = %u bytes \n", varSize, parSize, varSize + parSize);

  printf("\nelapsedTimeMedium_ALLOC_UPLOAD = %.3f(ms) \nelapsedTimeMedium_FREE_DOWNLOAD = %.3f(ms) \nelapsedTimeMedium_CALC = %.3f(ms)  \ntimePerStepMedium_CALC = %.3f(ms)", elapsedTimeMedium_ALLOC_UPLOAD, elapsedTimeMedium_FREE_DOWNLOAD, elapsedTimeMedium_CALC, timePerStepMedium_CALC);

  

  for (unsigned trial=0; trial<numTrials; trial++) {
  	printf("\n\nelapsedTime_CALC[%u] = %.3f(ms), total steps %u\n", trial, elapsedTime_CALC[trial], numStepsAllSims[trial]);
  	for (unsigned s=0; s<numSims; s++) printf("numSteps[%4u] = %u\n", s, numSteps[trial*numSims+s]);
  }

  sdkDeleteTimer(&timer_ALLOC_UPLOAD);
  sdkDeleteTimer(&timer_FREE_DOWNLOAD);
  sdkDeleteTimer(&timer_CALC);

  free(elapsedTime_CALC);
  free(numSteps);
  free(numStepsAllSims);

  fclose(stdout);

  cudaDeviceReset();
  cudaCheckErrors("device reset");

  exit(EXIT_SUCCESS);
}
