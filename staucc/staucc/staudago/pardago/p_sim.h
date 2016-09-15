#define INV_RNDMAX 1.0 / RAND_MAX
#define N_A  6.0221415e23L;

typedef  unsigned long int ULI;
typedef  long int LI;
typedef  long double LD;
struct entry {
	int iRows;
	ULI **uliMatr;
	struct entry *next;
};


/*
extern FILE **m_fPtr, **log_fPtr, **fs_fPtr, **rprob_fPtr;
extern ULI *uliRows, *uliColumns, uliIndexCol;
extern ULI **uliK_rule, **uliR_app;
extern ULI **uliOrder, **uliCritical, *uliRow_counter;
extern ULI uliMAXSEED;
extern LD **ldC_vector, *tau, **M, **Mbkp, **ldM_feed, **ldR_prob;
extern LD **ldR_prob_c, *ldTau_ssa;
extern LD **ldHOR, every, buffer_lines;

extern LD ***ldM_matrix;
extern LD *ldM_send, *ldM_recv;

extern LD **ldMolSize, *ldMembSize, *ldFS, *ldFSbkp;
extern int  **iTgt_vector, ***iTgt_matrix, iTgtcheck;
extern int ***iVar, *iFlag_SSA, **indexes, *iFlagStep, ***iVarSend;
extern int *checkFS;
extern unsigned ***unLeft_side, ***unRight_side;
extern struct entry **LSlist, **RSlist;
*/

struct globale {
  FILE **m_fPtr;
  FILE **log_fPtr;
  FILE **fs_fPtr;
  FILE **rprob_fPtr;

  ULI *uliRows;
  ULI *uliColumns;
  ULI *uliRow_counter;

  ULI **uliK_rule;
  ULI **uliR_app;
  ULI **uliOrder;
  ULI **uliCritical;

  LD *tau;
  LD *ldM_send;
  LD *ldM_recv;
  LD *ldTau_ssa;
  LD *ldMembSize;
  LD *ldFS;
  LD *ldFSbkp;
  
  LD **ldC_vector;
  LD **M;
  LD **Mbkp;
  LD **ldM_feed;
  LD **ldR_prob;
  LD **ldR_prob_c;
  LD **ldHOR;
  LD **ldMolSize;
  
  LD ***ldM_matrix;

  int *iFlag_SSA;
  int *iFlagStep;
  int* checkFS;
  
  int **indexes;
  int **iTgt_vector;
  
  int ***iTgt_matrix;
  int ***iVar;
  int ***iVarSend;

  unsigned ***unLeft_side;
  unsigned ***unRight_side;

  struct entry **LSlist;
  struct entry **RSlist;

  int iTgtcheck;
  ULI uliIndexCol;
  ULI uliMAXSEED;
  int every;
  int buffer_lines;
};


//void (*dpp_update)(int, gsl_rng *, int, struct globale*);
