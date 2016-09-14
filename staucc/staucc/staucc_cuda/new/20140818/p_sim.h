#define INV_RNDMAX 1.0 / RAND_MAX
#define N_A  6.0221415e23L;
#define LOGFILE

typedef  unsigned long int ULI;
typedef  long int LI;
typedef  long double LD;
struct entry {
	int iRows;
	ULI **uliMatr;
	struct entry *next;
};

struct globale {
  FILE **m_fPtr;
  FILE **log_fPtr; ///potrebbe essere uno solo
  FILE **fs_fPtr;
  //FILE **rprob_fPtr;

  FILE* flmatrix;
  FILE* frmatrix;  
  FILE* fcmatrix;   
  FILE* ftgtmatrix;   
  FILE* fm0;
  FILE* fmfeed;
  FILE* findexes;  
  FILE* fsize;
    
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
  
  LD **ldM_matrix;

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
  LD eps;
  
};
