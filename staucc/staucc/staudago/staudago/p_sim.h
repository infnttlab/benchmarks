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
  FILE **log_fPtr;
  //FILE **fs_fPtr;
  FILE **rprob_fPtr;

  ULI *uliRows;			//di solito diversa tra le membrane
  ULI *uliColumns;		//uguale per ogni membrana
  ULI *uliRow_counter;		//controlla il buffer ldM_matrix

  ULI **uliK_rule;
  ULI **uliR_app;
  ULI **uliOrder;
  ULI **uliCritical;

  LD *tau;
  LD *ldM_send;			//resettata ad ogni giro del main loop due volte
  //LD *ldM_recv;
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
  
  LD **ldM_matrix;		//buffer per la scrittura dell'evoluzione del sistema

  int *iFlag_SSA;
  int *iFlagStep;
  int* checkFS;
  
  int **indexes;
  int **iTgt_vector;
  
  int ***iTgt_matrix;
  int ***iVar;			//assegnata all'inizio e resta costante
  int ***iVarSend;		//idem

  unsigned ***unLeft_side;	//letto da file, resta costante
  unsigned ***unRight_side;	//letto da file, resta costante

  struct entry **LSlist;
  struct entry **RSlist;

  int iTgtcheck;
  ULI uliIndexCol;
  ULI uliMAXSEED;
  int every;
  int buffer_lines;
  int iotime;
};
