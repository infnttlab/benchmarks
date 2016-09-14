#define INV_RNDMAX 1.0 / RAND_MAX
#define N_A  6.0221415e23L;
#define LOGFILE
#define APPSIZE 20 //usato per lunghezza numero file input
typedef  unsigned long int ULI;
typedef  long int LI;
typedef  long double LD;

/*
struct entry {
	int iRows;
	ULI **uliMatr;
	struct entry *next;
};
*/
struct globale {

  ////restano costanti
  unsigned *unLeft_side;	
  unsigned *unRight_side;	
  int *iVar;		     
  int *iVarSend;	     

  int * side_index;

  ULI *uliRows;			//di solito diversa tra le membrane
  ULI *uliColumns;		//uguale per ogni membrana

  int **iTgt_vector;
  int *iTgt_matrix;

  LD **ldC_vector;
  LD **ldM_feed;
  LD **ldMolSize;
  LD *ldMembSize;


  ////ausiliari
  FILE **m_fPtr;		//file di output
  FILE **log_fPtr;

  LD **ldM_matrix;		//buffer per la scrittura dell'evoluzione del sistema
  ULI *uliRow_counter;		//controlla il buffer ldM_matrix


  ////realmente private e cambiano
  ULI **uliK_rule;
  ULI **uliCritical;

  LD *tau;
  LD *ldM_send;			//resettata ad ogni giro del main loop due volte
  LD *ldTau_ssa;
  LD *ldFS;
  LD *ldFSbkp;
        
  LD **M;
  LD **Mbkp;

  LD **ldR_prob;
  LD **ldR_prob_c;

  int *iFlag_SSA;
  int *iFlagStep;
  
  int iTgtcheck;
  ULI uliIndexCol;
  ULI uliMAXSEED;
  int every;
  int buffer_lines;
  int iotime;
  
  //int* checkFS; diventato inutile perche' lo controllo come output di chechnegfs  
  //ULI **uliR_app; sembra inutile , lo assegno ma non lo uso
  //int **indexes; inutile, letto e scritto nella lettura file e non piu' usato  
  //ULI **uliOrder; usata solo in get_HOR
  //LD **ldHOR; usata solo in dpp_step1
};
