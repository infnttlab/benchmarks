//-----------------------------------------------------------------------------
// Questo programma è stato creato da Rosa Brancaccio
//-----------------------------------------------------------------------------
// Dipartimento di Fisica
// Sezione Fisica Medica
// V.le Berti Pichat 6/2
// 40127 Bologna - ITALY
// 
// Tel.:   +39-051-2095132
// Tel.:   +39-051-2095134
// Fax.:   +39-051-2095047
// Email:  rossella_brancaccio@yahoo.it
// Web:    http://www.xraytomography.com
//-----------------------------------------------------------------------------
// REC_FUNCTION.h
// mio file di intestazione con funzioni di REC FUNCTIONS
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// variabili global per ricostruzione
//-----------------------------------------------------------------------------
extern float   *filter;
extern float   *row_re,*row_im;
extern float    *slice,*slice_rank;
extern float	*detector_weighing;
extern int	 h_start,h_end;
extern int	*x_start,*x_end;			// vettori con i punti x di partenza e arrivo in funzione di y
extern float	*row_modified;					// riga modificata per halfscan si allocherà solo se servirà
//-----------------------------------------------------------------------------
// per allocare row e filter è necessario calcolare la lunghezza degli array nel dominio di fourier
// è necessario avere almeno il doppio meno uno della dimensione dei punti da elaborare
// inoltre se si usa un array con lunghezza potenza di 2 la FFT lavora meglio per cui:
//-----------------------------------------------------------------------------
extern int CalcTheMinimumNumberPowerOfTwoAfter(int dim);
//-----------------------------------------------------------------------------
// pesatura per geometria di scanning: nel fan beam la pesatura dei punti del 
// rivelatore va fatta in funzione della geometria e si distingue in 
// rivelatori con raggi equispatiati o equiangolari (vedere il Kak)
// ora calcola il peso anche per feldkamp!!! la differenza sta nel fatto che dw diventa 2D
//-----------------------------------------------------------------------------
extern int CalculateWeighing_OnTheBaseOfCTGeometry(int detector_geometry,int FFT_start,int nrays,float sod,float pxsize,float pxcenter);
//-----------------------------------------------------------------------------
// calcolo del filtro da applicare nello spazio di Fourier a partire da RAMP e con BUTTER
// bisogna passargli il tipo di W2 che vogliamo applicare
//-----------------------------------------------------------------------------
extern int CalculateFourierFilter(int FFT_lenght, int FFT_start, double W2);
//---------------------------------------------------------
// calcolo del peso per equispaziati U --- W3 ---
// questa funzione alloca l'array Uquadro_teta_zero e lo calcola
//---------------------------------------------------------
extern int	CalcolateWeighing_EquispacedRays(float sod,float rxsize);
//-----------------------------------------------------------------------------
// funzione per retroproiezione filtrata CON constrain circle SOLO MATEMATICA
// N.B. PRIMA DI ACCEDERE A QUESTA FUNZIONE BISOGNA ALLOCARE SLICE ROW RE ROW IM E FILTER
// INOLTRE è NECESSARIO AVER CARICATO IN MEMORIA IMAGE DA ELABORARE E AVER CALCOLATO FILTER E PESI
// ------ FUNZIONE UNICA INDIPENDETEMENTE DA CONSTRAIN si usano sempre xstart e xend
// // **** 2013 07 19 ANCHE LOCAL CT **** SENZA HALFSCAN
//-----------------------------------------------------------------------------
extern int BackProjectFilteredDataOnSlice(int FFT_lenght,int FFT_start,float xs,float ys,float xd0,float yd0,float xd1,float yd1,int numBlocks,int threadsPerBlock,int dimBlock);
//-----------------------------------------------------------------------------
// calcolo del peso --- W3 --- vale 1/Uquadro per equispaziati
//-----------------------------------------------------------------------------
extern int CalculateWeight3(int detector_geometry,float l,double sod_quadro,double xs_rot_quadro,double ys_rot_quadro,float xs_rot,float ys_rot,float *W3_U);
//---------------------------------------------------------
// funzione per il calcolo della geometria necessaria alla ricostruzione FBP
//---------------------------------------------------------
extern void CalcAllGeometryData_FBP(float *xs,float *ys,float *xd0,float *yd0,float *xd1,float *yd1,int *y_start,int *y_end);
//---------------------------------------------------------
// unica funzione che in base a DataSet calcola il raggio r e fan angle
// calcoliamo fanangle e r con - shift 
//---------------------------------------------------------
extern void	Calculate_r_fi(int geometry,float *r,float *fi);
//-----------------------------------------------------------------------------
//  valore assoluto di un float
//-----------------------------------------------------------------------------
extern float Absolute(float val);
//-----------------------------------------------------------------------------
//	si compone pathname dai dati inseriti in data set
// se n=-1 allora dark se n=-2 allora bak 
//-----------------------------------------------------------------------------
extern void	ComposeFileNameFromDataSet(char *PathName);
//-----------------------------------------------------------------------------
// copia string_tocopy in string_target e poi la termina
//-----------------------------------------------------------------------------
extern void CopyStringToStringAndTerminateIt(char *string_target,char *string_tocopy);
//-----------------------------------------------------------------------------
// si calcolano max e min dell'immagine corrente in base al data set
//-----------------------------------------------------------------------------
extern void MaxMin(int dataset);
//-----------------------------------------------------------------------------
// si salva l'immagine in input come SDT nel path specificato e con le dimensioni specificate 
// ci sono due funzioni: questa è per i float
//-----------------------------------------------------------------------------
extern int Save_SDT_GenericFile_FLT(float *image,char *path,int width,int height);
//-----------------------------------------------------------------------------
// si salva l'SPR nel path specificato e con le dimensioni specificate 
// e con i paramteri specificati
// **** size_policy: SI SALVANO ANCHE PSIZE O RSIZE SE COSì SI VUOLE *************
//-----------------------------------------------------------------------------
extern int Save_SPR_GenericFile(char *path,int width,int height,int IMG_TYPE,int DATA_TYPE,int size_policy);  // IMG_TYPE "LINE 1 IMG_2D 2 IMG_3D 3" DATA_TYPE "1=U16, 5=U8, 3=FLT, 6=RGB"  
//-----------------------------------------------------------------------------
// FFT by numerical recipes in c
// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
// be an integer power of 2 (this is not checked for!).
//-----------------------------------------------------------------------------
extern void four1(float data[], unsigned long nn, int isign);
//-----------------------------------------------------------------------------
// fast fourier transform
//-----------------------------------------------------------------------------
extern void	FFT(float *row_re,float *row_im,int FFT_lenght);
//-----------------------------------------------------------------------------
// inverse fft: inverse fast fourier transform
//-----------------------------------------------------------------------------
extern void	InvFFT(float *row_re,float *row_im,int FFT_lenght);
//-----------------------------------------------------------------------------
// Stampa del dataset scelto nello standard output
//-----------------------------------------------------------------------------
extern void PrintData_DataSet(void);
//-----------------------------------------------------------------------------
//	si apre il file sct per leggere alcuni parametri
//  in ogni caso se il file NON esiste o se è stato premuto CANCEL si ritorna -1
//-----------------------------------------------------------------------------
extern int Open_SCT_File(void);
//-----------------------------------------------------------------------------
// si apre un file SDT
//-----------------------------------------------------------------------------
extern int Open_SDT_File(int n);
//-----------------------------------------------------------------------------
// Si apre un file spr e si restituiscono tutte le caratteristiche dell'immagine 
//-----------------------------------------------------------------------------
extern int OpenSPRFile(char *PathName);
//-----------------------------------------------------------------------------
// si salva l'immagine in memoria SDT
// se in path c'è tomostep allora si salva automaticamente la sequenza e SOVRASCRIVE
// altrimenti si usa il path specificato e si chiede prima di sovrascrivere
// ritorna -1 se anche con 99 tentativi non si è riusciti a salvare
// ritorna -2 se il tipo di file non è SDT salvabile
//-----------------------------------------------------------------------------
extern int Save_SDT_File(char *path);
//-----------------------------------------------------------------------------
// si salva il file SPR
// se in path c'è tomostep allora si salva automaticamente la sequenza e SOVRASCRIVE
// altrimenti si usa il path specificato e si chiede prima di sovrascrivere
//-----------------------------------------------------------------------------
extern int Save_SPR_File(char *path);
//-----------------------------------------------------------------------------
// free delle allocazioni di un dataset generico passato come valore int
// si dealloca la memoria solo se DataSet.new_dataset è FALSE!!!!
// per essere sicuri di questo si richiama prima il free e poi le inizializzazioni
//-----------------------------------------------------------------------------
extern void Free_AllocatedMemory_DataSet(int dataset);
//-----------------------------------------------------------------------------
// Inizializzazioni di un dataset generico passato come valore int
// si inizializza tutto TRANNE i puntatori da allocare
//-----------------------------------------------------------------------------
extern void Initialize_Data_DataSet(int dataset);
//-----------------------------------------------------------------------------
// Inizializzazioni di un dataset generico passato come valore int
// si inizializzano SOLO i puntatori da allocare
//-----------------------------------------------------------------------------
extern void Initialize_Pointers_DataSet(int dataset);
//-----------------------------------------------------------------------------
// si fa partire la ricostruzione con PREVIEW con formule EQUISPAZIATI
// qui si tracciano le rette ed è molto più veloce con subFUNZIONI
// INTERPOLAZIONE SUL DETECTOR ------- NUOVO ------------- SENZA PREVIEW DA SCRIVERE!!!!
// ------- SENZA PREVIEW ----------- CON CONSTRAIN CIRCLE -------------
// // **** 2013 07 19 ANCHE LOCAL CT **** ANCHE HALFSCAN
//-----------------------------------------------------------------------------
extern int	Reconstruction_FAN_BEAM_detector_interpolation_NOPreview(int start,int end,int numBlocks,int threadsPerBlock,int dimBlock);
//-----------------------------------------------------------------------------
















