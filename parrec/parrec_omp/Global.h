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
// GLOBAL.H
// mio file di intestazione con variabili globali, define, funzioni
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// miei include 
//-----------------------------------------------------------------------------
//#include "Global.h"
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// DEFTYPES.H
// mio file di intestazione con define e type
//-----------------------------------------------------------------------------

	//-----------------------------------------------------------------------------
	//  define
	//-----------------------------------------------------------------------------
	#define DUEPIGRECO		6.283185307179586476925286766558
	#define PIGRECO			3.141592653589793238462643383279
	#define PIGRECOMEZZI	1.5707963267948966192313216916395
	#define PIGRECOQUARTI	0.78539816339744830961566084581975

	//-----------------------------------------------------------------------------
	// define da cvi tolto
	//-----------------------------------------------------------------------------
	#define MAX_PATHNAME_LEN 260
	//typedef int                 ssize_t;
	//-----------------------------------------------------------------------------
	// define per la FFT e invFFT
	//-----------------------------------------------------------------------------
	#include <math.h>
	#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
	//-----------------------------------------------------------------------------
	// minimi e massimi dei vari tipi di bit o char
	//-----------------------------------------------------------------------------
	#define	MAX_STRING 9999  // dimensioni stringhe lunghe
	#define MEDIUM_STRING 999 // dimensioni stringhe medie
	#define SHORT2_STRING 99
	#define	SHORT_STRING 9   // dimensioni stringhe corte
	#define MAX_8BIT 255	 
//	#define MAX_FLOAT 99.999  // rossella ?????? come lo scelgo???? ora si trova in global.c
	#define MAX_16BIT 65536
	#define MAX_HIST 65536
	#define ZERO 0
	#define ZERO_F 0. //ora si trova in global.c come MIN_FLOAT
//	#define ZERO_F 0. ora si trova in global.c come MIN_FLOAT
	#define MAX_NUM_FILES 9999 // massimo diecimila file in una sequenza
	#define RGB_BYTE 3
	//-----------------------------------------------------------------------------
	#define NO_VALUE -99		// valori per la definizione di max e min in seq
	#define NO_VALUE_F -99.		// valori per la definizione di maxf e minf in seq
	//-----------------------------------------------------------------------------
	#define TRUE 1
	#define FALSE 0
	//-----------------------------------------------------------------------------
	#define YES 1
	#define NO 0
	//-----------------------------------------------------------------------------
	// data type delle immagini, 1 e 3 sono standard, 5 e 6 li ho definiti io
	// LLNL: 0-byte, 1-16-bit word, 2-32-bit integer, 3-32-bit floating point. 
	//-----------------------------------------------------------------------------
	#define U16 1   // corrispondente a LLNL 1
	#define U8 5   	// non più gestito come U16 
	#define FLT 3 	// corrispondente a LLNL 3
	#define RGB_ 6  	// attualmente non in uso ma serve per immagini a colori
#define DATATYPE_STRING "1=U16, 5=U8, 3=FLOAT, 6=RGB"   // così se cambio qui qualcosa cambia automaticamente ovunque
	//-----------------------------------------------------------------------------
	#define SDT 1
	#define TIF 2
	#define JPG 3
	#define BMP 4
#define FILETYPE_STRING "1=SDT 2=TIF 3=JPG 4=BMP" // idem sopra
	//-----------------------------------------------------------------------------
	#define LINE 1
	#define IMG_2D 2
	#define IMG_3D 3
	//-----------------------------------------------------------------------------
	// gestione file open save ed error in particolare di GetRadixName_FromFileName
	// e di GetExtension_Name_FromFileName
	//-----------------------------------------------------------------------------
	#define FNAME_NOT_VALID -2
	#define UNDERSCORE_NOT_FOUND -1
	#define NUM_NOT_FOUND -9
	//-----------------------------------------------------------------------------
	// define per scoprire se il file dei link c'è non c'è è completo o no
	//-----------------------------------------------------------------------------
	#define	NOT_KNOWN							0
	#define	FILE_EXISTS							1
	#define MAX_LINK							10 // massimo numero di link nel menu
	//-----------------------------------------------------------------------------
	// tutte le inzializzazioni di DataSet
	//-----------------------------------------------------------------------------
	#define NO_WIDTH 0
	#define NO_HEIGHT 0
	//-----------------------------------------------------------------------------
	// ora ci sono due DataSet quindi bisogna dire a quale data set ci si riferisce
	//-----------------------------------------------------------------------------
	#define DATASET1 0
	#define DATASET2 1
	//-----------------------------------------------------------------------------
	// define dell'istogramma grafico
	//-----------------------------------------------------------------------------
	#define HIST_WIDTH 500
	#define HIST_HEIGHT 200
	#define HIST_LEFT 30
	#define HIST_TOP 30
	//-----------------------------------------------------------------------------
	// define dell'immagine 1 grafica
	//-----------------------------------------------------------------------------
	#define IMAGE_WIDTH 100
	#define IMAGE_HEIGHT 100
	//-----------------------------------------------------------------------------
	// define dei flag di stretch
	//-----------------------------------------------------------------------------
	#define AUTO TRUE		   // stretch automatico a max e min trovati nell'immagine
	#define PERSONALIZED FALSE // usa i min e max inseriti nei due box dell'istogramma
	#define STRETCH TRUE	   // istogramma calcolato su immagine stretch
	#define ORIGINAL FALSE     // istogramma calcolato su immagine originale
	//-----------------------------------------------------------------------------
	// define di int TomographicType  // indica lo step tomografico, cioè il tipo 
	// 						di immagini può essere projection atenrad sinos recobj
	//-----------------------------------------------------------------------------
	#define PROJECTION 1   
	#define ATENRAD 2
	#define SINOS 3
	#define RECOBJ 4
	#define GENERIC 5
	#define NOT_CLASSIFIED 6	 // mi serve per immagini in sequenza non riconosciute come step tomografico
	//-----------------------------------------------------------------------------
	// è un tipo di DataSet.numberofimages
	#define SINGLE_IMAGE 1	     // mi serve per immagini singole che non fanno parte di una sequenza 
	#define NO_IMAGE 0
	//-----------------------------------------------------------------------------
	// define per distinguere bits bits_preview bits_preview2 si usano in SetColorInBits
	//-----------------------------------------------------------------------------
	#define BITS          1
	#define BITS_PREVIEW  2
	#define BITS_PREVIEW2 3
	#define IMAGE1			4
	#define IMAGE2			5
	//-----------------------------------------------------------------------------
	#define REVERT_TO_ATENRAD		6
	#define RINGO_Filter			7
	#define RINGO_FFT				8
	#define METAL_ARTIFACT_Filter	9
	#define METAL_ARTIFACT_FFT		10
	#define OUTLIER					11
	#define OUTLIER_NEW				12
	#define FFT_FILTER				13
	#define RINGO_NEW				14
	#define COMPLETE_SINOS			15
	//-----------------------------------------------------------------------------
	#define SPIKES					14	// contrassegno che quel punto è spikes
	#define NO_SPIKES				15	// contrassegno che quel punto NON è spikes
	#define SPIKES_RIPESCATO		16	// spikes in più
	#define SPIKES_ORIZZONTALE		17	// spikes in più in orizzontale
	#define SPIKES_DEVSMEDIA_SINGOLO	18	// ulteriori spikes non trovati cercati ricalcolando la media
	#define SPIKES_DEVSMEDIA_DOPPIO		19	// ulteriori spikes non trovati cercati ricalcolando la media
	#define SPIKES_DEVSMEDIA			20	// ulteriori spikes non trovati cercati ricalcolando la media
	#define FUORI_SINO				0.   // punti fuori con sinogramma fuori dall'immagine calccenter
	//-----------------------------------------------------------------------------
	#define SINISTRA	0   // mi servono in complete sinos
	#define DESTRA		1
	//-----------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------
	// DataSet.current = -9 se non è una sequenza (può anche essere singola immagine ma con u numero e una radice)
	#define NO_SEQ -9
	#define DARK -1
	#define IZERO -2
	//-----------------------------------------------------------------------------
	// define di size_policy per decidedere se salvare size in SPR
	//-----------------------------------------------------------------------------
	#define SAVE_PSIZE 0
	#define SAVE_RSIZE 1
	#define SAVE_ZERO  2
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	#define PROFILE_X 1
	#define PROFILE_Y 0
	//-----------------------------------------------------------------------------
	// UTILITY STEP
	//-----------------------------------------------------------------------------
	#define NO_STEP 0
	#define CROP 1
	#define BINNING 2
	#define TRESHOLD_GRAYVALUE 3
	#define TRESHOLD_GEOMETRIC 4
	#define STRETCH_SEQUENCE 5
	#define EDGE_SEARCH 6
	#define COLLATEH 7
	#define COLLATEV 8
	#define ROTATION 9
	#define MIRRORH 10
	#define MIRRORV 11
	#define ALGEBRAIC_OPERATION 12
	#define SUBTRACTION 13
	#define	SIGNAL_BACKGROUND 14
	#define FILTER_AVERAGE	15
	#define FILTER_MEDIAN	16
	#define FILTER_GAUSS	17
	#define DECIMATE_SEQUENCE 18
	//-----------------------------------------------------------------------------
	#define	RINGO_SINGOLO1 1		// null ring null
	#define	RINGO_SINGOLO2 2	
	#define	RINGO_DOPPIO1 3		// null ring null
	#define	RINGO_DOPPIO2 4	
	#define	RINGO_OK 5	
	#define	RINGO_POINT1 6
	#define	RINGO_POINT2 7
	#define	RINGO_POINTD1 8
	#define	RINGO_POINTD2 9
	//-----------------------------------------------------------------------------
	#define	NO_RINGO		0
	#define	RINGO_DEVS		1
	#define	RINGO_MEDIANA	2
	#define	SCARTO_MEDIO	3
	#define	ALTO			4
	#define	BASSO			5
	#define RINGO_DOPPIO	8
	#define RINGO_SINGOLO	9
	#define RINGO_TRIPLO	10
	//-----------------------------------------------------------------------------
	#define CONT 16   // contorni forse geo new
	#define CONT2 17   // contorni forse geo new
	//-----------------------------------------------------------------------------
	#define CANCEL 1		// valore di save policy
	#define OVERWRITE 2		// valore di save policy
	#define RENAME 3		// valore di save policy
	#define NO_SAVE 4		// valore di save policy
	//-----------------------------------------------------------------------------
	#define OUTLIER1 10
	#define OUTLIER2 11
	#define OUTLIER_FALSE 12   // falsi outlier
	//-----------------------------------------------------------------------------
	#define METAL_POINT_UP 13
	#define METAL_POINT_DOWN 14
	#define METAL_POINT 15
	#define MA_LINEAR_CORRECTION 1
	#define MA_EXPONENTIAL_CORRECTION 0
	//-----------------------------------------------------------------------------
	#define REFOBJ_BRIGHT 1
	#define REFOBJ_DARK 2
	//-----------------------------------------------------------------------------
	// halfscan defines
	//-----------------------------------------------------------------------------
	#define		HS_ONE_SIDE		0
	#define		HS_BOTH_SIDES	1
	#define		HS_LEFT			2
	#define		HS_RIGHT		3
	#define		HS_BOTH			4
	#define		HS_LINEAR_CORRECTION	1		// correzione del sinos per halfscan
	#define		HS_COSQUADRO_CORRECTION	0		// correzione del sinos per halfscan
	#define		HS_WEIGHTED_CORRECTION	2		// correzione del sinos per halfscan
	//-----------------------------------------------------------------------------
	#define DO_NOT_DEALLOCATE_DATASET2 0
	#define DEALLOCATE_DATASET2 1
	//-----------------------------------------------------------------------------
	// define di ricostruzione
	//-----------------------------------------------------------------------------
	#define PARALLEL_BEAM 		1
	#define FAN_BEAM 			2
	#define CONE_BEAM 			3
	#define EQUISPACED_RAYS		4   // peso W1 per FAN BEAM equispaziati
	#define EQUIANGULAR_RAYS	5   // peso W1 per FAN BEAM equiangolari
	#define FELDKAMP_CONEBEAM	6	// peso W1 per CONE BEAM FELDKAMP
	#define STEP				1.	// step di ogni quante proiezioni mostriamo la preview di reconstruction
	#define CUTSLICE			0	// di quanti pixel ridurre il cerchio di ricostruzione: xstart xend ystart yend
	//-----------------------------------------------------------------------------
	// define dei messaggi paralleli
	//-----------------------------------------------------------------------------
	#define NOMESSAGE				-1
	#define STOP					0
	#define	MKATENRAD				1
	#define MKSINOS					2
	#define MKRINGO					3
	#define MKOUTLIER				4
	#define MKMETALARTIFACT			5
	#define RECEIVE_DATASET			6
	#define RESET_DATASET			7
	#define STOP_MAKE_SINOS			11
	#define RECEIVE_SINOS_K			12
	#define WHEREAMI				13
	#define RECEIVE_IZERO_DARK		14
	#define REMAKE_ATENRAD			15
	#define STOP_MAKE_ATENRAD		16
	#define RECEIVE_PROJECTION		17
	#define RECEIVE_RANK_SAVE		18
	#define MKFFT_FILTER			19
	#define MKRECONSTRUCT			20   // ricostruzione parallela parallelizzando sugli angoli
	#define MKROTATION				21
	#define MKRECONSTRUCT_NEWMETHOD	22	 // ricostruzione parallela parallelizzando i pixel della slice NON PIù ATTIVO!!!!
	#define MKREVERT_TO_ATENRAD		23
	#define MKOUTLIERNEW			24
	#define MKRINGONEW				25
	//-----------------------------------------------------------------------------
	// messaggi paralleli per sequence utility
	//-----------------------------------------------------------------------------
	#define MKCROP						43
	#define MKTRESHGRAY_U8				44
	#define MKTRESHGRAY_U16				45
	#define MKTRESHGRAY_FLT				46
	#define MKBINNIG					47
	#define MKSTRETCH_U8				48
	#define MKSTRETCH_U16				49
	#define MKSTRETCH_FLT				50
	#define MKTRESHGEO					51
	#define STOP_MKUTIITY				52
	#define CONTINUE_MKUTILITY			53
	#define RECEIVE_IMAGE_K				54
	#define	TAG_RECEIVE_IMAGE_K			55
	#define	TAG_MKUTILIY_DONE			56
	#define TAG_MKUTILIY				57
	#define TAG_REMAKEUTILIY_OR_STOP	58
	//-----------------------------------------------------------------------------
	// define dei tag x i messaggi paralleli
	//-----------------------------------------------------------------------------
	#define	TAG_MKATRAD_DONE			1
	#define TAG_MAKE_ATRAD				2
	#define TAG_WHEREAMI				3
	#define TAG_REMAKE_OR_STOP_ATENRAD	4
	#define TAG_RECEIVE_PROJECTION		5
	#define TAG_RECEIVE_ATENRADTOSAVE	6
	#define TAG_RECEIVE_DATA_ON_ATENRADTOSAVE 7
	#define TAG_RECEIVE_SINOS_K			8
	#define TAG_RECEIVE_SINOSTOSAVE		9
	#define TAG_RECEIVE_DATA_ON_SINOSTOSAVE 10
	#define TAG_MKSINOS_DONE			11
	#define TAG_RECEIVE_H				12
	#define TAG_RECEIVE_Y				13		// ricostruzione parallela parallelizzando i pixel della slice NEWMETHOD
	#define TAG_RECEIVE_SLICE_PART		14
	//-----------------------------------------------------------------------------
	// altri define dei paralleli
	//-----------------------------------------------------------------------------
	#define	RANK0 0
	#define MAX_NODES 99
	//-----------------------------------------------------------------------------
	// definizione mie strutture
	//-----------------------------------------------------------------------------
	struct data_tomo		
	{
		int Number_of_Images; 	  // numero di immagini nella sequenza se =0 NO_IMAGE non ci sono immagini in memoria se =1 SINGLE_IMAGE una sola immagine e non una sequenza
		int Single_Image_inSequence;	// se è TRUE allora anche se è una singola immagine ha un numero altrimenti no
		int File_Type;			  // tipo di file SDT TIF JPG BMP
		int Data_Dim;			  // (1=linea, 2=bidimensionale, 3=tridimensionale)  
		int Data_Type;			  // tipo di immagini se U8, U16 o FLT  (U16=1 --- U8=5  --- FLT=3)
		int Width;				  // numero di colonne = larghezza
		int Height;				  // numero di righe = altezza
		int		current_image;	  // immagine corrente (se num img = 1 allora possibile solo 0) index from 0
		int		TomographicType;  // indica il tipo tomografico, cioè il tipo di immagini. Può essere projection atenrad sinos recobj
		int		StrangeName;	  // indica che il nome è diverso (non projection, atenrad,etc..) ma che è classificato come step tomografico
		int		Seq_start;		  // indica da che numero si parte nella sequenza
		int		Seq_end;		  // indica fino a che numero si arriva nella sequenza
		int		Dark;			  // presenza della dark
		int		Izero;			  // presenza della I0
		float pxsize;		  // dimensione del pixel X sul detector
		float pysize;		  // dimensione del pixel Y sul detector
		float rxsize;		  // dimensione del pixel X sull'oggetto
		float rysize;		  // dimensione del pixel Y sull'oggetto
		int min;	      // minimo di immagine U16
		int max;	      // massimo di immagine U16
		float min_f;	      		  // minimo di immagine FLOAT
		float max_f;	      		  // massimo di immagine FLOAT
		char	Path[MAX_STRING];		// percorso del file o della sequenza DIRECTORY
		char	Name[SHORT2_STRING];	// nome del file : se è una sequenza senza "_"
		int		Utility_Step_todo;	// indica se stiamo per fare uno degli step di utility
		int		CT_Step_todo;		// lo step da effettuare
		int		SCT_data_loaded;	// ci sono o no i dati sct?
		float	ringo_percent;		// parametro di ringo
		float	ringo_percent2;		// parametro di ringo 2
		float	outlier_percent;	// parametro per outlier
		int		Seq_min;		// il minimo dei valori su TUTTA la sequenza di U16
		int		Seq_max;		// il massimo dei valori su TUTTA la sequenza di U16
		float		Seq_min_f;		// il minimo dei valori su TUTTA la sequenza di FLOAT
		float		Seq_max_f;		// il massimo dei valori su TUTTA la sequenza di FLOAT
		int		new_dataset;		// mi serve per specificare in OpenSDT se deve allocare o meno evita allocazioni continue
		float  metalartifact_percent;
		float  metalartifact_pixel;
		float	metalartifact_down;
		int		metalartifact_correction;
		int		FFT_do_invFFT;   // se è TRUE (come di default) dopo fa anche l'antitrasformata se no no
		int		FFT_do_FeldkampWeight;  // se è TRUE applica anche il peso per Feldkamp (necessario per ricostruire)
	// ------------------------------------------------------------------------------------------------------
	// qualsiasi cambiamento nei dati SCT prevede che si definisca la nuova variabile in dataset (global.h), 
	// che venga letta in OpenSCT e salvata in SaveSCT e poi ovviamente usata dove serve!!!!!!!!
	// ------------------------------------------------------------------------------------------------------
		float	arange;		// SCT: range di acquisizione (es. 360 -360 180 ...)
		int		nslices;	// SCT: è il numero di righe del rivelatore ovvero di slices, cioè l'altezza del rivelatore o di projection o di atenrad o il numero di sinogrammi
		int		nangles;	// SCT: è il numero di angoli cioè il numero di step angolari ripresi 
		int		nrays;		// SCT: è il numero di colonne cioè la largehzza del rivelatore deve essere SEMPRE uguale a dataset.width
		float		sdd;		// SCT: è la distanza frala sorgente e rivelatore   mm
		float		odd;		// SCT: è la distanza frala oggetto e rivelatore (mm)
		float		sod;		// SCT: è la distanza frala sorgente e oggetto (mm)
		float	pxcenter;		// SCT: centro di rotazione in X		
		float	pzcenter;		// SCT: centro di rotazione in Y	
		float	xsource;	// mia considerazione, centro della sorgente si calcola da IZERO
		float	zsource;	// mia considerazione, centro della sorgente si calcola da IZERO
		int		rec_slice_circle_present;	// se è TRUE sono presenti in memoria x_start,x_end per la ricostruzione dentro al cerchio della slice
		int		rec_slice_circle_use;		// se è TRUE vanno usati x_start,x_end per la ricostruzione dentro al cerchio della slice
		char	recobj_name[SHORT2_STRING];	// nome del file della recobj: se è una sequenza senza "_"
		int		rec_try_pxcenter;      // se è TRUE allora il numero della sequenza corriponde a un centro diverso
		float	pxcenter_start;
		float	pxcenter_end;
		float	pxcenter_step;		
		int		localCT_status;		// se è TRUE si vuole ricostruire in localtomography
		int		localCT_x0;			// parametri localtomography
		int		localCT_y0;			// parametri localtomography
		int		localCT_x1;			// parametri localtomography
		int		localCT_y1;			// parametri localtomography
		int		use_pxcenter_vector;	// se si usa o meno il vettore dei centri
		float	pxcenter_vector_valS0;	// valore del centro nel sinogramma S0
		float	pxcenter_vector_valS1;	// valore del centro nel sinogramma S1
		int		pxcenter_vector_S0;		// numero del sinogramma S0
		int		pxcenter_vector_S1;		// numero del sinogramma S1
		float	ringonew_T_singoli;		// fattore moltiplicativo per determinare la soglia per ring singoli nell'algoritmo RingoNew
		float	ringonew_T_doppi;		// fattore moltiplicativo per determinare la soglia per ring doppi nell'algoritmo RingoNew
		float	ringonew_T_tripli;		// fattore moltiplicativo per determinare la soglia per ring tripli nell'algoritmo RingoNew
		int		halfscan;			// HALFSCAN: si ricostruisce in halfscan
		int		hs_left;			// HALFSCAN: allargamento a sinistra
		int		hs_right;			// HALFSCAN: allargamento a destra
		int		hs_corrtype;		// HALFSCAN: tipo di correzione
		int		REC_geometry;		// tipo di geometria: FAN_BEAM oppure CONE_BEAM
// !!!!! qualsiasi cosa aggiungi qui devi inizializzarla in dataset.c e iniviarla e riceverla nei due parallel !!!!!
		
		unsigned short 	*Image_U16;				// puntatore all'immagine U16
		float		  	*Image_FL;	 			// puntatore all'immagine FLOAT
		unsigned char	*Image_U8;				// puntatore all'immagine U8
		unsigned short 	*Image_U16_stretch;		// puntatore all'immagine U16 stretchata
		float		  	*Image_FL_stretch; 		// puntatore all'immagine FLOAT stretchata
		unsigned char	*Image_U8_stretch;		// puntatore all'immagine U8 stretchata
		unsigned short 	*Image_U16_Dark;		// puntatore all'immagine U16 DARK
		unsigned short 	*Image_U16_Izero;		// puntatore all'immagine U16 I0
		float		   	*Image_FL_Dark;			// puntatore all'immagine FLOAT DARK
		float		   	*Image_FL_Izero;		// puntatore all'immagine FLOAT I0
		float			*xfslice_start;			// puntatore xstart alle coordinate del cerchio della slice
		float			*xfslice_end;			// puntatore xend   alle coordinate del cerchio della slice
	};


//-----------------------------------------------------------------------------





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
// GLOBAL.H
// mio file di intestazione con variabili globali, define, funzioni
//-----------------------------------------------------------------------------
#ifndef GLOBAL_HEADER
	#define GLOBAL_HEADER
	//-----------------------------------------------------------------------------
	// Include CVI e C
	//-----------------------------------------------------------------------------
	#include <stdlib.h>
	#include <math.h>
	#include <stdio.h>
	#include <float.h>
	#include <string.h>

#include <time.h>



	//-----------------------------------------------------------------------------
	//  intestazioni di file sorgente
	//-----------------------------------------------------------------------------
	#include "REC_FUNCTIONS.h"
	//-----------------------------------------------------------------------------
	// variabili globali
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	extern float	MAX_FLOAT;		
	extern float	MIN_FLOAT;
	//-----------------------------------------------------------------------------
	extern int 				ImageControl_ID;		// handle del controllo immagine nel p_main
	extern struct data_tomo DataSet; 				// variabile globale per i parametri in memoria
	extern struct data_tomo DataSet2; 				// variabile globale per la seconda sequenza di parametri in memoria
	extern int 		   		current_bitmapID;  		// handle della bitmap da visualizzare
	extern int 		   		preview_bitmapID; 	// handle della bitmap preview da visualizzare
	extern unsigned char	*bits_preview2;  		// bits dell'immagine PREVIEW2 da visualizzare
	extern int	 		   	preview_bitmapID2;  		// handle della bitmap2 preview da visualizzare
	extern int				*Hist;	    			// puntatore dell'istogramma
	extern int			   	*Hist_height; 			// puntatore ai valori verticali dell'istogramma
	extern int				*Hist2;	    			// puntatore dell'istogramma
	extern int			   	*Hist2_height; 			// puntatore ai valori verticali dell'istogramma
	extern unsigned char	*bits;  				// bits dell'immagine da visualizzare
	extern unsigned char  	*bits_preview;  		// bits dell'immagine PREVIEW da visualizzare
	extern int				*Led_Array;				// controlli dei led
	extern int				menubar_handle;			// handle del menubar
	extern char				SCT_Data[MAX_STRING];	// tutto il file SCT
	extern int 				IMAGE1_LEFT;			
	extern int 				IMAGE1_TOP;
	extern int				*Min_array;				// array dei minimi di tutta la sequenza
	extern int				*Max_array;				// array dei massimi di tutta la sequenza
	extern float			*Minf_array;				// array dei minimi di tutta la sequenza
	extern float			*Maxf_array;				// array dei massimi di tutta la sequenza
	
	//-----------------------------------------------------------------------------
	extern int status_ANT;
	extern unsigned char	*bits_ant;
	extern int *old_step;
	extern int *new_step;
	extern int *xant,*yant;
	
	//-----------------------------------------------------------------------------
	extern int			n_pxcenter;   // mi serve per dire che indice devo salvare quando trycenter è attivo
	extern float		*try_pxcenter;	// vettore dei punti di centri provati
	//-----------------------------------------------------------------------------
#endif  /* GLOBAL_HEADER */
	
	
