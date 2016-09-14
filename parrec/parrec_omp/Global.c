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
//-----------------------------------------------------------------------------
#include "Global.h"
//-----------------------------------------------------------------------------
// variabili globali definite come extern nel global.h e qui inizializzate
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// variabili globali definite come extern nel global.h e qui inizializzate
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
float	MAX_FLOAT = 999.0; //max_f;		// rossella rossella ?????? come li scelgo questi due vaolri??                                     		
float	MIN_FLOAT = 0.0; //min_f;		// non possono essere il minimo e il massimo dell'immagine xchè quando stretcho non funziona più   
//-----------------------------------------------------------------------------
struct 			data_tomo DataSet; 			// variabile globale per tutti i parametri in memoria relativi alle proprietà dell'immagine o sequenza
struct 			data_tomo DataSet2; 		// variabile globale per tutti i parametri in memoria relativi alla seconda sequenza
int 		   	current_bitmapID=0;  		// handle della bitmap da visualizzare
int 		   	preview_bitmapID=0;  		// handle della bitmap preview da visualizzare
unsigned char  *bits_preview2=NULL;  		// bits dell'immagine PREVIEW2 da visualizzare
int 		   	preview_bitmapID2=0;  		// handle della bitmap2 preview da visualizzare
int			   *Hist=NULL;	    			// puntatore dell'istogramma
int			   *Hist_height=NULL; 			// puntatore ai valori verticali dell'istogramma
int			   *Hist2=NULL;	    			// puntatore dell'istogramma
int			   *Hist2_height=NULL; 			// puntatore ai valori verticali dell'istogramma
unsigned char  *bits=NULL;  				// bits dell'immagine da visualizzare
unsigned char  *bits_preview=NULL;  		// bits dell'immagine PREVIEW da visualizzare
int 			*Led_Array=NULL;			// controlli dei led per i ranks
int				menubar_handle=0;			// handle del menubar
char			SCT_Data[MAX_STRING]={'\0'};// tutto il file SCT
int 			IMAGE1_LEFT=30;			
int 			IMAGE1_TOP=110;
int				*Min_array=NULL;				// array dei minimi di tutta la sequenza
int				*Max_array=NULL;				// array dei massimi di tutta la sequenza
float			*Minf_array=NULL;				// array dei minimi di tutta la sequenza
float			*Maxf_array=NULL;				// array dei massimi di tutta la sequenza
float			*xfslice_start=NULL;					// puntatore xstart alle coordinate del cerchio della slice
float			*xfslice_end=NULL;					// puntatore xend   alle coordinate del cerchio della slice

//-----------------------------------------------------------------------------
int status_ANT=STOP;
unsigned char	*bits_ant=NULL;
int *old_step=NULL;
int *new_step=NULL;
int *xant=NULL,*yant=NULL;

//-----------------------------------------------------------------------------
int	n_pxcenter=0;   // mi serve per dire che indice devo salvare quando trycenter è attivo
float *try_pxcenter=NULL;	// vettore dei punti di centri provati

