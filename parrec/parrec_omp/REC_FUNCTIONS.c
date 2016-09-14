//-------------------------------------------------------------------------
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
// REC_FUNCTION.c
// mio file con le funzioni di ricostruzione per parallelizzazione
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
#include "Global.h"
#include <omp.h>

//-----------------------------------------------------------------------------
// variabili global per ricostruzione
//-----------------------------------------------------------------------------
double  *filter=NULL;						// FILTRO
double  *row_re=NULL,*row_im=NULL;			// le due parti reale e immaginaria della riga
float   *slice=NULL,*slice_rank=NULL;		// le slice da allocare e ricostruire
double	*detector_weighing=NULL;			// il peso del detector
int		h_start=0,h_end=0;					// la suddivisione degli angoli in ranks
int		*x_start=NULL,*x_end=NULL;			// vettori con i punti x di partenza e arrivo in funzione di y
double	*row_modified=NULL;					// riga modificata per halfscan si allocherà solo se servirà
//-----------------------------------------------------------------------------
// per allocare row e filter è necessario calcolare la lunghezza degli array nel dominio di fourier
// è necessario avere almeno il doppio meno uno della dimensione dei punti da elaborare
// inoltre se si usa un array con lunghezza potenza di 2 la FFT lavora meglio per cui:
//-----------------------------------------------------------------------------
int CalcTheMinimumNumberPowerOfTwoAfter(int dim) {

	int	N;
	float N_float;

	// si calcola il logaritmo in base e (naturale)
	N_float = log( (float)dim );

	// si trasforma in base 2
	N_float = N_float / log(2.);

	// si somma uno per arrotondare in eccesso (anche se venisse 10.2 vorremmo usare 11 perchè 10 non basterebbe!)
	N_float += 1.;

	// N è l'intero della potenza di due dell'intero di N_float
	N_float = (float) ( (int)N_float );
	N = (int) ( pow(2.,N_float) );

	return N;

}
//-----------------------------------------------------------------------------
// pesatura per geometria di scanning: nel fan beam la pesatura dei punti del 
// rivelatore va fatta in funzione della geometria e si distingue in 
// rivelatori con raggi equispatiati o equiangolari (vedere il Kak)
// ora calcola il peso anche per feldkamp!!! la differenza sta nel fatto che dw diventa 2D
//-----------------------------------------------------------------------------
int CalculateWeighing_OnTheBaseOfCTGeometry(int detector_geometry,int FFT_start,int nrays,float sod,float pxsize,float pxcenter) {

	switch(detector_geometry) {

		// pesatura per equispaziati SOLO FAN BEAM ---W1--- N.B. non dipende da h quindi facciamo i calcoli fuori dal ciclo!!!
		// tutti gli array dei valori che vanno nel dominio delle frequenze partono da 0 a FFT_Lenght
		// ma noi abbiamo dei valori del detector a partire da FFT_start
		// quindi k è l'indice degli array nel dominio delle frequenze e va da zero a FFT_lenght
		// mentre i è l'indice dei valori del detector e si trova a i=k-FFT_start e va da 0 a width-1
		// n è l'indice definito in modo che sia zero in pxcenter n=k-FFT_start-pxcenter
		// n va da -pxcenter a width - pxcenter ed è 0 in i=pxcenter e k=FFT_start+pxcenter
		// pesatura per equispaziati SOLO FAN BEAM è il valore di sod /sqrt(sod*sod+n*n*pxsize*pxsize) perchè EQUISPAZIATI 
		case(EQUISPACED_RAYS) : {

			double sod_quadro, pxsize_quadro,n;
			int k,width;

			width=DataSet.Width;
			sod_quadro = (double)sod*(double)sod;   
			pxsize_quadro = (double)pxsize*(double)pxsize;
			for(k=FFT_start; k<FFT_start+width; k++) {  // qui ho messo width al posto di nrays per halfscan
				n = (double)k-(double)FFT_start-(double)pxcenter; // l'ho trasformato in double
				detector_weighing[k-FFT_start] = sod / ( sqrt(sod_quadro+n*n*pxsize_quadro));
				// diminuiamo la memoria da allocare (dato che si tratta di un double) allocando width (e non FFT_lenght) e spostando poi k da FFT_start
			}
			return 0;
		}

		case(EQUIANGULAR_RAYS) : {
			
			// da sviluppare
			return 0;
		}

		// in questo caso il detector weighing è bidimensionale, dipende da x e da y del detector in funzione di pxcenter e pzcenter
		// ma non è dipendente da teta (angolo di rotazione) perchè sorgente e rivelatore si muovono solidalmente
		case(FELDKAMP_CONEBEAM) : {
			
			float	sod_f=(float)sod,pzcenter=DataSet.pzcenter,pysize=DataSet.pysize;
			float	sod_fq,kfq,jfq;
			int		j,k;
			int		width=DataSet.Width,height=DataSet.Height;
		
			// variabili per salvare dw
			//float	*prova=NULL;
			//char	fname[MAX_STRING];

			//prova = (float *)malloc(sizeof(float)*width*height+1);

			sod_fq = sod_f*sod_f; // ci serve per Feldkamp weight

			//printf("\n pxc=%f pyc=%f pxs=%f pys=%f",pxcenter,pzcenter,pxsize,pysize);
			//printf("\n sod_f=%f ",sod_f);

			for(j=0; j<height; j++) {
				// se si filtra per Feldkamp cambia il peso e va applicato nel ciclo perchè oltre che da x dipende anche da y
				jfq = (float)(j)-pzcenter;
				jfq = jfq*pysize;
				jfq = jfq*jfq;
				for(k=FFT_start; k<FFT_start+width; k++) {
					kfq = (float)(k)-(float)FFT_start-pxcenter;
					kfq = kfq*pxsize;
					kfq = kfq*kfq;
					detector_weighing[j*width+k-FFT_start] = (sod_f/sqrt(sod_fq+kfq+jfq));
					// diminuiamo la memoria da allocare (dato che si tratta di un double) allocando width (e non FFT_lenght) e spostando poi k da FFT_start
				}
			}

			return 0;
		}

	}

	return 0;
}
//-----------------------------------------------------------------------------
// calcolo del filtro da applicare nello spazio di Fourier a partire da RAMP e con BUTTER
// bisogna passargli il tipo di W2 che vogliamo applicare
//-----------------------------------------------------------------------------
int CalculateFourierFilter(int FFT_lenght, int FFT_start, double W2) {

	int k;
	double cutoff,delta_f,freq,butter;

	// inizializzazione: zero-padding
	for(k=0; k<FFT_lenght; k++) 
		filter[k]=0.;

	// in TEORIA: si aumenta linearmente il valore del filtro a partire dal centro in cui è zero // cutoff è per definizione 0.5
	// in realtà poichè la funzione FFT trasforma nello spazio di fourier con la frequenza di nyquist all'esterno e non al centro
	// invertiamo la definizione del filtro piuttosto che invertire due volte la projection nello spazio di fouirer
	// per tutti i dettagli vedi mio quaderno, il filtro qui implementato è 0.5 al centro (dove c'è Nyquist) e va a zero ai bordi

	// filtro di tipo ramp 
	cutoff = 0.5*W2;   // il filtro deve andare a 1/2pxsize in più c'è W2 che però ho messo a 1. altrimenti diventava pesantissimo!!!
	delta_f = 1.0 / (float)FFT_lenght;
	freq = 0.0;
	filter[0] = 0.;
	filter[FFT_lenght-1] = 0.;
	filter[FFT_start] = cutoff;

	// riempiamo il filtro di ramp
	for (k=1; k<=(FFT_start); k++) {  
		freq += delta_f;
		if (freq <= cutoff)										
			filter[k] = filter[FFT_lenght-k] = (double)k *delta_f*W2; // questo peso ci serve perchè il filtro vada anzichè a 0.5 a 0.5*1/pxsize
		else
			filter[k] = filter[FFT_lenght-k] = 0.0;
	}

	// filtro butter
	freq=0.;
	for (k=1; k<=(FFT_start); k++) {  
		freq += delta_f;
	    butter = 1.0 / sqrt(1.0 + pow((double)(freq/cutoff),10.));
	    if (butter<=0.0)
	       filter[k] = filter[FFT_lenght-k] = 0.0;
	    else
	       filter[k] = filter[FFT_lenght-k] = filter[k] * butter;
	}

	return 0;
}
//---------------------------------------------------------
// calcolo del peso per equispaziati U --- W3 ---
// questa funzione alloca l'array Uquadro_teta_zero e lo calcola
//---------------------------------------------------------
int	CalcolateWeighing_EquispacedRays(float sod,float rxsize) {

	int		y=0,width=DataSet.Width;
	float	val=0.,U,Uquadro;
	float	l;
	float	*Uquadro_teta_zero=NULL;			// il peso Uquadro

	// si alloca U
	free(Uquadro_teta_zero);		Uquadro_teta_zero = NULL;		Uquadro_teta_zero= (float *) malloc ( sizeof(float)*width+1);

	// si calcola l
	l = rxsize*((float)width/2.);

	// si calcola U di teta=0
	for(y=0; y<width; y++) {
		val = sod -( (float)y*rxsize-l);
		U = sod/val;
		Uquadro = U * U;
		Uquadro_teta_zero[y] = Uquadro;
	}

	printf("\n y\tU quadro di teta zero");
	for(y=0; y<width; y++) 
		printf("\n%d\t%f",y,Uquadro_teta_zero[y]);
	fflush(0);

	free(Uquadro_teta_zero);		Uquadro_teta_zero = NULL;	

	return 0;
}
//-----------------------------------------------------------------------------
// funzione per retroproiezione filtrata CON constrain circle SOLO MATEMATICA
// N.B. PRIMA DI ACCEDERE A QUESTA FUNZIONE BISOGNA ALLOCARE SLICE ROW RE ROW IM E FILTER
// INOLTRE è NECESSARIO AVER CARICATO IN MEMORIA IMAGE DA ELABORARE E AVER CALCOLATO FILTER E PESI
// // **** 2013 07 19 ANCHE LOCAL CT **** SENZA HALFSCAN
//-----------------------------------------------------------------------------
int BackProjectFilteredDataOnSlice(int FFT_lenght,int FFT_start,float xs,float ys,float	xd0,float yd0,float	xd1,float yd1) {

	// indici
	int		h,k,riga;
	int		x,y;
	float	xf,yf;
	// dati da SCT
	int		nrays,width;
	float	arange,nangles;
	float	pxsize,rxsize,sod;
	float	wmezzif;
	// variabili per stare dentro al cerchio e localtomography
	float	r=0.,fanangle=0.,l=0.;
	float	r_pixel=0.;
	int		y_start=0,y_end=0;
	// variabili per rotazione
	float	teta,teta_rad,teta_step,sin_teta,cos_teta;
	float	xs_rot,ys_rot;
	float	xd0_rot,yd0_rot,xd1_rot,yd1_rot;
	// variabili per rette
	float	detector_verticale=FALSE;   // se è TRUE vuol dire che B=0 e y=xd_rot
	float	ray_verticale = FALSE;		// se è TRUE vuol dire che b=0 e y=xf
	float	a,b;						// coefficienti della retta tra sorgente ruotata e punto xf yf
	float	A,B;						// coefficienti della retta del detector ruotato (da D0 a D1)
	float	xT,yT;						// coordinate del punto di intersezione fra retta sorgente e punto e retta detector
	float	xd_min,xd_max,yd_min,yd_max; // ci servono per vedere se c'è l'intersezione nel range del detector
	// variabili per interpolazione
	float	i_float,val;    // necessario per interpolare sul detector
	int		ia,ib;			// i due punti da interpolare
	float	dA,dB;			// distanze dal punto da interpolare
	float	d;				// distanza dal punto di intersezione xT e yT da D0 necessario per interpolazione
	// variabili per pesatura 
	double	U,U_quadro;						// pesatura per EQUISPAZIATI: rapporto fra sod e distanza della proiezione del punto della slice 
	double	xs_rot_quadro,ys_rot_quadro;	// sull'asse che va dalla sorgente all'origine
	float	val_pesato_U;
	double	sod_quadro;

	int num_thread;
	double start_time, run_time;

// per salvare slice parziale
//	char	path[MAX_STRING],path_SDT[MAX_STRING],path_SPR[MAX_STRING],txt1[SHORT2_STRING];
//	char	txt[MEDIUM_STRING];


	// dati SCT
	width = DataSet.Width;
	wmezzif = ( (float)width)/2.;   // width mezzif è width trasformato in float e diviso 2
	nrays	= DataSet.nrays;		// vale width a meno che non si faccia mezzo giro HALFSCAN
	sod		= DataSet.sod;
	pxsize	= DataSet.pxsize;					
	rxsize	= DataSet.rxsize;

	// N.B. nel mio sistema di riferimento la sorgente è in alto, il rivelatore è in basso
	// in questo modo sono consistente con me stessa ma rispetto a imgrec ho dovuto inserire il meno davanti a arange
	arange =  - DataSet.arange; 
	nangles = (float) DataSet.nangles;
	teta_step = arange / nangles; // teta step deve essere negativo altrimenti non cambierebbe il verso di rotazione

	// calcoliamo r e l in funzione di width e non nrays perchè nrays potrebbe essere metà
	l = wmezzif*rxsize; // l è la metà del lato del quadrato della slice in mm
	// qui non serve più il calcolo di fanangle e r!!!!!!!!!!!!
	// calcoliamo fanangle e r con - shift 
	Calculate_r_fi(FAN_BEAM,&r,&fanangle);
	r_pixel = r/rxsize;
	//printf("\nNUOVI  VALORI: fanangle=%f   gradi=%f   r=%f    r_pixel=%f",fanangle,fanangle*360./(2.*PIGRECO),r,r_pixel);
	// per stare dentro al cerchio r
	if(DataSet.localCT_status == FALSE ) {
		y_start = (int)(wmezzif - r_pixel)+CUTSLICE;
		y_end   = (int)(wmezzif + r_pixel)-CUTSLICE;
		// controllo y_start e y_end nei limiti
		y_start = (y_start >= 0    ) ? y_start : 0 ;
		y_start = (y_start < width ) ? y_start : width-1 ;
		y_end   = (y_end >= 0    ) ? y_end : 0 ;
		y_end   = (y_end < width ) ? y_end : width-1 ;
	} else {
		// localtomography
		y_start = DataSet.localCT_y0;
		y_end   = DataSet.localCT_y1;
	}
	// determinazione di x_start e x_end SENZA e CON localtomography
	// cerchio r
	// questa funzione alloca i vettori x_start e x_end e li riempe di valori
	ConstrainXstartXendIntoCircle_r(r,l,rxsize,width,y_start,y_end);


	// calcolo di sod quadro per pesatura su U
	sod_quadro = (double)sod*(double)sod; 


	for (num_thread=4;num_thread<=4;num_thread++){

	omp_set_num_threads(num_thread);
        start_time = omp_get_wtime();

	printf("%d threads\n",num_thread);

        // zero padding di slice: anche qui si usa width perchè nrays potrebbe essere la metà di width
        k=0;

 //      for(x=0; x<width; x++)
   //             for(y=0; y<width; y++)
     //                   slice[y*width+x]=0.;


	#pragma omp parallel for schedule(static)
	  for (x = 0; x < width*width; x++) slice[x] = 0.; 
		

	//---------------------------------------------------------
	// ciclo su teta in funzione di h
	//---------------------------------------------------------
	for(h=0; h<(int)nangles; h++) {  
	        teta = (float)h*teta_step;
		//	printf("%cteta %f of %f (h=%d of %d)",0x0D,teta,teta_step*(float)(nangles),h,(int)nangles); fflush(0);

		// teta è in gradi x trasformarlo in radianti si deve applicare la	teta rad = 2 pigreco * teta grad / 360.
		teta_rad = 2.* PIGRECO * teta / 360.;

		// zero padding di re e im
		#pragma omp parallel for schedule(static)
		for(k=0; k<FFT_lenght; k++) {
			row_re[k]=0.;
			row_im[k]=0.;
		}

		riga=h*width;


		// inseriamo dati di immagine in row_re centrando per avere negativi per fare FFT
		#pragma omp parallel for schedule(static)
		for(k=FFT_start; k<FFT_start+width; k++) {// qui ci va width e non nrays altrimenti in halfscan non fa bene il filtro
			row_re[k]=(double)DataSet.Image_FL[riga+k-FFT_start];
			// pesatura per equispaziati SOLO FAN BEAM ---W1---
			row_re[k] *= detector_weighing[k-FFT_start];
		}	
	
		// fast fourier transform
		FFT(row_re,row_im,FFT_lenght);

		// convoluzione con il filtro 
		#pragma omp parallel for schedule(static)
		for(k=0; k<FFT_lenght; k++){
			row_re[k]=row_re[k]*filter[k];
			row_im[k]=row_im[k]*filter[k];
		}		

		// inverse fft: inverse fast fourier transform
		InvFFT(row_re,row_im,FFT_lenght);

		// calcolo di sin teta e cos teta per non ripeterlo ogni volta
		sin_teta = sin(teta_rad);
		cos_teta = cos(teta_rad);
		
		// coordinate ruotate della sorgente rispetto all'origine 0,0
		xs_rot = xs*cos_teta-ys*sin_teta;
		ys_rot = xs*sin_teta+ys*cos_teta;

		// coordinate ruotate della sorgente al quadrato: ci serve per calcolare W3, le calcoliamo qui perchè non dipendono da i
		xs_rot_quadro = xs_rot*xs_rot;
		ys_rot_quadro = ys_rot*ys_rot;

		// calcoliamo d0 e d1 rot per tracciare le rette e interpolare sulla slice
		xd0_rot = xd0*cos_teta-yd0*sin_teta;
		yd0_rot = xd0*sin_teta+yd0*cos_teta;
		xd1_rot = xd1*cos_teta-yd1*sin_teta;	
		yd1_rot = xd1*sin_teta+yd1*cos_teta;


		xd_min = ( xd0_rot < xd1_rot ) ? xd0_rot : xd1_rot;
		xd_max = ( xd0_rot > xd1_rot ) ? xd0_rot : xd1_rot;
		yd_min = ( yd0_rot < yd1_rot ) ? yd0_rot : yd1_rot;
		yd_max = ( yd0_rot > yd1_rot ) ? yd0_rot : yd1_rot;

		// i coefficienti della retta del detector possono essere calcolati fuori perchè dipendono solo da teta
		// stessa cosa per xd,yd sia min che max				// retta y=Ax+B
		// trattiamo separatamente i due casi di detector non verticale e detector verticale
	 
		// Caso detector non verticale: si tracciano le rette

		if( xd0_rot != xd1_rot ) {
		  // detector_verticale=FALSE;   // si calcolano normalemente A e B
		  A = (yd0_rot-yd1_rot)/(xd0_rot-xd1_rot);
		  B = yd0_rot - xd0_rot*A;

		  // ---------------------------------------------
		  // ciclo su x e y INTERI per tracciare le rette
		  // ---------------------------------------------
		 #pragma omp parallel for collapse(2) private (xf,yf,xT,yT,a,b,ia,ib,dB,dA,d,i_float,val,U,U_quadro,val_pesato_U)
		  for(y=y_start; y<=y_end; y++) {
		        //                  yf = (float)y*rxsize -l;        // qui andrebbe + rxsize/2

                 // #pragma omp parallel for private (xf,xT,yT,a,b,ia,ib,dB,dA,d,i_float,val,U,U_quadro,val_pesato_U)
  
			for(x=y_start; x<=y_end; x++) {
//			  for(x=x_start[y]; x<=x_end[y]; x++) {
	    
			// calcolo di yf
			yf = (float)y*rxsize -l;	// qui andrebbe + rxsize/2
		 	// calcolo di xf
		      	xf = (float)x*rxsize -l;		// qui andrebbe + rxsize/2
		      
		      // ------ si tracciano le rette ------ //
		    
		      // calcoliamo la retta per il punto xf,yf xs_rot,ys_rot
		      // calcoliamo l'intersezione fra le due rette xT,yT, nell'ipotesi di detector non verticale
		      // la retta del detector è y=Ax+B
		      if( xs_rot == xf ) {
			//ray_verticale=TRUE;   // se è TRUE vuol dire che b=0 e y=xs_rot
			xT = xs_rot;
			yT = xT*A+B;
		      }	    else {
			//ray_verticale=FALSE;   // si calcolano normalmente a e b
			a = (ys_rot-yf)/(xs_rot-xf);
			b = yf - xf*a;
			xT = (b-B)/(A-a);
			yT = A*xT+B;
		      }
		    		    
		      // ora abbiamo calcolato xT,yT vediamo se si trova all'interno del detector
		      if( xT >= xd_min && xT <= xd_max && yT >= yd_min && yT <= yd_max) {
			// c'è intersezione calcoliamo i_float
			d = (xT-xd0_rot)*(xT-xd0_rot) + (yT-yd0_rot)*(yT-yd0_rot);
			d = sqrt(d);
			i_float = d/pxsize;
			
			// interpoliamo
			ia = (int) i_float;
			ib = ia+1;
			dB = i_float - (float)ia;
			dA = 1. - dB;
			
			// calcoliamo il valore interpolato da associare a questo punto
			// qui eventualmente saltare se //						val =0.;
			if( DataSet.Image_FL[riga+ia] > 0. && ib < nrays && DataSet.Image_FL[riga+ib] > 0. ) {
			  val = row_re[ia+FFT_start]*dA; 
			  val += row_re[ib+FFT_start]*dB; 
			  
			  // peso per la distanza dalla proiezione sulla retta della sorgente ---W3--- EQUISPAZIATI
			  U = (sod_quadro+xs_rot_quadro+ys_rot_quadro-2.*xs_rot*xf-2.*ys_rot*yf)/sod_quadro;  
			  U_quadro = U * U; 
			  val_pesato_U = val/U_quadro;
			  // fattore 1/N della FFT
 			  val_pesato_U /= (float) FFT_lenght;


			  
			  // associazione del punto nella slice_rank
			  slice[y*width+x] += val_pesato_U;
			}											
		      } // fine di controllo su xT e yT
		    } // fine di x
		  } // fine di y
		  
		  // Caso detector verticale: si tracciano le rette
		} else //if( xd0_rot == xd1_rot ), cioe' se il detector e' verticale
		  {
		    //detector_verticale=TRUE;   // se è TRUE vuol dire che B=0 e y=xd_rot
		    
		    // ---------------------------------------------
		    // ciclo su x e y INTERI per tracciare le rette
		    // ---------------------------------------------
		    #pragma omp parallel for collapse(2) private (xf, yf,xT,yT,a,b,ia,ib,dB,dA,d,i_float,val,U,U_quadro,val_pesato_U)
		    for(y=y_start; y<=y_end; y++) {
	                    for(x=y_start; x<=y_end; x++) {
				//for(x=x_start[y]; x<=x_end[y]; x++) {
			      	// calcolo di yf
			      	yf = (float)y*rxsize -l;	// qui andrebbe + rxsize/2
			      	// calcolo di xf
				xf = (float)x*rxsize -l;		// qui andrebbe + rxsize/2
		    
			// ------ si tracciano le rette ------ //
			
			// calcoliamo la retta per il punto xf,yf xs_rot,ys_rot
			// se il detector è verticale non può esserlo anche il raggio in quanto per definizione 
			// l'asse fra la sorgente e il detector è perpendicolare al detector e fan angle < 180
			//  ray_verticale=FALSE;   // si calcolano normalemente a e b
			a = (ys_rot-yf)/(xs_rot-xf);
			b = yf - xf*a;
						
			// calcoliamo l'intersezione fra le due rette xT,yT
			xT = xd0_rot;        // se è TRUE vuol dire che B=0 e y=xd_rot
			yT = xT*a+b;		 // se il detector è verticale non può esserlo anche il raggio in quanto per definizione 
			  // l'asse fra la sorgente e il detector è perpendicolare al detector e fan angle < 180
			
			// ora abbiamo calcolato xT,yT vediamo se si trova all'interno del detector
			if( xT >= xd_min && xT <= xd_max && yT >= yd_min && yT <= yd_max) {
			  // c'è intersezione calcoliamo i_float
			  d = (xT-xd0_rot)*(xT-xd0_rot) + (yT-yd0_rot)*(yT-yd0_rot);
			  d = sqrt(d);
			  i_float = d/pxsize;
			  
			  // interpoliamo
			  ia = (int) i_float;
			  ib = ia+1;
			  dB = i_float - (float)ia;
			  dA = 1. - dB;
			  
			  // calcoliamo il valore interpolato da associare a questo punto
			  // qui eventualmente saltare se //						val =0.;
			  if( DataSet.Image_FL[riga+ia] > 0. && ib < nrays && DataSet.Image_FL[riga+ib] > 0. ) {
			    val = row_re[ia+FFT_start]*dA; 
			    val += row_re[ib+FFT_start]*dB; 
			    
			    // peso per la distanza dalla proiezione sulla retta della sorgente ---W3--- EQUISPAZIATI
			    U = (sod_quadro+xs_rot_quadro+ys_rot_quadro-2.*xs_rot*xf-2.*ys_rot*yf)/sod_quadro;  
			    U_quadro = U * U; 
			    val_pesato_U = val/U_quadro;
	                    // fattore 1/N della FFT
                            val_pesato_U /= (float) FFT_lenght;
        
			    // associazione del punto nella slice_rank
			    slice[y*width+x] += val_pesato_U;
			  }											
			} // fine di controllo su xT e yT
		      } // fine di x
		    } // fine di y
		    

		  } // fine di if () else () per discriminare se era il caso di detector verticale
	
		
		
	} // fine di teta in funzione di h
	
	run_time = omp_get_wtime() - start_time;
	printf("%d threads, %f seconds\n",num_thread,run_time);
	} // fine di loop per uno specifico num_thread

	return 0;
}

//-----------------------------------------------------------------------------
// funzione per costringere x_start e x_end a stare dentro al cerchio r
// determinazione di x_start e x_end SENZA e CON localtomography
// questa funzione alloca i vettori x_start e x_end e li riempe di valori
//-----------------------------------------------------------------------------
void ConstrainXstartXendIntoCircle_r(float r,float l, float rxsize,int width,int y_start, int y_end) {

	int		xs,xe;
	float	radicale;
	int		y;
	float	yf;

	// bisogna allocare i due vettori
	if(!x_start)
		free(x_start);
	x_start = NULL;
	if(!x_end)
		free(x_end);
	x_end = NULL;
	x_start = (int *) malloc( sizeof(int)*width+1);
	x_end   = (int *) malloc( sizeof(int)*width+1);

//	printf("\n width=%d height=%d",width,DataSet.Height);
//	printf("\n pxcenter =%f",DataSet.pxcenter);
//	printf("\nr = %f   l=%f   (wmezzi=%d)",r,l,width/2);

	for(y=y_start; y<=y_end; y++) {

		// calcolo di yf
		yf = (float)y*rxsize -l;	// qui andrebbe + rxsize/2

		// cerchio r
		radicale = r*r-yf*yf;
		if( radicale > 0. ) {
			radicale = sqrt(radicale);
			xs = (int)( ( (l-radicale)/rxsize ) )+CUTSLICE;		// qui andrebbe + rxsize/2
			xe = (int)( ( (l+radicale)/rxsize ) )-CUTSLICE;		// qui andrebbe + rxsize/2
		} else {
			xs = (int)( l/rxsize )+CUTSLICE;
			xe = xs-CUTSLICE;
		}
		// controllo xs e xe nei limiti
		xs = (xs >= 0    ) ? xs : 0 ;
		xs = (xs < width ) ? xs : width-1 ;
		xe = (xe >= 0    ) ? xe : 0 ;
		xe = (xe < width ) ? xe : width-1 ;

		if(DataSet.localCT_status == TRUE ) {
			xs = (DataSet.localCT_x0 > xs) ? DataSet.localCT_x0 : xs;		// localtomography
			xe = (DataSet.localCT_x1 < xe) ? DataSet.localCT_x1 : xe;			// localtomography
		}

		// si memorizzano x_start e x_end calcolati
		x_start[y] = xs;
		x_end[y] = xe;

//		printf("\n%d\t%d\t%d",y,xs,xe);

	} // ciclo su y


	return;
}
//-----------------------------------------------------------------------------
// calcolo del peso --- W3 --- vale 1/Uquadro per equispaziati
//-----------------------------------------------------------------------------
int CalculateWeight3(int detector_geometry,float l,double sod_quadro,double xs_rot_quadro,double ys_rot_quadro,float xs_rot,float ys_rot,float *W3_U) {

	int		x,y;
	int		width=DataSet.Width;
	float	xf,yf;
	float	rxsize=DataSet.rxsize;
	float	U,U_quadro,uno_su_U_quadro;

	//printf("\n%s\n",DataSet.Name);
	
	for(y=0; y<width; y++) {

		// calcolo di yf
		yf = (float)y*rxsize -l;	// qui andrebbe + rxsize/2

		//printf("\n");

		// ciclo su x
		for(x=0; x<width; x++) {

			// calcolo di xf
			xf = (float)x*rxsize -l;		// qui andrebbe + rxsize/2

			// peso per la distanza dalla proiezione sulla retta della sorgente ---W3--- EQUISPAZIATI
			U = (sod_quadro+xs_rot_quadro+ys_rot_quadro-2.*xs_rot*xf-2.*ys_rot*yf)/sod_quadro;  
			U_quadro = U * U; 
			uno_su_U_quadro = 1. / U_quadro;

			W3_U[y*width+x]=uno_su_U_quadro ;
			//printf("%f\t",uno_su_U_quadro);
		}
	}


	return 0;
}
//---------------------------------------------------------
// funzione per il calcolo della geometria necessaria alla ricostruzione FBP
//---------------------------------------------------------
void CalcAllGeometryData_FBP(float *xs,float *ys,float *xd0,float *yd0,float *xd1,float *yd1,int *y_start,int *y_end) {

	int		width;
	float	wmezzif;
	float	pxsize,rxsize,sod,sdd,odd;
	int		nrays;
	float	pxcenter,shift,shift_pixel;
	float	arange,teta_step,nangles;
	// variabili per stare dentro al cerchio 
	float	r,fanangle,l;
	float	r_pixel;


	// mettiamo i valori che usiamo in delle variabili locali
	width = DataSet.Width;
	nrays = DataSet.nrays;
	odd = DataSet.odd;
	sod = DataSet.sod;
	sdd = DataSet.sdd;
	pxsize= DataSet.pxsize;					
	rxsize= DataSet.rxsize;
	pxcenter = DataSet.pxcenter;
	// width mezzif è width trasformato in float e diviso 2
	wmezzif = ( (float)width)/2.;
	// se pxcenter non è al centro del detector allora definiamo lo shift:
	shift_pixel = (pxcenter-wmezzif);  // è lo shift in pixel dal centro di width
	shift = (pxcenter-wmezzif)*pxsize; // è lo shift in mm

	// dati per la rotazione
	// N.B. nel mio sistema di riferimento la sorgente è in alto, il rivelatore è in basso
	// in questo modo sono consistente con me stessa e non si formano ondine strane
	// ma rispetto a imgrec ho dovuto inserire il meno davanti a arange
	arange =  - DataSet.arange; 
	nangles = (float) DataSet.nangles;
	teta_step = arange / nangles; // teta step deve essere negativo altrimenti non cambierebbe il verso di rotazione
	// teta_step_rad = 2.* PIGRECO * Absolute(teta_step) / 360.; // teta step rad è usato solo per moltiplicare come peso la proiezione e quindi deve essere positivo
	
	// calcoliamo r e l in funzione di width e non nrays perchè nrays potrebbe essere metà
	l = wmezzif*rxsize; // l è la metà del lato del quadrato della slice in mm

	// qui non serve più il calcolo di fi e r!!!!!!!!!!!!
	// calcoliamo fanangle e r con - shift 
	Calculate_r_fi(FAN_BEAM,&r,&fanangle);
	r_pixel = r/rxsize;
	//printf("\nNUOVI  VALORI: fi=%f   gradi=%f   r=%f    r_pixel=%f",fi,fi*360./(2.*PIGRECO),r,r_pixel);
	// per stare dentro al cerchio r
	if(DataSet.localCT_status == FALSE ) {
		*y_start = (int)(wmezzif - r_pixel)+CUTSLICE;
		*y_end   = (int)(wmezzif + r_pixel)-CUTSLICE;
		// controllo *y_start e *y_end nei limiti
		*y_start = (*y_start >= 0    ) ? *y_start : 0 ;
		*y_start = (*y_start < width ) ? *y_start : width-1 ;
		*y_end   = (*y_end >= 0    ) ? *y_end : 0 ;
		*y_end   = (*y_end < width ) ? *y_end : width-1 ;
	} else {
		// localtomography
		*y_start = DataSet.localCT_y0;
		*y_end   = DataSet.localCT_y1;
	}
	// determinazione di x_start e x_end SENZA e CON localtomography
	// cerchio r
	// questa funzione alloca i vettori x_start e x_end e li riempe di valori
	ConstrainXstartXendIntoCircle_r(r,l,rxsize,width,*y_start,*y_end);


	// calcoliamo le coordinate iniziali della sorgente che poi ruotando varieranno
	*xs = 0; // -shift;  // DARIMETTERE 0; 
	*ys = sod;   

	// coordinate x e y iniziale del detector che poi ruotando varierà 
	// calcoliamo d0 e d1 per tracciare le rette e interpolare sulla slice
	*xd0 = -wmezzif*pxsize-shift;
	*yd0 = - odd;
	*xd1 = +wmezzif*pxsize-shift;     // attenzione qui si presuppone che nrays==width
	*yd1 = - odd;					   // se si fa halfscan qui ci va nrays/2


	return;
}

//---------------------------------------------------------
// unica funzione che in base a DataSet calcola il raggio r e fan angle
// calcoliamo fanangle e r con - shift 
//---------------------------------------------------------
void	Calculate_r_fi(int geometry,float *r,float *fi) {

	float	fanangle,ray;
	float	wmezzif,pxsize,sdd,odd,sod;
	float	pxcenter;
	int		width;

	// dati da DataSet
	width = DataSet.Width;
	wmezzif = (float)width/2.;
	pxsize = DataSet.pxsize;
	sdd = DataSet.sdd;
	odd = DataSet.odd;
	sod = DataSet.sod;
	pxcenter = DataSet.pxcenter;

	if( geometry == FAN_BEAM ) {

		// calcolo dell'angolo di fan
		fanangle = atanf( wmezzif*pxsize/sdd );  // radianti
		ray = ((pxsize*((float)width-pxcenter)/(tanf(fanangle)))-odd) * sin(fanangle);

	} else if ( geometry == CONE_BEAM ) {

		// ---------
		// FELDKAMP
		// ---------
		fanangle = atanf( wmezzif*pxsize/sdd );  // radianti
		ray = sod * sinf(fanangle);

	} else { // PARALLEL BEAM

		fanangle =-99.;
		ray = -99.;
		return;
	}

	// in imgrec a destra r diminuisce a sinistra r aumenta
	// immisioni dati calcolati
	*r = ray;
	*fi = fanangle;

//	printf("\nfi=%f fi(grad)=%f ",fanangle,fanangle*360./(2.*PIGRECO));
//	printf("\nr=%f r_pixel=%f       (r fanbeam=%f  r pixel fanbeam =%f)",ray,ray/pxsize,((pxsize*((float)width-pxcenter)/(tanf(fanangle)))-odd) * sin(fanangle),(((pxsize*((float)width-pxcenter)/(tanf(fanangle)))-odd) * sin(fanangle))/pxsize);
	

	return;

}
//-----------------------------------------------------------------------------
//  valore assoluto di un float
//-----------------------------------------------------------------------------
float Absolute(float val) {

	if( val >= 0. )
		return val;
	else
		return -val;

}

//-----------------------------------------------------------------------------
//	si compone pathname dai dati inseriti in data set
// se n=-1 allora dark se n=-2 allora bak
//-----------------------------------------------------------------------------
void	ComposeFileNameFromDataSet(char *PathName) {
		
	char	txt[MAX_STRING],txt1[SHORT_STRING];
	
	txt[0]='\0';
	txt1[0]='\0';
	
	sprintf(txt,DataSet.Path);
	strcat(txt,DataSet.Name);  // qualsiasi sia il nome della sequenza qui viene immesso

	// se è una proiezione ci può essere una dark o una izero
	if( DataSet.TomographicType == PROJECTION ) {
		
		if(DataSet.current_image==DARK) {
			sprintf(txt1,"drk");
			strcat(txt,txt1); 	
		} else if(DataSet.current_image==IZERO) {
			sprintf(txt1,"bak");
			strcat(txt,txt1); 	
		} else if( DataSet.current_image != NO_SEQ ) {
			sprintf(txt1,"_%d",DataSet.current_image);
			strcat(txt,txt1); 	
		}

	// non è una proiezione ma è riconsciuta come parte di una sequenza con un numero
	} else if( DataSet.current_image != NO_SEQ ) {
		sprintf(txt1,"_%d",DataSet.current_image);
		strcat(txt,txt1);
	} 
		
	if( DataSet.File_Type == SDT)
		strcat(txt,".sdt");
	else if( DataSet.File_Type == JPG)
		strcat(txt,".jpg");
	else if( DataSet.File_Type == BMP)
		strcat(txt,".bmp");
	else if( DataSet.File_Type == TIF)
		strcat(txt,".tif");

	CopyStringToStringAndTerminateIt(PathName,txt);
	
	return;
	
	
}
//-----------------------------------------------------------------------------
// copia string_tocopy in string_target e poi la termina
//-----------------------------------------------------------------------------
void CopyStringToStringAndTerminateIt(char *string_target,char *string_tocopy) {

	int k=0,len=0;
	
	len = (int) strlen(string_tocopy);
	for(k=0; k<len; k++)
		string_target[k]=string_tocopy[k];
	string_target[len]='\0';

	return;
}
//-----------------------------------------------------------------------------
// si calcolano max e min dell'immagine corrente in base al dataset
//-----------------------------------------------------------------------------
void MaxMin(int dataset) {
	
	int x=0, y=0;
	int min=MAX_16BIT,max=ZERO,val=0;
	float min_f=FLT_MAX,max_f=-FLT_MIN,val_f=0;
	int		width,height;
	
	if(dataset == DATASET1) {

		width=DataSet.Width;
		height=DataSet.Height;

		if(DataSet.Data_Type == U16 ) {
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val=(int)DataSet.Image_U16[y*width+x];
					min = ( min < val ) ? min : val;
					max = ( max > val ) ? max : val;
				}
			}
			DataSet.min = min;
			DataSet.max = max;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			if( DataSet.Number_of_Images > SINGLE_IMAGE && Min_array != NULL && Max_array != NULL ) {
				Min_array[DataSet.current_image-DataSet.Seq_start]=min;
				Max_array[DataSet.current_image-DataSet.Seq_start]=max;
			}

		}
	
		if(DataSet.Data_Type == U8) {
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val=(int)DataSet.Image_U8[y*width+x];
					min = ( min < val ) ? min : val;
					max = ( max > val ) ? max : val;
				}
			}
			DataSet.min = min;
			DataSet.max = max;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			if( DataSet.Number_of_Images > SINGLE_IMAGE && Min_array != NULL && Max_array != NULL ) {
				Min_array[DataSet.current_image-DataSet.Seq_start]=min;
				Max_array[DataSet.current_image-DataSet.Seq_start]=max;
			}

		}
		
		if(DataSet.Data_Type == FLT) {	   // rossella float
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val_f=DataSet.Image_FL[y*width+x];
					min_f = ( min_f < val_f ) ? min_f : val_f;
					max_f = ( max_f > val_f ) ? max_f : val_f;
				}
			}
			DataSet.min_f = min_f;
			DataSet.max_f = max_f;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			if( DataSet.Number_of_Images > SINGLE_IMAGE && Minf_array != NULL && Maxf_array != NULL ) {
				Minf_array[DataSet.current_image-DataSet.Seq_start]=min_f;
				Maxf_array[DataSet.current_image-DataSet.Seq_start]=max_f;
			}
		
		}
	// fine di DATSET1
	} else {

		width=DataSet2.Width;
		height=DataSet2.Height;

		if(DataSet2.Data_Type == U16 ) {
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val=(int)DataSet2.Image_U16[y*width+x];
					min = ( min < val ) ? min : val;
					max = ( max > val ) ? max : val;
				}
			}
			DataSet2.min = min;
			DataSet2.max = max;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			//if( DataSet2.Number_of_Images > SINGLE_IMAGE && Min_array != NULL && Max_array != NULL ) {
			//	Min_array[DataSet2.current_image-DataSet2.Seq_start]=min;
			//	Max_array[DataSet2.current_image-DataSet2.Seq_start]=max;
			//}

		}
	
		if(DataSet2.Data_Type == U8) {
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val=(int)DataSet2.Image_U8[y*width+x];
					min = ( min < val ) ? min : val;
					max = ( max > val ) ? max : val;
				}
			}
			DataSet2.min = min;
			DataSet2.max = max;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			//if( DataSet2.Number_of_Images > SINGLE_IMAGE && Min_array != NULL && Max_array != NULL ) {
			//	Min_array[DataSet2.current_image-DataSet2.Seq_start]=min;
			//	Max_array[DataSet2.current_image-DataSet2.Seq_start]=max;
			//}

		}
		
		if(DataSet2.Data_Type == FLT) {	   // rossella float
		
			for(x=0; x<width; x++) {
				for(y=0; y<height; y++) {
					val_f=DataSet2.Image_FL[y*width+x];
					min_f = ( min_f < val_f ) ? min_f : val_f;
					max_f = ( max_f > val_f ) ? max_f : val_f;
				}
			}
			DataSet2.min_f = min_f;
			DataSet2.max_f = max_f;
			// se c'è una sequenza vale la pena inserire i dati calcolati in max min
			//if( DataSet2.Number_of_Images > SINGLE_IMAGE && Minf_array != NULL && Maxf_array != NULL ) {
			//	Minf_array[DataSet2.current_image-DataSet2.Seq_start]=min_f;
			//	Maxf_array[DataSet2.current_image-DataSet2.Seq_start]=max_f;
			//}
		
		}
	} // fine di dataset2
	return;
}

//-----------------------------------------------------------------------------
// si salva l'immagine in input come SDT nel path specificato e con le dimensioni specificate 
// ci sono due funzioni: questa è per i float
//-----------------------------------------------------------------------------
int Save_SDT_GenericFile_FLT(float *image,char *path,int width,int height) {
	
	FILE 	*f;
	size_t	err=0;
	size_t	count_byte=0;
	ssize_t		size=0;
	int		fsize=0;


	// il nome del file da salvare si trova in path
	// si salva il file SDT nel percorso stabilito

	// si scrive il file
	f = fopen(path,"wb");
	count_byte = sizeof(float)*width*height;
	fwrite(image,count_byte,1,f); 
	fclose(f);
	
	return 0;
}
//-----------------------------------------------------------------------------
// si salva l'SPR nel path specificato e con le dimensioni specificate 
// e con i paramteri specificati
// **** size_policy: SI SALVANO ANCHE PSIZE O RSIZE SE COSì SI VUOLE *************
//-----------------------------------------------------------------------------
int Save_SPR_GenericFile(char *path,int width,int height,int IMG_TYPE,int DATA_TYPE,int size_policy) {
	

	FILE 	*f;
	
	// Save spr file
	f = fopen(path, "w");
	if (!f)
	{
		printf("\nWarning! Unable to open spr selected file");
		return -1;
	}

	// si fa l'spr
	fprintf(f,"%d\n",IMG_TYPE);
	fprintf(f,"%d\n",width);
	fprintf(f,"0.000000\n");
	fprintf(f,"0.000000\n");
	fprintf(f,"%d\n",height);
	if( size_policy == SAVE_PSIZE ) {
		fprintf(f,"%.6f\n",DataSet.pxsize);
		fprintf(f,"%.6f\n",DataSet.pysize);
	} else if( size_policy == SAVE_RSIZE ) {
		fprintf(f,"%.6f\n",DataSet.rxsize);
		fprintf(f,"%.6f\n",DataSet.rysize);
	} else {
		fprintf(f,"0.000000\n");
		fprintf(f,"0.000000\n");
	}
	fprintf(f,"%d\n",DATA_TYPE);
	fclose(f);	
	
	return 0;
}
//-----------------------------------------------------------------------------
// FFT by numerical recipes in c
// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
// be an integer power of 2 (this is not checked for!).
//-----------------------------------------------------------------------------
void four1(float data[], unsigned long nn, int isign) {

	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;	//Double precision for the trigonomet
	float tempr,tempi;					//ric recurrences.

	n=nn << 1;
	j=1;
	for (i=1; i<n; i+=2) {		//This is the bit-reversal section of the //routine.
		if (j > i) { 
			SWAP(data[j],data[i]); //Exchange the two complex numbers.
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	// Here begins the Danielson-Lanczos section of the routine.
	mmax=2;
	while (n > mmax) {			//Outer loop executed log2 nn times.
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); //Initialize the trigonometric recurrence.
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {			//Here are the two nested inner loops.
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;				//This is the Danielson-Lanczos formula:
				tempr=wr*data[j]-wi*data[j+1]; 
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; //Trigonometric recurrence.
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
	return;
}
//-----------------------------------------------------------------------------
// fast fourier transform
//-----------------------------------------------------------------------------
void	FFT(double *row_re,double *row_im,int FFT_lenght) {
	
	float	*data=NULL;
	unsigned long nn;
	int		k;

	nn = FFT_lenght;
	data = (float *) malloc ( sizeof(float)*2*nn+1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		data[2*k+1]=(float)row_re[k];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		data[2*k+2]=(float)row_im[k];
	}

	// FFT by numerical recipes in c
	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
	// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
	// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
	// be an integer power of 2 (this is not checked for!).
	four1(data,nn,1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		row_re[k]=(double)data[2*k+1];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		row_im[k]=(double)data[2*k+2];
	}


	free(data);

	return;
}
//-----------------------------------------------------------------------------
// inverse fft: inverse fast fourier transform
//-----------------------------------------------------------------------------
void	InvFFT(double *row_re,double *row_im,int FFT_lenght) {
	
	float	*data=NULL;
	unsigned long nn;
	int		k;

	nn = FFT_lenght;
	data = (float *) malloc ( sizeof(float)*2*nn+1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		data[2*k+1]=(float)row_re[k];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		data[2*k+2]=(float)row_im[k];
	}

	// FFT by numerical recipes in c
	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
	// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
	// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
	// be an integer power of 2 (this is not checked for!).
	four1(data,nn,-1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		row_re[k]=(double)data[2*k+1];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		row_im[k]=(double)data[2*k+2];
	}

	free(data);

	return;
}
//-----------------------------------------------------------------------------
// Stampa del dataset scelto nello standard output
//-----------------------------------------------------------------------------
void PrintData_DataSet(void) {

		
		printf("\n-----------------------------------------------------------------------------"); fflush(0);         			
		printf("\n-----------PRINTING-OF-DATASET1----------------------------------------------"); fflush(0);         			
		printf("\n-----------------------------------------------------------------------------"); fflush(0);         			
		printf("\nDataSet.Number_of_Images 		 	 = %d  ",DataSet.Number_of_Images          ); fflush(0);         			
		printf("\nDataSet.Seq_start 		 		 = %d  ",DataSet.Seq_start 		           ); fflush(0);         			
		printf("\nDataSet.Seq_end 		 			 = %d  ",DataSet.Seq_end 		           ); fflush(0);         			
		printf("\nDataSet.File_Type 		 		 = %d  ",DataSet.File_Type 		           ); fflush(0);         			
		printf("\nDataSet.Data_Dim 		 		 	 = %d  ",DataSet.Data_Dim 		           ); fflush(0);         			
		printf("\nDataSet.Data_Type 		 		 = %d  ",DataSet.Data_Type 		           ); fflush(0);         			
		printf("\nDataSet.Width 			 		 = %d  ",DataSet.Width 			           ); fflush(0);         			
		printf("\nDataSet.Height 			 		 = %d  ",DataSet.Height 			       ); fflush(0);         			
		printf("\nDataSet.pxsize 	 			 = %f  ",DataSet.pxsize 	           ); fflush(0);         			
		printf("\nDataSet.pysize 	 			 = %f  ",DataSet.pysize 	           ); fflush(0);         			
		printf("\nDataSet.rxsize 		 		 = %f  ",DataSet.rxsize 		           ); fflush(0);         			
		printf("\nDataSet.rysize 		 		 = %f  ",DataSet.rysize 		           ); fflush(0);         			
		printf("\nDataSet.current_image	 		 	 = %d  ",DataSet.current_image	           ); fflush(0);         			
		printf("\nDataSet.min 			 			 = %d  ",DataSet.min 			           ); fflush(0);         			
		printf("\nDataSet.max 			 			 = %d  ",DataSet.max 			           ); fflush(0);         			
		printf("\nDataSet.min_f 			 		 = %f  ",DataSet.min_f 			           ); fflush(0);         			
		printf("\nDataSet.max_f 			 		 = %f  ",DataSet.max_f 			           ); fflush(0);         			
		printf("\nDataSet.TomographicType  		 	 = %d  ",DataSet.TomographicType           ); fflush(0);         			
		printf("\nDataSet.StrangeName 	 			 = %d  ",DataSet.StrangeName 	           ); fflush(0);         			
		printf("\nDataSet.Path			 			 = %s  ",DataSet.Path			           ); fflush(0);         			
		printf("\nDataSet.Name			 			 = %s  ",DataSet.Name			           ); fflush(0);         			
		printf("\nDataSet.Dark 			 		 	 = %d  ",DataSet.Dark 			           ); fflush(0);         			
		printf("\nDataSet.Izero 			 		 = %d  ",DataSet.Izero 			           ); fflush(0);         			
		printf("\nDataSet.CT_Step_todo 	 		 	 = %d  ",DataSet.CT_Step_todo 	           ); fflush(0);         			
		printf("\nDataSet.SCT_data_loaded	 		 = %d  ",DataSet.SCT_data_loaded	       ); fflush(0);         			
		printf("\nDataSet.arange 			 		 = %f  ",DataSet.arange 			       ); fflush(0);         			
		printf("\nDataSet.ringo_percent	 		 	 = %f  ",DataSet.ringo_percent	           ); fflush(0);         			
		printf("\nDataSet.ringo_percent2	 		 = %f  ",DataSet.ringo_percent2	           ); fflush(0);         			
		printf("\nDataSet.outlier_percent			 = %f  ",DataSet.outlier_percent		   ); fflush(0); 	   						
		printf("\nDataSet.Utility_Step_todo 		 = %d  ",DataSet.Utility_Step_todo 		   ); fflush(0); 						
		printf("\nDataSet.Seq_min					 = %d  ",DataSet.Seq_min				   ); fflush(0); 	   						
		printf("\nDataSet.Seq_max 					 = %d  ",DataSet.Seq_max 				   ); fflush(0); 		   			
		printf("\nDataSet.Seq_min_f 				 = %f  ",DataSet.Seq_min_f 				   ); fflush(0); 						
		printf("\nDataSet.Seq_max_f 				 = %f  ",DataSet.Seq_max_f 				   ); fflush(0); 						
		printf("\nDataSet.new_dataset 				 = %d  ",DataSet.new_dataset 			   ); fflush(0); 		   			
		printf("\nDataSet.metalartifact_percent 	 = %f  ",DataSet.metalartifact_percent 	   ); fflush(0); 						
		printf("\nDataSet.metalartifact_pixel 		 = %f  ",DataSet.metalartifact_pixel 	   ); fflush(0); 		   			
		printf("\nDataSet.metalartifact_down		 = %f  ",DataSet.metalartifact_down		   ); fflush(0); 						
		printf("\nDataSet.metalartifact_correction	 = %d  ",DataSet.metalartifact_correction  ); fflush(0);  					
		printf("\nDataSet.nslices					 = %d  ",DataSet.nslices				   ); fflush(0); 	   						
		printf("\nDataSet.nangles					 = %d  ",DataSet.nangles				   ); fflush(0); 	   						
		printf("\nDataSet.Image_U16					 = %d  ",DataSet.Image_U16				   ); fflush(0); 	   						
		printf("\nDataSet.Image_FL					 = %d  ",DataSet.Image_FL				   ); fflush(0); 	   						
		printf("\nDataSet.Image_U8					 = %d  ",DataSet.Image_U8				   ); fflush(0); 	   						
		printf("\n-----------------------------------------------------------------------------"); fflush(0);         			
																						

	return;
}


