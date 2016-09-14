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
// REC_FUNCTION.c
// mio file con le funzioni di ricostruzione per parallelizzazione
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
#include "Global.h"
#include <stdlib.h>
#include <cufft.h>

//-----------------------------------------------------------------------------
// variabili global per ricostruzione
//-----------------------------------------------------------------------------
float  *filter=NULL;						// FILTRO
float  *row_re=NULL,*row_im=NULL;
int  *x_start=NULL, *x_end=NULL;			// le due parti reale e immaginaria della riga
float   *slice=NULL,*slice_rank=NULL;		// le slice da allocare e ricostruire
float	*detector_weighing=NULL;			// il peso del detector
int		h_start=0,h_end=0;					// la suddivisione degli angoli in ranks
int		*d_x_start=NULL, *d_x_end=NULL;			// vettori con i punti x di partenza e arrivo in funzione di y
float	*row_modified=NULL;					// riga modificata per halfscan si allocherà solo se servirà
//-----------------------------------------------------------------------------
// per allocare row e filter è necessario calcolare la lunghezza degli array nel dominio di fourier
// è necessario avere almeno il doppio meno uno della dimensione dei punti da elaborare
// inoltre se si usa un array con lunghezza potenza di 2 la FFT lavora meglio per cui:
//-----------------------------------------------------------------------------


//--- dichiarazioni kernel ----//
__global__ void zeropaddingslice(float*, int);
__global__ void creation_xStartEnd(int*, int*, float, float, float, int, int, int, int, int, int);
__global__ void zp_d_signal(cufftComplex*, int);
__global__ void initialize_signal(cufftComplex*, float*, int , int , int , int, float*);
__global__ void convolution(cufftComplex*, float*, int);
__global__ void geometry_reconstruction(cufftComplex*, float*, float*, int*, int*, float, float, float, float, float, float, float, float, float, float, float, float, float, float, int, int, int, int, int, int, int, int, int, int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int);
//-------------------------------


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

			sod_fq = sod_f*sod_f; // ci serve per Feldkamp weight

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
int BackProjectFilteredDataOnSlice(int FFT_lenght,int FFT_start,float xs,float ys,float	xd0,float yd0,float	xd1,float yd1,int numBlocks,int threadsPerBlock,int dimBlock) {

	// indici
	int		h,riga;
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
	float	A,B;						// coefficienti della retta del detector ruotato (da D0 a D1)
	float	xd_min,xd_max,yd_min,yd_max; // ci servono per vedere se c'è l'intersezione nel range del detector
	// variabili per interpolazione
	float	xs_rot_quadro,ys_rot_quadro;	// sull'asse che va dalla sorgente all'origine
	//float	val_pesato_U;
	float	sod_quadro;
// per salvare slice parziale
//	char	path[MAX_STRING],path_SDT[MAX_STRING],path_SPR[MAX_STRING],txt1[SHORT2_STRING];
//	char	txt[MEDIUM_STRING];

	// dati SCT
	width = DataSet.Width;
	wmezzif = ((float)width)/2.;   // width mezzif è width trasformato in float e diviso 2
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

	// calcoliamo fanangle e r con - shift 
	Calculate_r_fi(FAN_BEAM,&r,&fanangle);
	r_pixel = r/rxsize;

	// calcolo di sod quadro per pesatura su U
	sod_quadro = sod*sod;

	//host:
	int dominio = width*width;
	//device:
	float *dev_slice;

	cudaMalloc((void**)&dev_slice, (dominio * sizeof(float)+1));//all'indirizzo di dev_slice alloco uno spazio di
																//	w*w*(dim di 1float su quella macchina in bytes)
	printf("\nnumBlocks = %d   threadsPerBlock = %d   dimBlock = %d\n",numBlocks,threadsPerBlock,dimBlock);

	zeropaddingslice<<<numBlocks, threadsPerBlock>>>(dev_slice, dominio);

// creazione di y_start/end e x_start/end[y]
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

	// x_start/end:
	int localCT_status = DataSet.localCT_status;
	int localCT_x0 = DataSet.localCT_x0;
	int localCT_x1 = DataSet.localCT_x1;

	//device:
	// per stare dentro al cerchio r
	if(!d_x_start)
		cudaFree(d_x_start);
	d_x_start = NULL;
	if(!d_x_end)
		cudaFree(d_x_end);
	d_x_end = NULL;

	cudaMalloc((void**)&d_x_start, (width * sizeof(int)+1));
	cudaMalloc((void**)&d_x_end, (width * sizeof(int)+1));
	
	creation_xStartEnd<<<numBlocks, threadsPerBlock>>>(d_x_start, d_x_end, r, l, rxsize, localCT_status, localCT_x0, localCT_x1, width, y_start, y_end);
	

//---------------------------------------------------------
// ciclo su teta in funzione di h
//---------------------------------------------------------

	for(h=0; h<(int)nangles; h++) {
		teta = (float)h*teta_step;

		// teta è in gradi x trasformarlo in radianti si deve applicare la	teta rad = 2 pigreco * teta grad / 360.
		teta_rad = 2.* PIGRECO * teta / 360.;

		cufftComplex *d_signal;
		cudaMalloc((void**)&d_signal, sizeof(cufftComplex)*(FFT_lenght+1));

		zp_d_signal<<<numBlocks, threadsPerBlock>>>(d_signal,FFT_lenght);

		int dimImageFL = DataSet.Width*DataSet.Height;
		riga=h*width;

		float *d_Image_FL;
		cudaMalloc((void**)&d_Image_FL, sizeof(float)*(dimImageFL+1));
		cudaMemcpy(d_Image_FL, DataSet.Image_FL, sizeof(float)*(dimImageFL+1), cudaMemcpyHostToDevice);

		float *d_detector_weighing;
		cudaMalloc((void**)&d_detector_weighing, sizeof(float)*(width+1));
		cudaMemcpy(d_detector_weighing, detector_weighing, sizeof(float)*(width+1), cudaMemcpyHostToDevice);

		cufftComplex *d_signal_image;
		cudaMalloc((void**)&d_signal_image, sizeof(cufftComplex)*(FFT_lenght+1));

		initialize_signal<<<numBlocks, threadsPerBlock>>>(d_signal, d_Image_FL,FFT_start, width, riga, FFT_lenght, d_detector_weighing);

		cudaFree(d_detector_weighing);

		//se la dimensione dell'array row_re (e row_im) è minore di 2^13 allora si fa l'fft del NR, altrimenti la cuFFT
		if (FFT_lenght<pow(2,13)){
			
			cufftComplex *h_signal = (cufftComplex *) malloc(sizeof(cufftComplex)*(FFT_lenght+1));
			cudaMemcpy(h_signal, d_signal,sizeof(cufftComplex)*(FFT_lenght+1), cudaMemcpyDeviceToHost);
			
			for(int ii=0; ii<FFT_lenght; ii++){
				row_re[ii] = h_signal[ii].x;
				row_im[ii] = h_signal[ii].y;
			}
			// fast fourier transform
			FFT(row_re,row_im,FFT_lenght);

			// convoluzione con il filtro 
			for(int k=0; k<FFT_lenght; k++){
				row_re[k]=row_re[k]*filter[k];
				row_im[k]=row_im[k]*filter[k];
			}
				
			// inverse fft: inverse fast fourier transform
			InvFFT(row_re,row_im,FFT_lenght);
			
			for(int ii=0; ii<FFT_lenght; ii++){
				h_signal[ii].x = row_re[ii];
				h_signal[ii].y = row_im[ii];
			}
			
			cudaMemcpy(d_signal, h_signal, sizeof(cufftComplex)*(FFT_lenght+1), cudaMemcpyHostToDevice);
			cudaFree(h_signal);
		}
		else {
			
			//mi preparo i piano per Fourier:
			cufftHandle plan;
			cufftPlan1d(&plan, FFT_lenght, CUFFT_C2C, 1);
	
			// Transformo
			cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE);
	
			float *d_filter;
			cudaMalloc((void**)&d_filter, (FFT_lenght+1) * sizeof(float));
			cudaMemcpy(d_filter, filter, (FFT_lenght+1) * sizeof(float), cudaMemcpyHostToDevice);
	
			//Convoluzione
			convolution<<<numBlocks, threadsPerBlock>>>(d_signal, d_filter, FFT_lenght);
	
			cudaFree(d_filter);
	
			//Antitrasformo
			cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD);
	
			cufftDestroy(plan);
		}

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

		// traccaire le rette:
		int caso = 0;

		// Caso detector non verticale: si tracciano le rette
		if( xd0_rot != xd1_rot ) {

			A = (yd0_rot-yd1_rot)/(xd0_rot-xd1_rot);
			B = yd0_rot - xd0_rot*A;

			caso = 1;
		  		}
		// Caso detector verticale: si tracciano le rette
		else{
			caso = 2;
			A = 0.;
			B = 0.;
		}

		//ora mi calcolo il lato minimo del quardrato entro cui deve essere inscritto il cerchio di ricostruzione (cioè il diametro del cerchio ovviamente):
		int side_square_circle = (y_end+1)-y_start; // visto che è il diametro potrebbe essere anche r_pixel*2--> controllare che valga per tutti i casi!!(localtomography)

		int dimgrid = (side_square_circle+(dimBlock-1))/dimBlock;

		dim3 sizeblock(dimBlock,dimBlock);
		dim3 sizegrid(dimgrid,dimgrid);

		geometry_reconstruction<<<sizegrid, sizeblock>>>(d_signal, dev_slice, d_Image_FL, d_x_start, d_x_end, sod_quadro, pxsize, r, rxsize, xs, ys, teta_rad, xd0, xd1, yd0, yd1, A, B, l, caso, nrays, FFT_lenght, FFT_start, h, localCT_status, localCT_x0, localCT_x1, width, y_start, y_end,sin_teta,cos_teta, xs_rot, ys_rot, xd0_rot, xd1_rot, yd0_rot,yd1_rot, xd_min,xd_max,yd_min, yd_max, xs_rot_quadro,  ys_rot_quadro,riga);

		cudaMemcpy(slice, dev_slice, (dominio * sizeof(float)+1), cudaMemcpyDeviceToHost);
		
		cudaFree(d_Image_FL);
		cudaFree(d_signal);

	} // fine di teta in funzione di h

	cudaFree(dev_slice);
	cudaFree(d_x_start);
	cudaFree(d_x_end);

	return 0;
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
	//size_t	err=0;
	size_t	count_byte=0;
	//ssize_t		size=0;
	//int		fsize=0;


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
void	FFT(float *row_re,float *row_im,int FFT_lenght) {
	
	float	*data=NULL;
	unsigned long nn;
	int		k;

	nn = FFT_lenght;
	data = (float *) malloc ( sizeof(float)*2*nn+1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		data[2*k+1]=row_re[k];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		data[2*k+2]=row_im[k];
	}

	// FFT by numerical recipes in c
	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
	// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
	// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
	// be an integer power of 2 (this is not checked for!).
	four1(data,nn,1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		row_re[k]=data[2*k+1];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		row_im[k]=data[2*k+2];
	}


	free(data);

	return;
}
//-----------------------------------------------------------------------------
// inverse fft: inverse fast fourier transform
//-----------------------------------------------------------------------------
void	InvFFT(float *row_re,float *row_im,int FFT_lenght) {
	
	float	*data=NULL;
	unsigned long nn;
	int		k;

	nn = FFT_lenght;
	data = (float *) malloc ( sizeof(float)*2*nn+1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		data[2*k+1]=row_re[k];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		data[2*k+2]=row_im[k];
	}

	// FFT by numerical recipes in c
	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
	// data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
	// data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
	// be an integer power of 2 (this is not checked for!).
	four1(data,nn,-1);

	for(k=0; k<FFT_lenght; k++ ) {
		// parte reale 1, 3 ... numeri dispari
		row_re[k]=data[2*k+1];
		// parte immaginaria 2,4...numeri pari a partire dal 2
		row_im[k]=data[2*k+2];
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


