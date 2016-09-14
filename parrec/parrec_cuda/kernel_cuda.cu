#include <cufft.h>
#include "Global.h"



//-----------------------------------------------------------

__global__ void zeropaddingslice(float *dev_slice, int dominio){

	int contatore = blockIdx.x*blockDim.x + threadIdx.x;
	while (contatore < dominio){

		dev_slice[contatore] = 0.;
		contatore += blockDim.x*gridDim.x;
	}
}
//-----------------------------------------------------------

__global__ void creation_xStartEnd(int *d_x_start, int *d_x_end, float r, float l, float rxsize, int localCT_status, int localCT_x0, int localCT_x1, int width, int y_start, int y_end){

	int	xs, xe;
	float yf, radicale;

	int IdThreadInGrid = blockIdx.x*blockDim.x + threadIdx.x;
	int y = y_start + IdThreadInGrid; /*questo ci vuole perchè se no in ogni thread y partirebbe sempre da y_strat mentre ora
										th 0: y=y_start + 0
										th 1: y=y_start + 1
										th 2: y=y_start + 2
										th 3: y=y_start + 3
										     [...]        ...emulando così il for che non parte da zero */

//-----------------------------------------------------------------------------
// funzione per costringere x_start e x_end a stare dentro al cerchio r
// determinazione di x_start e x_end SENZA e CON localtomography
// questa funzione alloca i vettori x_start e x_end e li riempe di valori
//-----------------------------------------------------------------------------
	while(y <= y_end){

		yf = (float)y*rxsize -l;

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

		if(localCT_status == TRUE ) {
			xs = (localCT_x0 > xs) ? localCT_x0 : xs;		// localtomography
			xe = (localCT_x1 < xe) ? localCT_x1 : xe;			// localtomography
		}

		// si memorizzano x_start e x_end calcolati
		d_x_start[y] = xs;
		d_x_end[y] = xe;

		y += blockDim.x*gridDim.x;
	}
}
//-----------------------------------------------------------

__global__ void zp_d_signal(cufftComplex *d_signal, int FFT_lenght){

	int IdThreadInGrid = blockIdx.x*blockDim.x + threadIdx.x;

	int contatore = IdThreadInGrid;
        while (contatore < FFT_lenght){

                d_signal[contatore].x = 0.;
                d_signal[contatore].y = 0.;
                contatore += blockDim.x*gridDim.x;
        }
}
//-----------------------------------------------------------
__global__ void initialize_signal(cufftComplex *d_signal, float *d_Image_FL, int FFT_start, int width, int riga, int FFT_lenght, float *d_detector_weighing){

	int IdThreadInGrid = blockIdx.x*blockDim.x + threadIdx.x;

	int k = FFT_start + IdThreadInGrid;
	while(k<FFT_start+width){

		d_signal[k].x = d_Image_FL[riga+k-FFT_start];
		d_signal[k].x *= d_detector_weighing[k-FFT_start];
		k += blockDim.x*gridDim.x;
	}
}
//-----------------------------------------------------------

__global__ void convolution(cufftComplex *d_signal, float *d_filter, int FFT_lenght){

	int contatore = blockIdx.x*blockDim.x + threadIdx.x;
	while (contatore < FFT_lenght){

		d_signal[contatore].x = d_signal[contatore].x * d_filter[contatore];
		d_signal[contatore].y= d_signal[contatore].y * d_filter[contatore];

		contatore += blockDim.x*gridDim.x;
	}
}
//-----------------------------------------------------------

__global__ void geometry_reconstruction(cufftComplex *d_signal, float *dev_slice, float *d_Image_FL, int *d_x_start, int *d_x_end, float sod_quadro, float pxsize, float r, float rxsize, float xs, float ys, float teta_rad, float xd0, float xd1, float yd0, float yd1, float A, float B, float l, int caso, int nrays, int FFT_lenght, int FFT_start, int h, int localCT_status, int localCT_x0, int localCT_x1, int width, int y_start, int y_end, float sin_teta, float cos_teta, float xs_rot,float ys_rot,float xd0_rot,float xd1_rot,float yd0_rot,float yd1_rot,float xd_min,float xd_max,float yd_min,float yd_max, float xs_rot_quadro, float ys_rot_quadro, int riga){

	int iy = blockIdx.y*blockDim.y + threadIdx.y;
	int ix = blockIdx.x*blockDim.x + threadIdx.x;

	int y = y_start + iy;
	int x = d_x_start[y] + ix;

	int ia, ib;
	float yf, xf, xT, yT, a, b, d, i_float, val, val_pesato_U, dA, dB;
	float U, U_quadro;


	if (y > y_end || x > d_x_end[y])
		return;

	// calcolo di yf
	yf = (float)y*rxsize -l;	// qui andrebbe + rxsize/2

	// calcolo di xf
	xf = (float)x*rxsize -l;		// qui andrebbe + rxsize/2

	// ------ si tracciano le rette ------ //

	// calcoliamo la retta per il punto xf,yf xs_rot,ys_rot
	// calcoliamo l'intersezione fra le due rette xT,yT, nell'ipotesi di detector non verticale
	// la retta del detector è y=Ax+B
	if( caso == 1 && xs_rot == xf ) {
		//ray_verticale=TRUE;   // se è TRUE vuol dire che b=0 e y=xs_rot
		xT = xs_rot;
		yT = xT*A+B;
	}
	else {
		//ray_verticale=FALSE;   // si calcolano normalmente a e b
		a = (ys_rot-yf)/(xs_rot-xf);
		b = yf - xf*a;
		if (caso == 1){
			xT = (b-B)/(A-a);
			yT = A*xT+B;
		}
		else {
			xT = xd0_rot;        // se è TRUE vuol dire che B=0 e y=xd_rot
			yT = xT*a+b;		 // se il detector è verticale non può esserlo anche il raggio in quanto per definizione
		}
	}

	// ora abbiamo calcolato xT,yT vediamo se si trova all'interno del detector
	if( xT >= xd_min && xT <= xd_max && yT >= yd_min && yT <= yd_max) {
		// c'è intersezione calcoliamo i_float
		d = (xT-xd0_rot)*(xT-xd0_rot) + (yT-yd0_rot)*(yT-yd0_rot);
		d = sqrt(d);
		i_float = d/pxsize;
		
		if( i_float >=0. && i_float < nrays) {  // qui ci va nrays x halfscan	
			// interpoliamo
			ia = (int) i_float;
			ib = ia+1;
			dB = i_float - (float)ia;
			dA = 1. - dB;
	
			// calcoliamo il valore interpolato da associare a questo punto
			// qui eventualmente saltare se //						val =0.;
			if( d_Image_FL[riga+ia] > 0. && ib < nrays && d_Image_FL[riga+ib] > 0. ) {
				val = d_signal[ia+FFT_start].x * dA;
				val += d_signal[ib+FFT_start].x * dB;
	
				// peso per la distanza dalla proiezione sulla retta della sorgente ---W3--- EQUISPAZIATI
				U = (sod_quadro+xs_rot_quadro+ys_rot_quadro-2.*xs_rot*xf-2.*ys_rot*yf)/sod_quadro;
				U_quadro = U * U;
				val_pesato_U = val/U_quadro;
				// fattore 1/N della FFT
				val_pesato_U /= (float) FFT_lenght;
	
				// associazione del punto nella slice_rank
				dev_slice[y*width+x] += val_pesato_U;
			}
		} // fine di se i_float valori trovati sono validi
	} // fine di controllo su xT e yT
}
//-----------------------------------------------------------



