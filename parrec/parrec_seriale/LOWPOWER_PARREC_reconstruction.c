

//-----------------------------------------------------------------------------
// Questo programma � stato creato da Rosa Brancaccio
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
// PAR_REC_RB.c : definisce il punto di ingresso dell'applicazione console.
//-----------------------------------------------------------------------------
#include <math.h>
#include "Global.h"


//-----------------------------------------------------------------------------
// ********************* MAIN *************************************************
//-----------------------------------------------------------------------------
int main(int argc, char* argv[]) {

	int		r;
	int		start,end;

	// Inizializzazioni dei parametri dell'immagine o della sequenza 
	// solo inizializzazioni e nessuna deallocazione solo qui per costringerlo si pone newdataset TRUE
	DataSet.new_dataset = TRUE;
	// si dealloca la memoria solo se DataSet.new_dataset � FALSE!!!!
	// per essere sicuri di questo si richiama prima il free e poi le inizializzazioni
	Free_AllocatedMemory_DataSet(DATASET1);
	// per evitare deallocazioni di null pointer si inizializzano qui tutti i dati ma non si fanno free
	Initialize_Data_DataSet(DATASET1);
	// si mettono a zero i puntatori di dataset
	Initialize_Pointers_DataSet(DATASET1);
	
	// informazioni per il set di dati da ricostruire (1 sola immagine)
	DataSet.new_dataset = TRUE;
	sprintf(DataSet.Name,"sinos\0");
	DataSet.Seq_start = atoi(argv[1]);
	DataSet.Seq_end = atoi(argv[2]);
	DataSet.current_image = DataSet.Seq_start;
	DataSet.Number_of_Images = (DataSet.Seq_end-DataSet.Seq_start)+1;
	DataSet.File_Type = SDT;

	char name_sinos_spr[999], u_nimage[9];
	//costruisco il nome del sinos da aprire per leggere i dati (uso start, tanto i dati valgono per tutti i sinos della stessa tomografia)
	sprintf(name_sinos_spr,"sinos");
	sprintf(u_nimage,"_%d",DataSet.Seq_start);
	strcat(name_sinos_spr,u_nimage);
	strcat(name_sinos_spr,".spr");

	FILE *fd33;
	int i;
	double read_spr_file[8];
	fd33=fopen(name_sinos_spr, "r");
	if( fd33==NULL ) {
			perror("Errore in apertura del file");
			exit(1);
	}
	for(i=0; i<8; i++)
			fscanf(fd33, "%lf", &read_spr_file[i]);
	fclose(fd33);

	DataSet.Width = (int)read_spr_file[1];
	DataSet.Height = (int)read_spr_file[4];
	DataSet.Data_Type = FLT;
	printf("\nSequence name is <<%s>> and there are <<%d>> images",DataSet.Name,DataSet.Number_of_Images);
	printf("\nImages have %d width and %d height",DataSet.Width,DataSet.Height);

	// aprimao sdt	
	r = Open_SDT_File(DataSet.current_image);
	if( r < 0 )		return;

	// si vede se c'� il file sct e nel caso si carica	
	r = Open_SCT_File();
	if( r < 0 )		return;

	// unica funzione per ricostruire, usa in ogni caso x_start e x_end ma se � stato calcolato il cerchio della slice 
	// questi valori cambieranno, altrimenti varranno rispettivamente 0 e nrays-1
	start = DataSet.Seq_start;
	end = DataSet.Seq_end;
	r = Reconstruction_FAN_BEAM_detector_interpolation_NOPreview(start,end);
	
	// free
	Free_AllocatedMemory_DataSet(DATASET1);
	Initialize_Data_DataSet(DATASET1);
	Initialize_Pointers_DataSet(DATASET1);

   return 0;

}
	
//-----------------------------------------------------------------------------
//	si apre il file sct per leggere alcuni parametri
//  in ogni caso se il file NON esiste o se � stato premuto CANCEL si ritorna -1
//-----------------------------------------------------------------------------
int Open_SCT_File(void) {
	
	char	fname[MAX_PATHNAME_LEN];
	FILE	*file=NULL;
	char	txt[MEDIUM_STRING],txt1[MEDIUM_STRING],buffer[MAX_STRING],buf_errors[MAX_STRING];
	int		r=-1,fileLen,k=0;
	float	val_f=0.;
	int		val=0,err=0;
	char 	*pch;  
	float	m,rxsize,rysize;
	int		rank=0;
	
	// ------------------------------------------------------------------------------------------------------
	// qualsiasi cambiamento nei dati SCT prevede che si definisca la nuova variabile in dataset (global.h), 
	// che venga letta in OpenSCT e salvata in SaveSCT e poi ovviamente usata dove serve!!!!!!!!
	// ------------------------------------------------------------------------------------------------------



	// si vuole cercare automaticamente se c'� gi�
	// un file SCT che vada bene cio� che abbia il nome della sequenza 
	
	sprintf(fname,DataSet.Path);
	
	// si compone il filename
	if(DataSet.StrangeName == TRUE)
		strcat(fname,DataSet.Name);
	else {
		switch(DataSet.TomographicType)  {	//	#define PROJECTION 1   //	#define ATENRAD 2	//	#define SINOS 3	
											//	#define RECOBJ 4	   //	#define GENERIC 5	//	#define NOT_CLASSIFIED 6	
			case(PROJECTION) : { 	strcat(fname,"projection");	break;	}
			case(ATENRAD) : {		strcat(fname,"atenrad");	break;	}
			case(SINOS) : {			strcat(fname,"sinos");		break;	}
			case(RECOBJ) : {		strcat(fname,"recobj");		break;	}
			case(GENERIC) : {		strcat(fname,DataSet.Name);	break;	}
			case(NOT_CLASSIFIED) : {strcat(fname,DataSet.Name);	break;	}
		}
	}

	strcat(fname,".sct");
	
	// Open file
	file = fopen(fname, "r");
	if (!file) 
		return -1;
		

	// se si arriva qui vuol dire che si carica un nuovo file
	fseek(file, 0, SEEK_END); // Get file length
	fileLen=ftell(file);
	fseek(file, 0, SEEK_SET);	
	fread(SCT_Data, fileLen, 1, file); 
	fclose(file);
	
	
	// ora si vanno a leggere i dati che ci servono
	// ogni volta che si trovano dei dati non concordanti si mettono in un buffer
	// e alla fine della lettura viene avvertito l'utente 
	sprintf(buf_errors,"\0");
	
	// arange
	pch = strstr (SCT_Data,"-arange");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+7,99);
		buffer[99] = 0;
		DataSet.arange = (float) atof(buffer);
    }
	// nangles
	pch = strstr (SCT_Data,"-nangles");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+8,99);
		buffer[99] = 0;
		val = atoi(buffer);
  		if( DataSet.TomographicType == PROJECTION || DataSet.TomographicType == ATENRAD ) {
			if (val != DataSet.Number_of_Images) {
				sprintf(txt1,"\n Nangles in the SCT file is %5d and I calculated %5d",val,DataSet.Number_of_Images);  		
  				strcat(buf_errors,txt1);
			}
		} else if ( DataSet.TomographicType == SINOS ) {
  			if (val != DataSet.Height) {
				sprintf(txt1,"\n Nangles in the SCT file is %5d and I calculated %5d",val,DataSet.Height);  		
				strcat(buf_errors,txt1);
  			}
		}
		DataSet.nangles = val;
	}
	// nslices
	pch = strstr (SCT_Data,"-nslices");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+8,99);
		buffer[99] = 0;
		val = atoi(buffer);
		if( DataSet.TomographicType == PROJECTION || DataSet.TomographicType == ATENRAD ) {
			if (val != DataSet.Height ) {
				sprintf(txt1,"\n Nslices in the SCT file is %5d and I calculated %5d",val,DataSet.Height);  		
				strcat(buf_errors,txt1);
			}
		} 
		DataSet.nslices = val;
	}
	
	// nrays
	pch = strstr (SCT_Data,"-nrays");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+6,99);
		buffer[99] = 0;
		val = atoi(buffer);
		if (val != DataSet.Width) {  
			sprintf(txt1,"\n Nrays in the SCT file is %5d and I calculated %5d",val,DataSet.Width);  		
			strcat(buf_errors,txt1);
		} 
		DataSet.nrays = val;
	}
	// sdd
	pch = strstr (SCT_Data,"-sdd");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+4,99);
		buffer[99] = 0;
		DataSet.sdd = (float) atof(buffer);
	}
	// sod
	pch = strstr (SCT_Data,"-sod");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+4,99);
		buffer[99] = 0;
		DataSet.sod = (float) atof(buffer);
	}
	// odd
	pch = strstr (SCT_Data,"-odd");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+4,99);
		buffer[99] = 0;
		DataSet.odd = (float) atof(buffer);
		if( (DataSet.odd + DataSet.sod) != DataSet.sdd ) {
			sprintf(txt1,"\n The sum of odd=%f sod=%f in the SCT file is %f and sdd=%f ",DataSet.odd,DataSet.sod,DataSet.sdd);  		
			strcat(buf_errors,txt1);
		}

	}
	// pxsize
	pch = strstr (SCT_Data,"-pxsize");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+7,99);
		buffer[99] = 0;
		DataSet.pxsize = (float) atof(buffer);
	}
	// pysize
	pch = strstr (SCT_Data,"-pysize");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+7,99);
		buffer[99] = 0;
		DataSet.pysize = (float) atof(buffer);
	}
	// pxcenter
	pch = strstr (SCT_Data,"-pxcenter");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+9,99);
		buffer[99] = 0;
		DataSet.pxcenter = (float) atof(buffer);
	}
	// pzcenter
	pch = strstr (SCT_Data,"-pzcenter");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+9,99);
		buffer[99] = 0;
		DataSet.pzcenter = (float) atof(buffer);
	}
	// rxsize
	pch = strstr (SCT_Data,"-rxsize");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+7,99);
		buffer[99] = 0;
		// prima di memorizzarlo controllare la magnificazione
		// calcolo della magnificazione
		m = DataSet.sdd /DataSet.sod;
		// calcolo degli rxsize e rysize
		rxsize = DataSet.pxsize / m; 
		val_f = (float) atof(buffer);
		if( Absolute(val_f - rxsize) > 0.000000999 ) {  // ci interessano solo le prime sei cifre decimali
			sprintf(txt1,"\n Rxsize in the SCT file is %f and I calculated %f",val_f,rxsize);  		
			strcat(buf_errors,txt1);
		}
		DataSet.rxsize = val_f;
	}
	// rysize
	pch = strstr (SCT_Data,"-rysize");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+7,99);
		buffer[99] = 0;
		// prima di memorizzarlo controllare la magnificazione
		// calcolo della magnificazione
		m = DataSet.sdd /DataSet.sod;
		// calcolo degli rysize e rysize
		rysize = DataSet.pysize / m; 
		val_f = (float) atof(buffer);
		if( Absolute(val_f - rysize) > 0.000000999 ) {  // ci interessano solo le prime sei cifre decimali
			sprintf(txt1,"\n Rysize in the SCT file is %f and I calculated %f",val_f,rysize);  		
			strcat(buf_errors,txt1);
		}
		DataSet.rysize = val_f;
	}
	// trypxcenter
	pch = strstr (SCT_Data,"-trypxcenter");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+12,99);
		buffer[99] = 0;
		val = atoi(buffer);
		DataSet.rec_try_pxcenter = val;
	}
	pch = strstr (SCT_Data,"-trypxcenter_start");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+18,99);
		buffer[99] = 0;
		DataSet.pxcenter_start = (float) atof(buffer);
	}
	pch = strstr (SCT_Data,"-trypxcenter_end");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+16,99);
		buffer[99] = 0;
		DataSet.pxcenter_end = (float) atof(buffer);
	}
	pch = strstr (SCT_Data,"-trypxcenter_step");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+17,99);
		buffer[99] = 0;
		DataSet.pxcenter_step = (float) atof(buffer);
	}
	// local CT			localtomography
	pch = strstr (SCT_Data,"-localCT");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+8,99);
		buffer[99] = 0;
		val = atoi(buffer);
		DataSet.localCT_status = val;
	}
	pch = strstr (SCT_Data,"-localCT_X0");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.localCT_x0 = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-localCT_Y0");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.localCT_y0 = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-localCT_X1");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.localCT_x1 = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-localCT_Y1");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.localCT_y1 = (float) atoi(buffer);
	}
	// vettore dei centri
	pch = strstr (SCT_Data,"-usecentvec");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.use_pxcenter_vector = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-centvec_S0");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.pxcenter_vector_S0 = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-centvec_S1");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+11,99);
		buffer[99] = 0;
		DataSet.pxcenter_vector_S1 = (float) atoi(buffer);
	}
	pch = strstr (SCT_Data,"-centvec_valS0");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+14,99);
		buffer[99] = 0;
		DataSet.pxcenter_vector_valS0 = (float) atof(buffer);
	}
	pch = strstr (SCT_Data,"-centvec_valS1");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+14,99);
		buffer[99] = 0;
		DataSet.pxcenter_vector_valS1 = (float) atof(buffer);
	}
	// halfscan
	pch = strstr (SCT_Data,"-halfscan");
	if( pch != NULL ) {
  	  	strncpy (buffer,pch+9,99);
		buffer[99] = 0;
		DataSet.halfscan = (float) atoi(buffer);
	} 
	

	// fine lettura dati
	
	// output gestione errori
	if( strlen(buf_errors) > 0 ) {
		sprintf(txt,"\n I found some errors reading the SCT files.\nPlease correct the following discrepancies:");
		strcat(txt,buf_errors);
		sprintf(txt1,"\nI am rank %d :",rank);
		printf(txt1);
		printf(txt);
		fflush(0);
	}
	
	// rossella rossella SCT devi leggere anche gli altri dati	
	
	DataSet.SCT_data_loaded=TRUE;

	return 0;	
}

//-----------------------------------------------------------------------------
// si apre un file SDT
//-----------------------------------------------------------------------------
int Open_SDT_File(int n) { 

	FILE 	*file;
	unsigned long fileLen;
	int 	return_value=0;
	size_t	err=0;
	char	PathName[MAX_STRING];
	

	DataSet.current_image = n; // qui va cambiato se dark o bak allora va rimesso il current precedente al di fuori di questa funzione
	
	// si ottiene il nome del file da aprire
	ComposeFileNameFromDataSet(PathName);  

	// ho messo in composefilename DRK per dark	e BAK per I0
	
	
	// Open file
	file = fopen(PathName, "rb");
	if (!file)	{
		if( n == DARK) { 
			DataSet.Dark=FALSE;
			return_value=-1;       
			goto error; 
		} else if ( n == IZERO) { 
			DataSet.Izero=FALSE;
			return_value=-1;       
			goto error; 
		} else {
			printf("\nUnable to open selected file!\n%s",PathName);
			fflush(0);
			return_value=-1;       
			goto error;
		}
	}
	
	
	// Allocate memory  on the basis of the data type
	// DATA TYPE 	1 = unsigned short;	 	3 = floating point   
	
	// si allocano e si caricano le immagini a seconda del tipo di dati
	// e a seconda che si tratti di una dark, Izero oppure una normale projection
	
	if (DataSet.Data_Type == U16 ) {
		
		if( n == DARK) {
		
			// l'allocazione NON si fa solo la prima volta

			// Get file length (solo la prima volta!!)
			fseek(file, 0, SEEK_END);
			fileLen=ftell(file);
			fseek(file, 0, SEEK_SET);
				
				// allocazione
				free(DataSet.Image_U16_Dark);
				DataSet.Image_U16_Dark=NULL;
				DataSet.Image_U16_Dark=(unsigned short *)malloc(fileLen*(sizeof(unsigned short))+1);		
				if (!DataSet.Image_U16_Dark)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for dark image of type unsigned short");
        			fclose(file);
					return_value=-1;       
					DataSet.Dark=FALSE;
					goto error;
				}
			
	
			// Read file contents into buffer
			err = fread(DataSet.Image_U16_Dark, fileLen, 1, file);
			fclose(file);
			if(err == 0) {
				DataSet.Dark=FALSE;
				return_value=-1;       
			} else {
				DataSet.Dark=TRUE;
				return_value=0;
			}
			
		} else if ( n == IZERO ) {
			
			// l'allocazione NON si fa solo la prima volta

			// Get file length (solo la prima volta!!)
			fseek(file, 0, SEEK_END);
			fileLen=ftell(file);
			fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_U16_Izero);
				DataSet.Image_U16_Izero=NULL;
				DataSet.Image_U16_Izero=(unsigned short *)malloc(fileLen*(sizeof(unsigned short))+1);		
				if (!DataSet.Image_U16_Izero)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for Izero image of type unsigned short");
        			fclose(file);
					return_value=-1;       
					DataSet.Izero=FALSE;
					goto error;
				}
			

			// Read file contents into buffer
			err = fread(DataSet.Image_U16_Izero, fileLen, 1, file);
			fclose(file);
			if(err == 0) {
				DataSet.Izero=FALSE;
				return_value=-1;       
			} else {
				DataSet.Izero=TRUE;
				return_value=0;
			}

		} else {			// FILE SDT U16
				
			// l'allocazione si fa solo la prima volta
			if( DataSet.new_dataset == TRUE ) {

				// Get file length (solo la prima volta!!)
				fseek(file, 0, SEEK_END);
				fileLen=ftell(file);
				fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_U16);
				DataSet.Image_U16=NULL;
				DataSet.Image_U16=(unsigned short *)malloc(fileLen*(sizeof(unsigned short))+1);		
				if (!DataSet.Image_U16)	{
   					printf("Warning,Memory error!\nI am not able to allocate memory for 1 image of type unsigned short");
        			fclose(file);
					return_value=-1;       
					goto error;
				}
				DataSet.new_dataset = FALSE;
			} else
				fileLen = DataSet.Width*DataSet.Height*sizeof(unsigned short);
	
			// Read file contents into buffer
			err = fread(DataSet.Image_U16, fileLen, 1, file);
			fclose(file);
			if(err == 0)
				return_value=-1;       
			else
				MaxMin(DATASET1);
		}
		
	} else if (DataSet.Data_Type == FLT ) {
		
		if( n == DARK) {
			
			// l'allocazione NON si fa solo la prima volta
				
				// Get file length (solo la prima volta!!)
				fseek(file, 0, SEEK_END);
				fileLen=ftell(file);
				fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_FL_Dark);
				DataSet.Image_FL_Dark=NULL;
				DataSet.Image_FL_Dark=(float *)malloc(fileLen*(sizeof(float))+1);		
				if (!DataSet.Image_FL_Dark)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for dark image of type float");
        			fclose(file);
					return_value=-1;       
					DataSet.Dark=FALSE;
					goto error;
				}
			
	
			// Read file contents into buffer
			err = fread(DataSet.Image_FL_Dark, fileLen, 1, file);
			fclose(file);
			if(err == 0) {
				DataSet.Dark=FALSE;
				return_value=-1;       
			} else {
				DataSet.Dark=TRUE;
				return_value=0;
			}
			
		} else if ( n == IZERO ) {
			
			// l'allocazione NON si fa solo la prima volta

				// Get file length (solo la prima volta!!)
				fseek(file, 0, SEEK_END);
				fileLen=ftell(file);
				fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_FL_Izero);
				DataSet.Image_FL_Izero=NULL;
				DataSet.Image_FL_Izero=(float *)malloc(fileLen*(sizeof(float))+1);		
				if (!DataSet.Image_FL_Izero)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for Izero image of type float");
        			fclose(file);
					return_value=-1;       
					DataSet.Izero=FALSE;
					goto error;
				}
			
			
			// Read file contents into buffer
			err = fread(DataSet.Image_FL_Izero, fileLen, 1, file);
			fclose(file);
			if(err == 0) {
				DataSet.Izero=FALSE;
				return_value=-1;       
			} else {
				DataSet.Izero=TRUE;
				return_value=0;
			}
			
		} else {			// FILE SDT FLT

			// l'allocazione si fa solo la prima volta
			if( DataSet.new_dataset == TRUE ) {

				// Get file length (solo la prima volta!!)
				fseek(file, 0, SEEK_END);
				fileLen=ftell(file);
				fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_FL); // rossella memory	davanti a tutti i malloc metti il free
				DataSet.Image_FL=NULL;
				DataSet.Image_FL=(float *)malloc(fileLen*(sizeof(float))+1);		
				if (!DataSet.Image_FL)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for 1 image of type float");
        			fclose(file);
					return_value=-1;       
					goto error;
				}
				DataSet.new_dataset = FALSE;
			} else
				fileLen = DataSet.Width*DataSet.Height*sizeof(float);

			// Read file contents into buffer
			err = fread(DataSet.Image_FL, fileLen, 1, file);
			fclose(file);
			if(err == 0)
				return_value=-1;       
			else
				MaxMin(DATASET1);
		}


	} else if (DataSet.Data_Type == U8 ) {

		// rossella rossella al momento non gestiamo dark n� I0 a 8 bit

			// l'allocazione si fa solo la prima volta
			if( DataSet.new_dataset == TRUE ) {

				// Get file length (solo la prima volta!!)
				fseek(file, 0, SEEK_END);
				fileLen=ftell(file);
				fseek(file, 0, SEEK_SET);

				// allocazione
				free(DataSet.Image_U8);
				DataSet.Image_U8=NULL;
				DataSet.Image_U8=(unsigned char *)malloc(fileLen*(sizeof(unsigned char))+1);		
				if (!DataSet.Image_U8)	{
					printf("\nWarning, Memory error!\nI am not able to allocate memory for 1 image of type unsigned char");
        			fclose(file);
					return_value=-1;       
					goto error;
				}
				DataSet.new_dataset = FALSE;
			} else
				fileLen = DataSet.Width*DataSet.Height*sizeof(unsigned char);
	
			// Read file contents into buffer
			err = fread(DataSet.Image_U8, fileLen, 1, file);
			fclose(file);
			if(err == 0)
				return_value=-1;       
			else
				MaxMin(DATASET1);
		//}
		
	} else {
		printf("\nWarning, Data type not recognized");
        fclose(file);
		return_value=-1;       
		goto error;
	}
	
	
error:
	


	return return_value;
}
//-----------------------------------------------------------------------------
// Si apre un file spr e si restituiscono tutte le caratteristiche dell'immagine 
//-----------------------------------------------------------------------------
int OpenSPRFile(char *PathName) {
	
	char		SprPathName[MAX_PATHNAME_LEN]={'\0'};
	char    	*buf=NULL,*buf_rid=NULL;
	int			k=0,return_value=0;
	int  		endline_position=0,temp_position=0,point_position=0,len=0;
	FILE 		*file=NULL;
	unsigned long fileLen=0;
	int			data_dim=0,data_type=0,width=0,height=0;
	float		pxsize=0.0,pysize=0.0,rxsize=0.0,rysize=0.0;
	
	
	
	// per prima cosa costruiamo il nome del file da leggere
	strcpy(SprPathName,PathName);
	point_position = (int) strcspn(SprPathName,".");
	len = (int) strlen(SprPathName);
	if ( point_position == 0 || len == 0) {
		printf("\nWarning! File extension not recognized");
		return_value=-1;
		goto error;
	}
	SprPathName[point_position+1]='\0';
	strcat(SprPathName,"spr");

	
	
	// si legge tutto il file e poi si cercano i vari parametri

	
	// Open spr file
	file = fopen(SprPathName, "r");
	if (!file)
	{
		printf("\nWarning! Unable to open spr selected file");
		return_value=-1;
		goto error;
	}
	
	// Get file length
	fseek(file, 0, SEEK_END);
	fileLen=ftell(file);
	fseek(file, 0, SEEK_SET);

	// Allocate memory
	free(buf);
	buf=NULL;
	buf=(char *)malloc(fileLen+1);		
	if (!buf)
	{
		printf("\nWarning! Memory error!\nI am not able to allocate memory");
        fclose(file);
		return_value=-1;
		goto error;
	}

		
	// Read file contents into buffer
	fread(buf, fileLen, 1, file);
	fclose(file);
	
	// si termina buf
	buf[fileLen]='\0';
	
	// Allocate memory for buf_rid
	free(buf_rid);
	buf_rid=NULL;
	buf_rid=(char *)malloc(fileLen+1);		
	if (!buf_rid)
	{
		printf("\nWarning Memory error!\nI am not able to allocate memory");
        fclose(file);
		return_value=-1;
		goto error;
	}

	// si termina buf
	for(k=0; k<(int)fileLen; k++)
		buf_rid[k]=' ';
	buf_rid[fileLen]='\0';

	// ora si leggono dentro buf i vari valori letti

	// esempio di file spr tipico
	// 2				data_dim
	// 1092				width
	// 0.000000			pxsize
	// 0.000000			pysize
	// 736				height
	// 0.000000			rxsize
	// 0.000000			rysize
	// 1				data_type

	// si legge data_dim il tipo di immagine   (1=linea, 2=bidimensionale, 3=tridimensionale)
	endline_position = (int) strcspn(buf,"\n");
	data_dim = atoi ( &buf[0] ); // atoi ( &buf[endline_position-1] );			   tanto � sempre e solo un numero, il primo!!!
	temp_position = endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	buf_rid[fileLen-temp_position+1]='\0';
	
	// si legge width
	endline_position = (int) strcspn(buf_rid,"\n");
	width = atoi ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge il primo double px-detector
	endline_position = (int) strcspn(buf_rid,"\n");
	pxsize = atof ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge il primo double py-detector
	endline_position = (int) strcspn(buf_rid,"\n");
	pysize = atof ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge height
	endline_position = (int) strcspn(buf_rid,"\n");
	height = atoi ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge il primo double px-object
	endline_position = (int) strcspn(buf_rid,"\n");
	rxsize = atof ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge il primo double py-object
	endline_position = (int) strcspn(buf_rid,"\n");
	rysize = atof ( buf_rid );
	temp_position += endline_position;
	for (k=temp_position+1; k<=(int)fileLen; k++)
		buf_rid[k-temp_position]=buf[k];
	
	// si legge data_type		 // 1 = unsigned short 3 = floating point
	endline_position = (int) strcspn(buf_rid,"\n");
	data_type = atoi ( buf_rid );
	
	DataSet.Data_Dim = data_dim;
	DataSet.Data_Type = data_type;
	DataSet.Width = width;
	DataSet.Height = height;
//	DataSet.pxsize= pxsize;
//	DataSet.pysize= pysize;
//	DataSet.rxsize = rxsize;
//	DataSet.rysize = rysize;
// per il momento si leggono i psize e gli rsize solo da sct
// bisognerebbe valutare se non ci sono nell'sct di farli leggere qui...


error:
	
	free(buf);
	free(buf_rid);
	buf=NULL;
	buf_rid=NULL;
	
	
	return return_value;
}
//-----------------------------------------------------------------------------
// si salva l'immagine in memoria SDT
// se in path c'� tomostep allora si salva automaticamente il file in memoria e SOVRASCRIVE
// altrimenti si usa il path specificato e si chiede prima di sovrascrivere
// ritorna -1 se anche con 99 tentativi non si � riusciti a salvare
// ritrna -2 se il tipo di file non � SDT salvabile
//-----------------------------------------------------------------------------
int Save_SDT_File(char *path) {

	FILE 	*f;
	size_t	err=0;
	char	PathName[MAX_STRING];
	char	txt[MAX_STRING],txt1[SHORT_STRING];
	int		type=0;
	size_t	count_byte=0;
	ssize_t		size=0;
	int		fsize=0;

	// si ottiene il nome del file da salvare
	// a seconda che si tratti di salvataggio automatico o selettivo
	err = strcmp(path,"atenrad");
	if(err == 0 ) {
		txt[0]='\0';
		txt1[0]='\0';
		sprintf(txt,DataSet.Path);
		strcat(txt,"atenrad");
		sprintf(txt1,"_%d",DataSet.current_image);
		strcat(txt,txt1);
		strcat(txt,".sdt");
		CopyStringToStringAndTerminateIt(PathName,txt);
		type = FLT;
	} else {
		err = strcmp(path,"sinos");
		if(err == 0 ) {
			txt[0]='\0';
			txt1[0]='\0';
			sprintf(txt,DataSet.Path);
			strcat(txt,"sinos");
			sprintf(txt1,"_%d",DataSet.current_image);
			strcat(txt,txt1);
			strcat(txt,".sdt");
			CopyStringToStringAndTerminateIt(PathName,txt);
			type = FLT;
		} else {
			err = strcmp(path,"recobj");
			if(err == 0 ) {
				txt[0]='\0';
				txt1[0]='\0';
				sprintf(txt,DataSet.Path);
				strcat(txt,"recobj");
				sprintf(txt1,"_%d",DataSet.current_image);
				strcat(txt,txt1);
				strcat(txt,".sdt");
				CopyStringToStringAndTerminateIt(PathName,txt);
				type = FLT;
			} else {   // caso generico, si usa il pathname scelto dall'utente
				CopyStringToStringAndTerminateIt(PathName,path);
				type = DataSet.Data_Type;
			}
		}
	} 
		
	// si salva il file SDT nel percorso stabilito

	// si creano file binari con dentro i dati
	if ( type == FLT ) {
		// si scrive il file
		f = fopen(PathName,"wb");
		count_byte = sizeof(float)*DataSet.Width*DataSet.Height;
		fwrite(DataSet.Image_FL,count_byte,1,f); 
		fclose(f);
	} else if ( type == U16 ) { 
		// si scrive il file
		f = fopen(PathName,"wb");
		count_byte = sizeof(unsigned short)*DataSet.Width*DataSet.Height;
		fwrite(DataSet.Image_U16,count_byte,1,f); 
		fclose(f);
	} else if ( type == U8) { 
		// si scrive il file
		f = fopen(PathName,"wb");
		count_byte = sizeof(unsigned char)*DataSet.Width*DataSet.Height;
		fwrite(DataSet.Image_U8,count_byte,1,f); 
		fclose(f);
	}
	
 return 0;
	
}
//-----------------------------------------------------------------------------
// si salva il file SPR
// se in path c'� tomostep allora si salva automaticamente la sequenza e SOVRASCRIVE
// altrimenti si usa il path specificato e si chiede prima di sovrascrivere
//-----------------------------------------------------------------------------
int Save_SPR_File(char *path) {

	FILE 	*f;
	size_t	err=0;
	char	txt[MAX_STRING],txt1[SHORT_STRING];
	int		type=0;
	char	SprPathName[MAX_PATHNAME_LEN]={'\0'};
	
	
	// si ottiene il nome del file da salvare
	// a seconda che si tratti di salvataggio automatico o selettivo
	err = strcmp(path,"atenrad");
	if(err == 0 ) {
		txt[0]='\0';
		txt1[0]='\0';
		sprintf(txt,DataSet.Path);
		strcat(txt,"atenrad");
		sprintf(txt1,"_%d",DataSet.current_image);
		strcat(txt,txt1);
		strcat(txt,".spr");
		CopyStringToStringAndTerminateIt(SprPathName,txt);
		type = FLT;
	} else {
		err = strcmp(path,"sinos");
		if(err == 0 ) {
			txt[0]='\0';
			txt1[0]='\0';
			sprintf(txt,DataSet.Path);
			strcat(txt,"sinos");
			sprintf(txt1,"_%d",DataSet.current_image);
			strcat(txt,txt1);
			strcat(txt,".spr");
			CopyStringToStringAndTerminateIt(SprPathName,txt);
			type = FLT;
		} else {
			err = strcmp(path,"recobj");
			if(err == 0 ) {
				txt[0]='\0';
				txt1[0]='\0';
				sprintf(txt,DataSet.Path);
				strcat(txt,"recobj");
				sprintf(txt1,"_%d",DataSet.current_image);
				strcat(txt,txt1);
				strcat(txt,".spr");
				CopyStringToStringAndTerminateIt(SprPathName,txt);
				type = FLT;
			} else {   // caso generico, si usa il SprPathName scelto dall'utente
				CopyStringToStringAndTerminateIt(SprPathName,path);
				type = DataSet.Data_Type;
			}
		}
	} 

	// Save spr file
	f = fopen(SprPathName, "w");
	if (!f)
	{
		printf("\n Warning! Unable to save spr selected file");
		return -1;
	}

	// si fa l'spr
	fprintf(f,"2\n");
	fprintf(f,"%d\n",DataSet.Width);
	fprintf(f,"0.000000\n");
	fprintf(f,"0.000000\n");
//	if( strcmp(path,"sinos") == 0 ) // se stiamo scrivendo i sinogrammi height � il numero di atenrad oppure il range angolare rossella rossella
	if( DataSet.CT_Step_todo == SINOS )
		fprintf(f,"%d\n",DataSet.Seq_end-DataSet.Seq_start+1);
	else
		fprintf(f,"%d\n",DataSet.Height);
	fprintf(f,"0.000000\n");
	fprintf(f,"0.000000\n");
	fprintf(f,"%d\n",type);
	fclose(f);	

	return 0;
}
//-----------------------------------------------------------------------------
// free delle allocazioni di un dataset generico passato come valore int
// si dealloca la memoria solo se DataSet.new_dataset � FALSE!!!!
// per essere sicuri di questo si richiama prima il free e poi le inizializzazioni
//-----------------------------------------------------------------------------
void Free_AllocatedMemory_DataSet(int dataset) {

	switch(dataset) {
		
		case(DATASET1) : {
			
			// si dealloca la memoria solo se DataSet.new_dataset � FALSE!!!!
			if( DataSet.new_dataset == FALSE ) {
		
				switch( DataSet.Data_Type) {
					case(U8) : {
						free(DataSet.Image_U8);
						DataSet.Image_U8 = NULL;
						break;
					}
					case(U16) : {
						free(DataSet.Image_U16);
						DataSet.Image_U16 = NULL;
						if(DataSet.Izero == TRUE) {
							free(DataSet.Image_U16_Izero);
							DataSet.Image_U16_Izero = NULL;
						}
						if(DataSet.Dark == TRUE) {
							free(DataSet.Image_U16_Dark);
							DataSet.Image_U16_Dark = NULL;
						}
						break;
					}
					case(FLT) : {
						free(DataSet.Image_FL);
						DataSet.Image_FL = NULL;
						if(DataSet.Izero == TRUE) {
							free(DataSet.Image_FL_Izero);
							DataSet.Image_FL_Izero = NULL;
						}
						if(DataSet.Dark == TRUE) {
							free(DataSet.Image_FL_Dark);
							DataSet.Image_FL_Dark = NULL;
						}
						break;
					}
				} // fine dello switch per sapere cosa si deve eventualmente deallocare
			}
			break;
		}

		case(DATASET2) : {

			// si dealloca la memoria solo se DataSet.new_dataset � FALSE!!!!
			if( DataSet2.new_dataset == FALSE ) {
				
				switch( DataSet2.Data_Type) {
					case(U8) : {
						free(DataSet2.Image_U8);
						DataSet2.Image_U8 = NULL;
						break;
					}
					case(U16) : {
						free(DataSet2.Image_U16);
						DataSet2.Image_U16 = NULL;
						break;
					}
					case(FLT) : {
						free(DataSet2.Image_FL);
						DataSet2.Image_FL = NULL;
						break;
					}
				} // fine dello switch per sapere cosa si deve eventualmente deallocare
			}
			break;
		}
			
	
	}
	
	return;
}
//-----------------------------------------------------------------------------
// Inizializzazioni di un dataset generico passato come valore int
// si inizializza tutto tranne i puntatori da allocare
//-----------------------------------------------------------------------------
extern void Initialize_Data_DataSet(int dataset) {

	switch(dataset) {

		case(DATASET1) : {

			DataSet.Number_of_Images	= NO_IMAGE;
			DataSet.Single_Image_inSequence = FALSE;
			DataSet.Seq_start 			= 0;
			DataSet.Seq_end 			= 0;
			DataSet.File_Type 			= 0;
			DataSet.Data_Dim 			= IMG_2D;
			DataSet.Data_Type 			= U16;
			DataSet.Width 				= NO_WIDTH;
			DataSet.Height 				= NO_HEIGHT;
			DataSet.pxsize 		= ZERO_F;
			DataSet.pysize 		= ZERO_F;
			DataSet.rxsize 			= 0.0;
			DataSet.rysize 			= 0.0;
			DataSet.current_image 		= NO_SEQ;
			DataSet.min 				= 0 ;
			DataSet.max 				= MAX_16BIT;
			DataSet.min_f 				= -FLT_MAX;
			DataSet.max_f 				= FLT_MAX;
			DataSet.TomographicType 	= NOT_CLASSIFIED;
			DataSet.StrangeName 		= FALSE;
			sprintf(DataSet.Path,"\0");
			sprintf(DataSet.Name,"\0");
			DataSet.Dark 				= FALSE;
			DataSet.Izero 				= FALSE;
			DataSet.CT_Step_todo 		= NOT_CLASSIFIED;
			DataSet.SCT_data_loaded		=FALSE;	// ci sono o no i dati sct?
			DataSet.ringo_percent		= 0.;
			DataSet.ringo_percent2		= 0.;
			DataSet.outlier_percent		= 0.3;
			DataSet.Utility_Step_todo 	= NO_STEP;
			DataSet.Seq_min 			= NO_VALUE;
			DataSet.Seq_max 			= NO_VALUE;
			DataSet.Seq_min_f 			= NO_VALUE_F;
			DataSet.Seq_max_f 			= NO_VALUE_F;
			DataSet.new_dataset 		= TRUE; // � un nuovo data set e quindi si rialloca
			DataSet.metalartifact_percent 	= 0.;
			DataSet.metalartifact_pixel 	= 0.;
			DataSet.metalartifact_down		= 0.;
			DataSet.metalartifact_correction=-1;
			DataSet.FFT_do_invFFT		= TRUE;
			DataSet.FFT_do_FeldkampWeight = TRUE;
			DataSet.arange 			= ZERO_F;
			DataSet.nslices			= ZERO;
			DataSet.nangles			= ZERO;
			DataSet.nrays			= ZERO;
			DataSet.sdd				= ZERO_F;
			DataSet.odd				= ZERO_F;
			DataSet.sod				= ZERO_F;
			DataSet.pxcenter		= ZERO_F ;			
			DataSet.pzcenter		= ZERO_F ;		
			DataSet.xsource			= ZERO_F ;		
			DataSet.zsource			= ZERO_F ;	
			sprintf(DataSet.recobj_name,"recobj_ok");	// nome del file della recobj: se � una sequenza senza "_"
			DataSet.rec_try_pxcenter = FALSE;
			DataSet.pxcenter_start   = 0.;
			DataSet.pxcenter_end     = 0.;
			DataSet.pxcenter_step    = 0.;		
			DataSet.localCT_status=FALSE;		// se � TRUE si vuole ricostruire in localtomography
			DataSet.localCT_x0=0;			// parametri localtomography
			DataSet.localCT_y0=0;			// parametri localtomography
			DataSet.localCT_x1=0;			// parametri localtomography
			DataSet.localCT_y1=0;			// parametri localtomography
			DataSet.use_pxcenter_vector=FALSE;			// se si usa o meno il vettore dei centri
			DataSet.pxcenter_vector_valS0=0.;	// valore del centro nel sinogramma S0
			DataSet.pxcenter_vector_valS1=0.;	// valore del centro nel sinogramma S1
			DataSet.pxcenter_vector_S0=0;		// numero del sinogramma S0
			DataSet.pxcenter_vector_S1=0;		// numero del sinogramma S1
			DataSet.ringonew_T_singoli=1.;		// fattore moltiplicativo per determinare la soglia per ring singoli nell'algoritmo RingoNew
			DataSet.ringonew_T_doppi=1.;		// fattore moltiplicativo per determinare la soglia per ring doppi nell'algoritmo RingoNew
			DataSet.ringonew_T_tripli=1.;		// fattore moltiplicativo per determinare la soglia per ring tripli nell'algoritmo RingoNew
			DataSet.halfscan=FALSE;			// si ricostruisce in half scan
			DataSet.hs_left=0;			// HALFSCAN: allargamento a sinistra
			DataSet.hs_right=0;			// HALFSCAN: allargamento a destra
			DataSet.hs_corrtype=HS_LINEAR_CORRECTION;		// HALFSCAN: tipo di correzione
			DataSet.REC_geometry=FAN_BEAM;		// tipo di geometria: FAN_BEAM oppure CONE_BEAM
			break;
		}

		case(DATASET2) : {

			DataSet2.Number_of_Images	= NO_IMAGE;
			DataSet2.Single_Image_inSequence = FALSE;
			DataSet2.Seq_start 			= 0;
			DataSet2.Seq_end 			= 0;
			DataSet2.File_Type 			= 0;
			DataSet2.Data_Dim 			= IMG_2D;
			DataSet2.Data_Type 			= U16;
			DataSet2.Width 				= NO_WIDTH;
			DataSet2.Height 			= NO_HEIGHT;
			DataSet2.pxsize 		= ZERO_F;
			DataSet2.pysize 		= ZERO_F;
			DataSet2.rxsize 			= 0.0;
			DataSet2.rysize 			= 0.0;
			DataSet2.current_image 		= NO_SEQ;
			DataSet2.min 				= 0 ;
			DataSet2.max 				= MAX_16BIT;
			DataSet2.min_f 				= -FLT_MAX;
			DataSet2.max_f 				= FLT_MAX;
			DataSet2.TomographicType 	= NOT_CLASSIFIED;
			DataSet2.StrangeName 		= FALSE;
			sprintf(DataSet2.Path,"\0");
			sprintf(DataSet2.Name,"\0");
			DataSet2.Dark 				= FALSE;
			DataSet2.Izero 				= FALSE;
			DataSet2.CT_Step_todo 		= NOT_CLASSIFIED;
			DataSet2.SCT_data_loaded	=FALSE;	// ci sono o no i dati sct?
			DataSet2.ringo_percent		= 0.;
			DataSet2.ringo_percent2		= 0.;
			DataSet2.outlier_percent	= 0.3;
			DataSet2.Utility_Step_todo 	= NO_STEP;
			DataSet2.Seq_min 			= NO_VALUE;
			DataSet2.Seq_max 			= NO_VALUE;
			DataSet2.Seq_min_f 			= NO_VALUE_F;
			DataSet2.Seq_max_f 			= NO_VALUE_F;
			DataSet2.new_dataset 		= TRUE; // � un nuovo data set e quindi si rialloca
			DataSet2.metalartifact_percent 	= 0.;
			DataSet2.metalartifact_pixel 	= 0.;
			DataSet2.metalartifact_down		=0.;
			DataSet2.metalartifact_correction=-1;
			DataSet2.FFT_do_invFFT		= TRUE;
			DataSet2.FFT_do_FeldkampWeight = TRUE;
			DataSet2.arange 			= ZERO_F;
			DataSet2.nslices			= ZERO;
			DataSet2.nangles			= ZERO;
			DataSet2.nrays				= ZERO;
			DataSet2.sdd				= ZERO_F;
			DataSet2.odd				= ZERO_F;
			DataSet2.sod				= ZERO_F;
			DataSet2.pxcenter			= ZERO_F ;			
			DataSet2.pzcenter			= ZERO_F ;			
			DataSet2.xsource			= ZERO_F ;		
			DataSet2.zsource			= ZERO_F ;		
			DataSet2.rec_slice_circle_present	= FALSE;	// se � TRUE sono presenti in memoria x_start,x_end per la ricostruzione dentro al cerchio della slice
			DataSet2.rec_slice_circle_use		= TRUE;	// se � TRUE vanno usati x_start,x_end per la ricostruzione dentro al cerchio della slice
			sprintf(DataSet2.recobj_name,"recobj_ok");	// nome del file della recobj: se � una sequenza senza "_"
			DataSet2.rec_try_pxcenter = FALSE;
			DataSet2.pxcenter_start   = 0.;
			DataSet2.pxcenter_end     = 0.;
			DataSet2.pxcenter_step    = 0.;		
			DataSet2.localCT_status=FALSE;		// se � TRUE si vuole ricostruire in localtomography
			DataSet2.localCT_x0=0;			// parametri localtomography
			DataSet2.localCT_y0=0;			// parametri localtomography
			DataSet2.localCT_x1=0;			// parametri localtomography
			DataSet2.localCT_y1=0;			// parametri localtomography
			DataSet2.use_pxcenter_vector=FALSE;			// se si usa o meno il vettore dei centri
			DataSet2.pxcenter_vector_valS0=NO_VALUE_F;	// valore del centro nel sinogramma S0
			DataSet2.pxcenter_vector_valS1=NO_VALUE_F;	// valore del centro nel sinogramma S1
			DataSet2.pxcenter_vector_S0=NO_VALUE;		// numero del sinogramma S0
			DataSet2.pxcenter_vector_S1=NO_VALUE;		// numero del sinogramma S1
			DataSet2.ringonew_T_singoli=1.;		// fattore moltiplicativo per determinare la soglia per ring singoli nell'algoritmo RingoNew
			DataSet2.ringonew_T_doppi=1.;		// fattore moltiplicativo per determinare la soglia per ring doppi nell'algoritmo RingoNew
			DataSet2.ringonew_T_tripli=1.;		// fattore moltiplicativo per determinare la soglia per ring tripli nell'algoritmo RingoNew
			DataSet2.halfscan=FALSE;			// si ricostruisce in half scan
			DataSet2.hs_left=0;			// HALFSCAN: allargamento a sinistra
			DataSet2.hs_right=0;			// HALFSCAN: allargamento a destra
			DataSet2.hs_corrtype=HS_LINEAR_CORRECTION;		// HALFSCAN: tipo di correzione
			DataSet2.REC_geometry=FAN_BEAM;		// tipo di geometria: FAN_BEAM oppure CONE_BEAM
			break;
		}

	}
	return;
	
}
//-----------------------------------------------------------------------------
// Inizializzazioni di un dataset generico passato come valore int
// si inizializzano SOLO i puntatori da allocare
//-----------------------------------------------------------------------------
void Initialize_Pointers_DataSet(int dataset) {


	switch(dataset) {

		case(DATASET1) : {
        
			DataSet.Image_U16=NULL;				// puntatore all'immagine U16
			DataSet.Image_FL =NULL; 			// puntatore all'immagine FLOAT
			DataSet.Image_U8=NULL;				// puntatore all'immagine U8
			DataSet.Image_U16_stretch=NULL;		// puntatore all'immagine U16 stretchata
			DataSet.Image_FL_stretch =NULL; 	// puntatore all'immagine FLOAT stretchata
			DataSet.Image_U8_stretch=NULL;		// puntatore all'immagine U8
			DataSet.Image_U16_Dark=NULL;		// puntatore all'immagine U16 DARK
			DataSet.Image_U16_Izero=NULL;		// puntatore all'immagine U16 I0
			DataSet.Image_FL_Dark =NULL;		// puntatore all'immagine FLOAT DARK
			DataSet.Image_FL_Izero =NULL;		// puntatore all'immagine FLOAT I0
			break;
		}

		case(DATASET2) : {
			
			DataSet2.Image_U16=NULL;				// puntatore all'immagine U16
    		DataSet2.Image_FL =NULL; 			// puntatore all'immagine FLOAT
			DataSet2.Image_U8=NULL;				// puntatore all'immagine U8
			DataSet2.Image_U16_stretch=NULL;		// puntatore all'immagine U16 stretchata
			DataSet2.Image_FL_stretch =NULL; 	// puntatore all'immagine FLOAT stretchata
			DataSet2.Image_U8_stretch=NULL;		// puntatore all'immagine U8
			DataSet2.Image_U16_Dark=NULL;		// puntatore all'immagine U16 DARK
			DataSet2.Image_U16_Izero=NULL;		// puntatore all'immagine U16 I0
			DataSet2.Image_FL_Dark =NULL;		// puntatore all'immagine FLOAT DARK
			DataSet2.Image_FL_Izero =NULL;		// puntatore all'immagine FLOAT I0
			break;
		}
	}
	
	return;
		
}
//-----------------------------------------------------------------------------
// si fa partire la ricostruzione con PREVIEW con formule EQUISPAZIATI
// qui si tracciano le rette ed � molto pi� veloce con subFUNZIONI
// INTERPOLAZIONE SUL DETECTOR ------- NUOVO ------------- SENZA PREVIEW DA SCRIVERE!!!!
// ------- SENZA PREVIEW ----------- CON CONSTRAIN CIRCLE -------------
// // **** 2013 07 19 ANCHE LOCAL CT **** ANCHE HALFSCAN
//-----------------------------------------------------------------------------
int	Reconstruction_FAN_BEAM_detector_interpolation_NOPreview(int start,int end) {

	int		i=0;
	int		width,new_width,new_height;
	float	wmezzif;
	float	pxsize,rxsize,sod,sdd,odd;
	float	r=0.,fanangle=0.,l=0.;
	// variabili per stare dentro al cerchio e localtomography
	float	r_pixel=0.;
	int		y_start=0,y_end=0;
	int		x_start=0,x_end=0;

	float	xs,ys;
	float	xd0,yd0,xd1,yd1;

	int		nrays;
	float	pxcenter,shift,shift_pixel;

	int		k;

	int		FFT_lenght,FFT_start;
//	double	n; // n � definito in modo che sia zero in pxcenter
	double	W2; // peso ---W2--- *0.5

	float	W4=0.;  // peso --- W4 ---
	
	float 	min_f=0.,max_f=0.,val_f=0.;
	int		image=0;
	char	path[MAX_STRING],path_SDT[MAX_STRING],path_SPR[MAX_STRING],txt1[SHORT2_STRING];
	int		progress_bar=0;
	char	txt[MEDIUM_STRING];

	// par
	int		rank=0,size=0;
	double	mytime=0.;
	double	totmytime=0.;

	printf("\nCalculating and saving recostructed image...");
	printf("\nREC: fan beam detector interpolation FBP");  // poi metteremo i valori usati
	printf("\nREC: start=%d end=%d (tot=%d)",start,end,end-start+1);
	printf("\nREC: pxcenter=%f arange=%f nangles=%d",DataSet.pxcenter,DataSet.arange,DataSet.nangles);
	if(DataSet.localCT_status == TRUE) {
		printf(txt,"\nREC: Local Tomography data X0=%d Y0=%d X1=%d Y1=%d",DataSet.localCT_x0,DataSet.localCT_y0,DataSet.localCT_x1,DataSet.localCT_y1);
	}
	
	// la slice � quadrata di lato width
	width = DataSet.Width;
	new_width = width;
	new_height = width;
	
	// allochiamo slice 
	free(slice);			slice = NULL;			slice = (float *) malloc ( sizeof(float)*width*width+1);

	if( slice == NULL ) {
		printf("\n Rank 0 : I cannot allocate memory for reconstruction");
		return -1;
	}

	nrays  = DataSet.nrays; // vale width a meno che non si faccia mezzo giro HALFSCAN

	// per allocare row e filter � necessario calcolare la lunghezza degli array nel dominio di fourier
	// � necessario avere almeno il doppio meno uno della dimensione dei punti da elaborare
	// inoltre se si usa un array con lunghezza potenza di 2 la FFT lavora meglio per cui:
	FFT_lenght = CalcTheMinimumNumberPowerOfTwoAfter(2*width-1);

	printf("\nCalculated FFT_lenght is %d",FFT_lenght);
	
	free(row_re);				row_re = NULL;				row_re = (float *) malloc ( sizeof(float)*FFT_lenght+1);
	free(row_im);				row_im = NULL;				row_im = (float *) malloc ( sizeof(float)*FFT_lenght+1);
	free(filter);				filter = NULL;				filter = (float *) malloc ( sizeof(float)*FFT_lenght+1);
	free(detector_weighing);	detector_weighing = NULL;	detector_weighing = (float *) malloc ( sizeof(float)*width+1);
	// --------------------------------------------------

	if( row_re == NULL || row_im == NULL || filter == NULL || detector_weighing == NULL ) {
		printf("\n Rank 0 : I cannot allocate memory for reconstruction");
		return -1;
	}

	// mettiamo i valori che usiamo in delle variabili locali
	odd = DataSet.odd;
	sod = DataSet.sod;
	sdd = DataSet.sdd;
	pxsize= DataSet.pxsize;					
	rxsize= DataSet.rxsize;
	pxcenter = DataSet.pxcenter;
	// width mezzif � width trasformato in float e diviso 2
	wmezzif = ( (float)width)/2.;
	// se pxcenter non � al centro del detector allora definiamo lo shift:
	shift_pixel = (pxcenter-wmezzif);  // � lo shift in pixel dal centro di width
	shift = (pxcenter-wmezzif)*pxsize; // � lo shift in mm

	// qui non serve pi� il calcolo di fi e r!!!!!!!!!!!!
	// calcoliamo fanangle e r con - shift 
	Calculate_r_fi(FAN_BEAM,&r,&fanangle);

	printf("\nCalculated geometry: pxsize=%f rxsize=%f pxcenter=%f shift=%f",pxsize,rxsize,pxcenter,shift);
	printf("\nCalculated geometry: r=%f fi=%f",r,fanangle);
	printf("\nCalculated geometry: sdd=%f sod=%f odd=%f",sdd,sod,odd);


	r_pixel = r/rxsize;
	//printf("\nNUOVI  VALORI: fi=%f   gradi=%f   r=%f    r_pixel=%f",fi,fi*360./(2.*PIGRECO),r,r_pixel);
	
	// calcoliamo le coordinate iniziali della sorgente che poi ruotando varieranno
	xs = 0; // -shift;  // DARIMETTERE 0; 
	ys = sod;   

	// coordinate x e y iniziale del detector che poi ruotando varier� 
	// calcoliamo d0 e d1 per tracciare le rette e interpolare sulla slice
	xd0 = -wmezzif*pxsize-shift;
	yd0 = - odd;
	xd1 = +wmezzif*pxsize-shift;     // attenzione qui si presuppone che nrays==width
	yd1 = - odd;					   // se si fa halfscan qui ci va nrays/2

	// coordinate del centro dell'oggetto: con il nuovo metodo mi serve solo per la visualizzazione perch� poi non lo user� in quanto nulle

	// dati per FFT: abbiamo gi� calcolato la lunghezza necessaria per la FFT ora calcoliamo il punto dove inserire
	// la proiezione che si trover� al centro dell'array nello spazio delle frequenze
	FFT_start = FFT_lenght/2;

	// pesatura per equispaziati SOLO FAN BEAM ---W1--- N.B. non dipende da h quindi facciamo i calcoli fuori dal ciclo!!!
	CalculateWeighing_OnTheBaseOfCTGeometry(EQUISPACED_RAYS,FFT_start,nrays,sod,pxsize,pxcenter);

	// calcolo del filtro da applicare nello spazio di Fourier a partire da RAMP e con BUTTER (W2)
	W2=1.; // peso ---W2--- dovrebbe essere: W2=0.5/(double)pxsize	//attenzione che questo filtro � troppo pesante! va meglio con W2=1.
	CalculateFourierFilter(FFT_lenght, FFT_start, W2);

	// calcolo del peso per equispaziati U --- W3 --- N.B. il calcolo di U � dentro a BackProjectFilteredDataOnSlice
	//CalcolateWeighing_EquispacedRays(sod,rxsize);

	// calcolo del peso nel caso di halfscan andrebbe messo qui e moltiplicato per W1 cos� non devi n� memorizzarlo n� ricalcolarlo ogni volta
	printf("\nCalculated filter and weights");


// si cicla su image
//******************************************************************************
for(image=start; image<=end; image++) {	
//******************************************************************************
	
	// se si stanno elaborando pi� file (o anche uno solo ma non c'� preview) bisogna caricare il file da ricostruire
	Open_SDT_File(image);

	//---------------------------------------------------------
	// ciclo su teta in funzione di h
	//---------------------------------------------------------
	printf("\nBackprojection started, please wait...\n");
	BackProjectFilteredDataOnSlice(FFT_lenght,FFT_start,xs,ys,xd0,yd0,xd1,yd1);
	printf("\nBackprojection finished");

	// calcolo di W4 
	W4 = pxsize;	//rxsize*(float)DataSet.nangles/PIGRECO;	// secondo DAN

	// si ripulisce l'immagine, si tolgono valori negativi aumenta il contrasto si pesa per ---W4---
	for(k=0; k<width*width; k++) {
		val_f = slice[k];
		if ( val_f < 0. ) 
			slice[k] = 0.;
		else
			slice[k] = val_f*W4;   // --- W4 ---
	}
	
	// anche NON si tratta di una preview vogliamo salvare il file
	// si ricava il PATH nome dell'immagine corretta da salvare sia SDT che SPR
	sprintf(path,DataSet.Path);
	strcat(path,DataSet.recobj_name);
	sprintf(txt1,"_%d",image);
	if( DataSet.rec_try_pxcenter == TRUE)
		sprintf(txt1,"_%d",n_pxcenter);
	strcat(path,txt1);
			
	// si salva l'SDT
	sprintf(path_SDT,path);
	strcat(path_SDT,".sdt");
	Save_SDT_GenericFile_FLT(slice,path_SDT,new_width,new_height);
		
	// si salva l'SPR
	sprintf(path_SPR,path);
	strcat(path_SPR,".spr");
	Save_SPR_GenericFile(path_SPR,new_width,new_height,IMG_2D,FLT,SAVE_RSIZE);

	printf("\nSlice saved");


//******************************************************************************
} // fine ciclo su image
//******************************************************************************
	

	printf("\nFinished reconstructing the sequence");

	// si libera la memoria allocata
	if(!slice)
		free(slice);
	slice = NULL;
	if(!row_re)
		free(row_re);
	row_re = NULL;
	if(!row_im)
		free(row_im);
	row_im = NULL;
	if(!filter)
		free(filter);
	filter = NULL;
	if(!detector_weighing)
		free(detector_weighing);
	detector_weighing = NULL;
	// --------------------------------------------------
	// HALFSCAN
	// --------------------------------------------------
	if(!row_modified)
		free(row_modified);				
	row_modified = NULL;	
	// --------------------------------------------------

	return 0;

}

