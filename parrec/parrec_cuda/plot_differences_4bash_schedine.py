#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
# quello scritto sopra serve per la e accentata
from __future__ import division #questo bisogna metterlo per avere 5/2=2.5 perchè se no lui calcola: 5/2=2
#import matplotlib
import numpy as np # questo mi serve per poter decodificare un raw file
#matplotlib.use("Agg")
from pylab import *
import math
import os
import sys

#----------------------------------------------------------------------------
#            I parte: caricamento e decodificazione dei dati
#----------------------------------------------------------------------------


nameSPR = sys.argv[1]
namefileSDTref = sys.argv[2]
namefileSDTx = sys.argv[3]


# leggo il file .spr della slice che mi interessa per salvarmi le informazioni per la decodifica corretta:
dataSPR = np.loadtxt(nameSPR) #.spr file non sono altro che file di testo txt

width = int(dataSPR[1]) #larghezza: num di pixel in ogni riga
height = int(dataSPR[4]) #altezza: numero delle colonne

reference = np.fromfile(namefileSDTref, dtype='<f', sep="") #carico il file .sdt e lo decodifico con questa funzione
sliceSDT = np.fromfile(namefileSDTx, dtype='<f', sep="")

  #sintassi di np.fromfile((a),(b),(c))
					# (a): apre il file specificato
					# (b): tipo di dati dell'array ritornato. Per i file binari serve per le dimensioni e
					#      l'ordine dei byte degli elementi nel file.
					#  --> Nel nostro caso: <f = little endian 32 bit floating point single precision
					# (c): Separatore tra le voci se il file è un file di testo.
					#  --> Empty separatore ("") significa che il file deve essere trattato come binario.
					#      Spaces (" ") nel separatore zero o più caratteri di spaziatura.
					#      Un separatore che consiste solo di spazi deve corrispondere ad almeno uno spazio bianco. 

#----------------------------------------------------------------------------
#    II parte: calcolo differenza e deviazione standard sulla differenza
#----------------------------------------------------------------------------

xi_reference = reference[:]
xi_sliceSDT = sliceSDT[:]

xid_differenza = abs(xi_reference - xi_sliceSDT)

N = len(xid_differenza)


average = np.average(xid_differenza)
sigma = np.std(xid_differenza, dtype=np.float64)

f = open('schedine_average_diffSlice.txt', 'a')
f.write('%0.10f' % average)
f.write("\n")
f.close()


print '%0.10f' % sigma




#-----------------------------------------------------------------------------


