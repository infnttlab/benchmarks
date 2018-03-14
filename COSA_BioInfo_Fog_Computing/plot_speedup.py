#!/usr/bin/python

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import os.path
import matplotlib.ticker as ticker

def printHelp():
    print '\nUsage: python ' + sys.argv[0] + ' [n-dataset] [data-csv1] [data-csv2] ... [data-csvN]'
    print 'Requirements: - [n-dataset] is INT'
    print '              - csv filename as <hostname>_perf.csv, the delimiter is the comma\n'
    sys.exit()

if len(sys.argv) < 2:
    printHelp()
else:
    if (sys.argv[1] == '--help') or (sys.argv[1] == '-h'):
        printHelp()
    elif len(sys.argv) < 3:
        printHelp()

try:
    dataset_n = int(sys.argv[1])
except ValueError:
    print 'ERROR: Bad first value.'
    printHelp()
  #      pass  # it was a string, not an int.

# ho tanti file quante sono le macchine su cui ho fatto il test,
# e poi passo il numero di dataset usati.
# quindi posso ricavarmi la sucessione dei core usati di volta in volta facendo:
# n-tot-righe-file / n-dataset = n-righe-per-ogni-dataset, dove ogni riga rappresenta 
# l'aumento del numero di core usati di una potenza di due.
# es. n-tot-righe-file = 25; n-dataset = 5; ==> 25/5=5, 5 righe per ogni dataset
# ==> 1a r: 1core, 2a r: 2core, 3a r: 4core, 4a r: 8core, 5a r: 16core usati.

host_n = len(sys.argv) - 2

# for each csv:
# host dbName maxTh nFile n-bp t_deepnV [bp/s] t_kraken [bp/s]

# sara' dentro un ciclo for che cicla sui csv presi come n di argomenti in input..
for hN in range(0, host_n):
    input_index = hN + 2
    csv_file = sys.argv[input_index]
    if not os.path.exists(csv_file):
        print 'ERROR: file '+ csv_file + ' not exist!'
        printHelp()

    print '>>> ' + csv_file

    data_name = []
    deepnano_bps = []
    kraken_bps = []
    row_n = 0

    with open(csv_file) as csvF:
        next(csvF) # skip the first row, the header
        csvReader = csv.reader(csvF, delimiter=',')
        for row in csvReader:
            data_name.append(row[1])
            deepnano_bps.append(float(row[5]))
            kraken_bps.append(float(row[6]))
            row_n = row_n + 1

    # ciclo for in cui ogni iterazione rappresenta un dataset diverso
    # per ognuno di questi devo fare una linea di speedup sul plot

    # Two subplots, sharing X axis
    f, axarr = plt.subplots(2, sharex=True)

    step_n = row_n / dataset_n
    x = [2**n for n in range(0, step_n)]
    for i in range(0, dataset_n):
        deepnano_seq = deepnano_bps[step_n*i]
        kraken_seq = kraken_bps[step_n*i]
        deepnano_speedup = []
        kraken_speedup = []
        for count in range(0, step_n):
            index = (i*step_n)+count
            deepnano_speedup.append(float(deepnano_bps[index]/deepnano_seq))
            kraken_speedup.append(float(kraken_bps[index])/float(kraken_seq))

        name = data_name[step_n*i]
        axarr[0].plot(x, deepnano_speedup, label=name) # top
        axarr[1].plot(x, kraken_speedup) # bottom

    host = csv_file[:-9]
    d_title = 'Deepnano Speedup on ' + host
    k_title = 'Kraken Speedup on ' + host
    axarr[0].set_title(d_title, fontsize='14')
    axarr[1].set_title(k_title, fontsize='14')
    axarr[1].set_xticks(x, minor=False)
    axarr[1].set_xlabel('CPU Cores Used', fontsize='14')
    axarr[0].set_ylabel('Speedup', fontsize='14')
    axarr[1].set_ylabel('Speedup', fontsize='14')
    axarr[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    figTitle = 'speedup_deepnano_kraken_' + host + '.png'
    plt.savefig(figTitle, bbox_inches="tight")
