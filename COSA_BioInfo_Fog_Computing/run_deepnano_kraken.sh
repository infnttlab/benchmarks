#!/bin/bash

if [[ ( -z "$2" ) ]]
then
        echo -e "\\nUsage: $0  [n-proc]  [db-name] [n-file]\\n"
        exit 1
fi

if [[ ( $1 == "--help" ) || ( $1 == "-h" ) ]]
then
        echo -e "\\nUsage: $0  [n-proc]  [db-name] [n-file]\\n"
        exit 0
fi

host=$(hostname -s)

echo " "
echo "This script takes the time-performances of deepnano and kraken using ${2} dataset and minikraken database."
echo "Times will be saved into \"${host}_perf.csv\" file, output files into \"outputs\" dir."
echo " "

export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib


if ! which bc > /dev/null; then
        sudo apt-get install bc
fi

spinner() {
        local i sp n
        sp='/-\|'
        n=${#sp}
        echo -n ' '
        while sleep 0.1;
        do
                printf "\b${sp:i++%n:1}"
        done
}

upload_db_in_ram() {
        printf " * Uploading db on RAM, please wait... "
        spinner &
        cat ../minikraken_20141208/database.* > /dev/null
        PID=$!
        disown $PID
        #kill "$!" # kill the spinner
        kill -9 $PID # kill the spinner
        printf '\n'
}

mkdir -p outputs

pathDatasets=/root/fog_computing


dbName=${2}
nf=${3}
nFile=$(find "../${dbName}"_1_"${nf}" -name *.fast5 | wc -l)

maxTh=$1
##################################################################################
echo " * deepnano..."
for th in $(seq ${maxTh})
do
        var=$RANDOM
        dir=${pathDatasets}/${dbName}_${th}_${nFile}
        echo "   - Dir: ${dir}, files: ${nFile}, proc: ${th}"
        OMP_NUM_THREADS=1  THEANO_FLAGS=base_compiledir=/tmp/${var}/theano.NOBACK python  basecall_no_metrichor.py  --directory ${dir}  --output outputs/${th}_${nFile}_files_db500GB.fasta 1>deepnano_${th}_${nFile}.txt &
done
wait

echo " "
upload_db_in_ram

echo " * kraken-db..."
for th in $(seq ${maxTh})
do
        ../kraken-0.10.5-beta//kraken  --db ../minikraken_20141208/  outputs/${th}_${nFile}_files_db500GB.fasta> outputs/sequences_${th}_${nFile}.kraken 2>kdb_${th}_${nFile}.txt &
done
wait

echo " * kraken-lab..."
for th in $(seq ${maxTh})
do
        ../kraken-0.10.5-beta//kraken-translate  --db ../minikraken_20141208/  outputs/sequences_${th}_${nFile}.kraken > outputs/sequences_${th}_${nFile}.labels &
done
wait
##################################################################################

echo " * Computing performance..."

t_bp=0
t_deepnV_bps=0.00
t_kraken_bps=0.00
for th in $(seq ${maxTh})
do
        # tempo deeepnano solo sequenziamento senza tempo iniz rete neurale:
        cat deepnano_${th}_${nFile}.txt | grep vero > deepnano_${th}_${nFile}_tempo.txt
        time_deepV=$( cat deepnano_${th}_${nFile}_tempo.txt | awk '{print $1;}' )

        # n basi:
        sed -r 's/\r//g' kdb_${th}_${nFile}.txt > kdb_${th}_${nFile}_tmp.txt #tolgo ^M
        mv kdb_${th}_${nFile}_tmp.txt kdb_${th}_${nFile}.txt
        cat kdb_${th}_${nFile}.txt | grep processed > kdb_${th}_${nFile}_1r.txt # seleziono la prima riga
        bp=$( cat kdb_${th}_${nFile}_1r.txt | awk '{print $4;}' | sed -r 's/\(//g' ) # n di basi processate (bp)

        # bp/s deepnano:
        deepnV_bps=$( printf '%.*f\n' 2 $( echo $bp / $time_deepV | bc -l) )

        #Mbp/min kraken:
        MbpM=$( cat  kdb_${th}_${nFile}_1r.txt |sed 's/.*Kseq\/m, \([\.0-9]*\).*/\1/'|head -n1 )
        #MbpM=$( cat kdb_${th}_${nFile}_1r.txt | awk '{print $15;}' | sed -r 's/\)//g' )
        # bp/s kraken:
        kraken_bps=$( printf '%.*f\n' 2 $( echo "$MbpM * 1000000 / 60" | bc) )

        # incremento variabili:
        t_bp=$((${t_bp} + $bp))
        t_deepnV_bps=$(echo ${t_deepnV_bps} + $deepnV_bps | bc)
        t_kraken_bps=$(echo ${t_kraken_bps} + $kraken_bps | bc)

        rm  kdb_${th}_${nFile}.txt  deepnano_${th}_${nFile}.txt  deepnano_${th}_${nFile}_tempo.txt  kdb_${th}_${nFile}_1r.txt
        rm -r outputs
done

if [ ! -f ${host}_perf.csv ]; then
    echo -e "#host,\\tdbName,\\tmaxTh,\\tnFile,\\tn-bp,\\tt_deepnV_bps,\\tt_kraken_bps" >> ${host}_perf.csv
fi

echo -e "${host}\\t${dbName}\\t${maxTh}\\t${nFile}\\t${t_bp}\\t${t_deepnV_bps}\\t${t_kraken_bps}" >> ${host}_perf.csv
echo "/-----------------------------------------------------/"
echo "Done."
