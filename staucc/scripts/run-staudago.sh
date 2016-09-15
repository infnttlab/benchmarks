#!/bin/bash

display_usage(){
        echo -e "\nUsage: bash $0 <COMPIL-VERSION>  <N-IT>  <N-MEM>  <T-MAX>  <EVERY>
Es. bash $0 gcc4.8.2 10 81 1 100\n"
}


if [[ ( $1 == "--help" ) || $1 == "-h" || ( -z $1 ) || ( -z $2 ) || ( -z $3 ) || ( -z $4 ) || ( -z $5 ) ]]
then
        display_usage
else

        #gcc -w -Wall -g -fno-inline-functions -o no_out -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm
        gcc -w -Wall -g -O2 -fno-inline-functions -o O2_out -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm
        gcc -w -Wall -g -O3 -fno-inline-functions -o O3_out -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm
        gcc -w -Wall -g -Ofast -fno-inline-functions -o Ofast_out -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm

        ls -lhtr

        iterations=$2
        membrane=$3
        timeMax=$4
        every=$5
        vgcc=$1

        for exe in *_out
        do
                echo "### Running" $exe "..."

                for i in `seq 1 $iterations`;
                do
                        START=$(date +%s.%3N)
                        ./${exe}
                        END=$(date +%s.%3N)

                        parrec_times[i]=$(echo "$END - $START" | bc)
                        # echo "${parrec_times[i]}" >> ../output/allTimes_${nameEXE}_1slice_${iterations}it.txt
                done

                #TEMPI:
                sum_t=0
                #um_qt=0
                for k in "${parrec_times[@]}"
                do
                        sum_t=$(echo "$sum_t+$k" | bc -l)
                done

                # media di tutti gli N tempi:
                av_times=$(echo "$sum_t/$iterations" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')

                # errore sulla media (deviazione standard):
                #er_av_t=$(echo "sqrt($rad_t)" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')

                #------------------------------------------------------
                #SALVARE SUL FILE TXT:

                av_times=${av_times/./,}
                echo "$HOSTNAME $vgcc $membrane $timeMax $every $exe $av_times $iterations" >> times_staudago.txt

                rm -r output/ $exe
        done
fi
