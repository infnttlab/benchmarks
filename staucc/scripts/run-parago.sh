#!/bin/bash

display_usage(){
        echo -e "\nUsage: bash $0  <COMPIL-VERSION>  <N-PROC-PAR>  <N-IT>  <MEMB>  <T-MAX>  <EVERY>
Es. bash $0 gcc4.8.2 4 10 81 1 100\n"
}


if [[ ( $1 == "--help" ) || $1 == "-h" || ( -z $1 ) || ( -z $2 ) || ( -z $3 ) || ( -z $4 ) || ( -z $5 ) || ( -z $6 ) ]]
then
        display_usage
else

        mpicc -Wall  -w -g -o no_out -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
        mpicc -Wall -O2 -w -g -o O2_out -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
        mpicc -Wall -O3 -w -g -o O3_out -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
        mpicc -Wall -Ofast -w -g -o Ofast_out -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm

        iterations=$3
        membrane=$4
        timeMax=$5
        every=$6
        vgcc=$1
        procPar=$2

        for exe in *_out
        do
                echo "### Running" $exe "..."

                for i in `seq 1 $iterations`;
                do
                        START=$(date +%s.%3N)
                        mpirun -np $procPar ./${exe}
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
                echo "$HOSTNAME $vgcc $membrane $timeMax $every $procPar $exe $av_times $iterations" >> times_parago.txt

                rm -r output/ $exe
        done
fi
