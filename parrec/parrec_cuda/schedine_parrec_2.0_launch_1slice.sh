#!/bin/bash

#cpufreq-set --governor performance
#echo 852000000 > /sys/kernel/debug/clock/override.gbus/rate

size_slice=${BASH_ARGV[*]}

#module load compilers/cuda-5.5
/usr/local/cuda-6.0/bin/nvcc -O3 -arch=sm_32 -o parrec_2.0_${size_slice}.out *.cu -lm -lcufft

#rm -f *.txt

start=400
end=400
iterations=5


#----------------------------------
#INIZIALIZZAZIONE DEI 2 FILE TXT:

#primo file txt (riassuntivo di tutto)
#echo "size_slice = ${size_slice}    start = $start    end = $end   iterations = $iterations" >> schedine_times_parrec2.0_${size_slice}_1slice_${iterations}it.txt
#echo "[numBlock  threadsPerBlock  dimBlock  av_sigma  er_sigma  av_times(s)  er_times(s)]" >> schedine_times_parrec2.0_${size_slice}_1slice_${iterations}it.txt
#echo "" >> schedine_times_parrec2.0_${size_slice}_1slice_${iterations}it.txt


#secondo file txt contentente tutti i tempi di tutte le iterazioni con le relative sigma
#echo "size_slice = ${size_slice}    start = $start    end = $end   iterations = $iterations" >> schedine_all_values_parrec2.0_${size_slice}_1slice_${iterations}it.txt
#echo "[numBlock  threadsPerBlock  dimBlock  av_sigma  er_sigma  av_times(s)  er_times(s)]" >> schedine_all_values_parrec2.0_${size_slice}_1slice_${iterations}it.txt
#echo "" >> schedine_all_values_parrec2.0_${size_slice}_1slice_${iterations}it.txt


#----------------------------------
#CORPO DEL PROGRAMMA:
for numBlock in 16 
do
	for threadsPerBlock in 32
	do
		for dimBlock in 16
		do
			echo "$numBlock   $threadsPerBlock   $dimBlock" >> schedine_all_values_parrec2.0_${size_slice}_1slice_${iterations}it.txt
			for i in `seq 1 $iterations`;
			do
				START=$(date +%s.%3N)
				./parrec_2.0_${size_slice}.out $start $end $numBlock $threadsPerBlock $dimBlock
				END=$(date +%s.%3N)
				parrec_times[i]=$(echo "$END - $START" | bc)
				sigma[i]=`python plot_differences_4bash_schedine.py recobj_ok_400.spr slice_reference_${size_slice}_400.sdt recobj_ok_400.sdt`
				echo "                                           ${sigma[i]}   ${parrec_times[i]}" >> schedine_all_values_parrec2.0_${size_slice}_1slice_${iterations}it.txt
			done
			
			#------------------------------------------------------------------------------------------------------
			# Calcolo la MEDIA e la DEVIAZIONE STANDARD(cioè l'errore sulla media) sia per i tempi sia per la sigma
			#------------------------------------------------------------------------------------------------------

			#TEMPI:
			
			sum_t=0
			sum_qt=0
			for k in "${parrec_times[@]}"
			do
				sum_t=$(echo "$sum_t+$k" | bc -l)
			        sum_qt=$(echo "$sum_qt+($k*$k)" | bc -l)
			done
			rad_t=$(echo "$sum_qt/$iterations-($sum_t/$iterations)*($sum_t/$iterations)" | bc -l)

			# media di tutti gli N tempi:
			av_times=$(echo "$sum_t/$iterations" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')

			# errore sulla media (deviazione standard):
			er_av_t=$(echo "sqrt($rad_t)" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')

			#------------------------------------------------------
			#SIGMA:
			
			sum_s=0
                        sum_qs=0
                        for j in "${sigma[@]}"
                        do
                                sum_s=$(echo "$sum_s+$j" | bc -l)
                                sum_qs=$(echo "$sum_qs+($j*$j)" | bc -l)
                        done
                        rad_s=$(echo "$sum_qs/$iterations-($sum_s/$iterations)*($sum_s/$iterations)" | bc -l)

                        # media di tutte le N sigma:
                        av_sigma=$(echo "$sum_s/$iterations" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')

                        # errore sulla media (deviazione standard):
                        er_av_s=$(echo "sqrt($rad_s)" | bc -l | awk -F"," -v OFS="," ' { for(i=0;NF- i++;){sub("[.]*0+ *$","",$i)};$1=$1 }1 '| sed 's/^\./0./')
			#------------------------------------------------------
			#SALVARE SUL FILE TXT:

			echo "$numBlock   $threadsPerBlock   $dimBlock   $av_sigma   $er_av_s   $av_times   $er_av_t" >> schedine_times_parrec2.0_${size_slice}_1slice_${iterations}it.txt

                done
        done
done

#rm parrec_2.0_${size_slice}.out

#mv *.txt ../../../schedine_risultati_test/


#-------------------------------------------------------------------------------------------------------------------------------------------




