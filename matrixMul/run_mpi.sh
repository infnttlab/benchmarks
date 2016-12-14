#!/bin/bash
echo " "
for sizeMtx in 128 512 2048
do
        for warmup in 0 1
        do
                for proc in 2 4 8 16
                do
                        if [ ${warmup} -eq 0 -a ${proc} -eq 2 ]; then
                                echo mpirun -np ${proc} --allow-run-as-root --mca btl ^openib ./mtxMul.mpi.O3 -rA=${sizeMtx} -cA=${sizeMtx} -cB=${sizeMtx} -w=${warmup}
                        fi
                        mpirun -np ${proc} --allow-run-as-root --mca btl ^openib ./mtxMul.mpi.O3 -rA=${sizeMtx} -cA=${sizeMtx} -cB=${sizeMtx} -w=${warmup}
                done
                echo " "
        done
        echo "---------------------------------------------------------"
done
