#!/bin/bash

rm -f mtxMul.seq.sub.O3
g++ -I../../common/inc -mfloat-abi=hard -O3 -o matrixMul-seq.o -c matrixMul_seq_submtx.c
g++ -O3 -mfloat-abi=hard -Xlinker --dynamic-linker=/lib/ld-linux-armhf.so.3 -o mtxMul.seq.sub.O3 matrixMul-seq.o
rm -f matrixMul-seq.o