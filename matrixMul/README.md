#matrixMul

Benchmarks: [Google Drive Benchmarks](https://docs.google.com/spreadsheets/d/1j6MXqHPtD__qOVB4mXoK4K4jAJX0GHGzMdnMNwUD_-M/edit#gid=0)
Nella cartella deve essere presente il file: `helper_string.h`.

### Sequenziale
```sh
g++ -O3 matrixMul_seq.c -o mtxMul.seq.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB
         -w=WarmUpData
         -v=Verbose
```
_Default_: A = (512,512) B = (512,512); WARMUP = 0; VERBOSE = 0  
#### _Troubleshooting:_
**1. `Segmentation fault` con dimensioni >= 1000**  
Aumentare la dimensione dello stack: `ulimit -s 20000`

---
### CUDA
_NB:_ Se non trova il compilatore `nvcc` aggiungerci il rispettivo path (es. `/usr/local/cuda-x.y/bin/`)
```sh
nvcc -O3 matrixMul_cuda.cu -lm -o mtxMul.cuda.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB
         -db=DimBlock                             | DimBlock(in threads): [1-32], block(DimBlock, DimBlock)
         -w=WarmUpData
         -v=Verbose
```  
_Default:_ DIM_BLOCK = 16

---
### OpenACC
* [Documentations](http://www.openacc.org/node/1)

_NB:_ The code to run needs CUDA libraries and gcc version >= 5.0
```sh
g++ matrixMul_oacc.c -fopenacc -foffload="-O3" -O3 -o mtxMul.oacc.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB
         -w=WarmUpData
         -v=Verbose
```

---
### OMP
```sh
g++ -O3 -fopenmp matrixMul_omp.c -o mtxMul.omp.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB
         -p=Threads
         -w=WarmUpData
         -v=Verbose
```
_Default:_ THREADS = 2 

---
### MPI
```sh
mpic++ -O3 matrixMul_mpi.c -o mtxMul.mpi.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB divisibile per num.processi
         -w=WarmUpData
         -v=Verbose
```




