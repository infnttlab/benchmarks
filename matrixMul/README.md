#matrixMul

Benchmarks: [Google Drive Benchmarks](https://docs.google.com/spreadsheets/d/1j6MXqHPtD__qOVB4mXoK4K4jAJX0GHGzMdnMNwUD_-M/edit#gid=0)  
Nella cartella deve essere presente il file: `helper_string.h`.

### Sequenziale
#### A. Vanilla
```sh
g++ -O3 matrixMul_seq.c -o mtxMul.seq.O3

Usage:   -rA=RowsA     -cA=ColumnsA  -cB=ColumnsB | matrix(row,col), ColumnsA = RowsB
         -w=WarmUpData
         -v=Verbose
```
_Default_: A = (512,512) B = (512,512); WARMUP = 0; VERBOSE = 0  

_Troubleshooting:_  
**1.** `Segmentation fault` con dimensioni >= 1000  
Aumentare la dimensione dello stack: `ulimit -s 20000`

#### B. Sub-matrix
Per compilare:
```sh
./compile_seq_sub.sh

Usage: -size=SizeMatrix  -subSize=SizeSubMatrix
       -perf=Performance -v=Verbose

Note:  Size matrix must be a multiple of the size sub-matrix.
```
_Default_: SizeMatrix = 512; SizeSubMatrix = 2  

#### C. BLAS lib
Installazione:  
[Link-BLAS](http://www.netlib.org/blas/)  
[Link-tutorial-install](http://matrixprogramming.com/2008/01/matrixmultiply)  
```sh
# la libreria è stata già installata sulla jetson-k1-01 e si trova in /root/BLAS-3.5.0
wget http://www.netlib.org/blas/blas.tgz
tar zxvf blas.tgz
cd BLAS-3.5.0
make
# error: gfortran: Command not found
# apt-get install gfortran
mv blas_LINUX.a libblas.a #oppure: ar rv libblas.a *.o
```
Nel codice bisogna aggiungerci:  
```sh
#include <stdio.h>  
#include <stdlib.h>
```

Compilare:
```sh
g++ matrixMul_seq_blas.c -L/root/BLAS-3.5.0/libblas.a -lblas -O3 -o mtxMul.seq.blas.O3

# error: /usr/bin/ld: cannot find -lblas
# apt-get install libblas-dev

Usage: -size=SizeMatrix(d:512)
       -perf=Performance(d:0)   -v=Verbose(d:0)
```
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




