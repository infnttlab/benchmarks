#matrixMul

### Sequenziale
```sh
gcc -O3 matrixMul_seq.c -o mtxMul.seq.O3
```
__Uso:__
```sh
./mtxMul.seq.O3 <ROW_A> <COL_A> <COL_B>
```
_Dove:_ MATRIX(ROW, COL) e COL_A == ROW_B  
_Default:_ A = (512,512); B = (512,512)  
#### _Troubleshooting:_
**1. `Segmentation fault` con dimensioni >= 1000**  
Aumentare la dimensione dello stack: `ulimit -s 20000`

---
### CUDA
_NB:_ Se non trova il compilatore `nvcc` aggiungerci il rispettivo path (es. `/usr/local/cuda-x.y/bin/`)
```sh
nvcc -O3 matrixMul_cuda.cu -lm -o mtxMul.cuda.O3
```  
__Uso:__
```sh
./mtxMul.cuda.O3 <ROW_A> <COL_A> <COL_B> <DIM_BLOCK>
```
_Dove:_  DIM_BLOCK: [1-32]; BLOCK(dimBlock, dimBlock)  
_Default:_ DIM_BLOCK = 16

---
### OpenACC
* [Documentations](http://www.openacc.org/node/1)

_NB:_ The code to run needs CUDA libraries and gcc version >= 5.0
```sh
gcc matrixMul_oacc.c -fopenacc -foffload="-O3" -O3 -o mtxMul.oacc.O3
```
__Uso:__
```sh
./mtxMul.oacc.O3 <ROW_A> <COL_A> <COL_B>
```

---
### OMP
```sh
gcc -O3 -fopenmp matrixMul_omp.c -o mtxMul.omp.O3
```
__Uso:__
```sh
./mtxMul.omp.O3 <ROW_A> <COL_A> <COL_B> <THREADS>
```
_Default:_ THREADS = 2 
