#matrixMul

### Sequenziale
```sh
gcc matrixMul_seq.c 
```
Uso:
```sh
./a.out <ROW_A> <COL_A> <COL_B>
```
Dove: MATRIX(ROW, COL) e COL_A == ROW_B  
Default: A = (512,512); B = (512,512)  
#### _Troubleshooting:_
**1. `Segmentation fault` con dimensioni >= 1000**  
Aumentare la dimensione dello stack: `ulimit -s 20000`

---
### CUDA
```sh
/usr/local/cuda-6.5/bin/nvcc matrixMul_cuda.cu -lm -o a.out
```
Uso:
```sh
./a.out <ROW_A> <COL_A> <COL_B> <DIM_BLOCK>
```
Dove:  DIM_BLOCK: [1-32]; BLOCK(dimBlock, dimBlock)  
Default: DIM_BLOCK = 16
