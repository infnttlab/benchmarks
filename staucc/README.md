# STAUCC
### Il Repository:
In questo repository si trovano quattro cartelle:
- `staucc`: cartella in cui si trovano i codici
- `doc`: contiene i paper di riferimento
- `input`: cartella in cui si possono trovare gli zip di tre input con dimesioni rispettivamente 36, 81 e 4096 membrane
- `scripts`: script bash creati per automatizzare il calcolo dei benchmarks

### Il Codice:
Sono presenti tre versioni del codice:

1. sequenziale (nella cartella _staucc_seqpar/staudago/_)
2. parallela (nella cartella _pardago/_ dentro a quella sequenziale)
3. CUDA (nella cartella _staucc_cuda/_)

Nella cartella `staucc_seqpar/` si trovano le seguenti cose:
- _input/_
- _pardago/_
- `compila.sh`: 
- `io_file.c`
- `io_file.h`
- `main.c`
- `manage_prob.c`
- `manage_prob.h`
- `p_sim.h`
- `scrivioutput.cc`
- `valgrind.txt`

### Compilazione:
Se non sono ancora installate, installare le librerie gsl:
```sh
apt-cache search libgsl
apt-get install libgsl0-dev
```

Sequenziale (nella cartella `staucc_seqpar/staudago/`):
```sh
gcc -w -Wall -g -O2 -fno-inline-functions -o staudago_out -I/usr/include/gls -L/usr/lib64 *.c -ldl -lgsl -lgslcblas -lm
./staudago_out
```

Parallela (nella cartella `pardago/`):
```sh
mpicc -Wall -g -o pardago_out -I/usr/include/gls -L/usr/lib64 *.c -lgsl -lgslcblas -lm
mpirun -np 4 ./pardago_out
```

CUDA (nella cartella `staucc_cuda/new/`):
```sh
[/usr/local/cuda-6.5/bin/]nvcc  -I/usr/local/cuda-6.5/samples/common/inc -gencode arch=compute_20,code=sm_20 staudpp.cu -o staudpp_out
./staudpp_out
```
_NB_: I comandi per la compilazione si trovano nel file `compila.sh` nelle rispettive cartelle e in `staudpp.cu` per CUDA.

### Cambiare i parametri:
I parametri che possono essere variati sono:
- time-max (s) (`input/time_max`)
- every (`input/every` - il primo numero)

Per controllare il numero di mebrane di quel set di dati andare a vedere in `input/numMembranes.txt`.
