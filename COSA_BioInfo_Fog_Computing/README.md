## Install Deepnano and Minikraken DB:

Run the script `install_deepnano_kraken.sh`

```sh
./install_deepnano_kraken.sh
```

At the end of its executiond will be created a `fog_computing/` directory and into this:
- `deepnano/` forlder that contains all script to run Deepnano software and the scripts that you can use to do tests.
- `jellyfish-1.1.11/` folder, it's a necessary library
- `kraken-0.10.5-beta/`
- `minikraken_20141208/`

## Dataset Used:

Abbiamo selezionato 32 file da ciascun dataset per formare un mini-dataset e processare quello scalando sui core.

| DB           |   DB size  | n file tot |  Size tot.  |
| ------------ | :--------: | :--------: | :---------: |
| Ecoli_2D     |    2.1G    |     32     |     50M     |
| Equal_v5_2D  |    1.2G    |     32     |     54M     |
| Equal_v6_2D  |    2.2G    |     32     |     46M     |
| Maeru_2D     |    1.1G    |     32     |     63M     |
| Pfluor_2D    |    1.8G    |     32     |     51M     |
| Rare_v6_2D   |    1.6G    |     32     |     53M     |
| Selong_2D    |    367M    |     32     |     62M     |
| Staggered_2D |    4.8G    |     32     |     44M     |

## Schede usate:

| Host     |   N.core  | 
| -------- | :-------: | 
| xeond4   |     8     | 
| avoton2  |     8     | 
| J4205-01 |     4     | 
| n3700-01 |     4     | 
| n3710-03 |     4     | 


```sh
    
#DB used:
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Ecoli_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Equal_v5_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Equal_v6_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Maeru_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Pfluor_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Rare_v6_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Selong_2D.tar.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/Staggered_2D.tar.gz

cd data/
ls -1
    
dataset4test
Ecoli_2D
Equal_v5_2D
Equal_v6_2D
Maeru_2D
Pfluor_2D
Rare_v6_2D
readme.txt
Selong_2D
split_db.sh
split_db.sh
Staggered_2D
    
# Lo script di seguito crea un db di 32 file provenienti dal db Equal_v5_2D (in questo caso)
# questo db verrà poi diviso in sotto cartelle contenenti un numero uguale di file 
# in modo tale che si possa processare il db usando un numero di core crescente.
# Queste sotto cartelle sono raggruppate in cartelle (contenute in dataset4test/) 
# in base al numero di core usati, quindi:
#
# nomeDir   nSubDir   nFileXdir
# th1        1           32        
# th2        2           16
# th4        4           8
# th8        8           4

./split_db.sh 8 32 Equal_v5_2D/
    
# Ora si va sulla macchina sulla quale far girare il test (es. xeond4)
cd fog_computing/deepnano
    
root@xeond4:~/fog_computing/deepnano# ls -1

align_2d
align_2d.cc
basecall_no_metrichor_devel.py
basecall_no_metrichor.py
basecall.py
helpers.py
helpers.pyc
LICENSE
nets_data
r9
README.md
rnn_fin.py
rnn_fin.pyc
run_all_test.sh
run_deepnano_kraken.sh
run_single_db.sh
scheduler.py
training
    
# Si fa girare "run-all-test.sh" che a sua volta chiama "run-single-db.sh" che chiama "running_in_the_fog_MCORE_flag.sh"
./run-all-test.sh
    
# run-single-db.sh: ./run-single-db.sh  [n-max-core]  [n-file-tot]  [dir]
# es: ./run-single-db.sh 8 32 Ecoli_2D
#
# running_in_the_fog_MCORE_flag.sh: ./running_in_the_fog_MCORE_flag.sh  [n-proc]  [db-name] [n-file]
# es: ./running_in_the_fog_MCORE_flag.sh  8  Ecoli_2D   4
    
# verrà creato il file "performance_bps_MCORE.txt" con le performance, i cui dati saranno ordinati secondo lo schema:
# host    dbName  maxTh   nFile   n-bp    t_deepnV [bp/s] t_kraken [bp/s]
```
