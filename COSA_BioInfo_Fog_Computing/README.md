## How to install Deepnano and Minikraken DB

Run the script `install_deepnano_kraken.sh`

```sh
./install_deepnano_kraken.sh
```

A folder named `fog_computing/` will be created, which contains:
- `deepnano/` with all the scripts needed to run Deepnano software do tests.
- `jellyfish-1.1.11/` with a needed library
- `kraken-0.10.5-beta/`
- `minikraken_20141208/`

## How to run
### Deepnano:

```sh
OMP_NUM_THREADS=1 python  basecall_no_metrichor.py  --directory  [dir-fast5-file]  --output [output-name].fasta
```
`basecall_no_metrichor.py`, originally from Deepnano, has been modified in order to output execution time.

### Kraken:

```sh
./kraken --db ../minikraken_20141208/ [deepnano-output].fasta
```
### Dedicated scripts:

```sh
# To execute Deepnano and then Kraken on n_proc cores simultaneously (n_proc must be equal to the number of folders in DB)
# Folder names in DB must be: <db-name>_<index>_<n-file-in-this-dir>
# Example: if you want to use four processes to process your dataset with 32 files in total,
# you have to split it in four folders in which there are 8 files each. So you'll have:
# Ecoli_2D_1_8/; Ecoli_2D_2_8/; Ecoli_2D_3_8/; Ecoli_2D_4_8/ 
./run_deepnano_kraken.sh  [n-proc]  [db-name] [n-file]
Example: ./run_deepnano_kraken.sh  8  Ecoli_2D   4

# Scheduler for executing multiple instances of Deepnano and Kraken in parallel, feeding a queue and preparing .js file for AWS
python scheduler.py -path INPUT_PATH -dataset DATASET_NAME -n MAX_NUM_CORE
Example: python scheduler.py -p /mnt/avoton/fog/data/prova/splitDB/split8 -d Ecoli -n 8 

# To create a mini-database and to further split it in folders.
# n-file files are equally distributed among different sub-folders,
# which are then grouped in different folders according to the number of cores
# used in the tests, i.e. th1, th2, th4, th8 for 1, 2, 4, 8 cores respectively.
# The structure in /dataset4test is the following:
#
# folderName   nSubFolder   nFileXFolder
# th1            1           32        
# th2            2           16
# th4            4           8
# th8            8           4
./split_db.sh  [dest-path]   [n-max-core]  [n-file]  [db-name]
Example: ./split_db.sh /mnt/avoton/fog/data/tt/ 8 32 Equal_v5_2D

```

## Our tests:
### Dataset:

The mini-dataset we processed with our SoCs is composed of 32 files selected from each of the following datasets:

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

which can be downloaded from ftp://ftp.sra.ebi.ac.uk/vol1/ERA742/ERA742349/oxfordnanopore_native/


### SoCs:

| Host     |   N.core  | 
| -------- | :-------: | 
| xeond4   |     8     | 
| avoton2  |     8     | 
| J4205-01 |     4     | 
| n3700-01 |     4     | 
| n3710-03 |     4     | 


```sh
    
### Example of usage:
# Use split_db.sh to create a mini-database of 32 files for Equal_v5_2D.
./split_db.sh 8 32 Equal_v5_2D/
    
# Change working directory to where the applications are installed, e.g.
cd fog_computing/deepnano

# Execute 4 instances of Deepnano -> Kraken and then upload to AWS:
./run_deepnano_kraken.sh  4  Equal_v5_2D   4

# Execute 4 instances of Deepnano -> Kraken and then upload to AWS:
./python scheduler.py -p /mnt/avoton/fog/data/splitDB/split16 -d Equal_v5_2D -n 4
```
