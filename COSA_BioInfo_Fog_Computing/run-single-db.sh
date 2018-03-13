#!/bin/bash

if [[ ( $1 == "--help" ) || ( $1 == "-h" ) ]]
then
        echo -e "\\nUsage: $0  [n-max-core]  [n-file-tot]  [dir]\\n"
        exit 0
fi

if [[ ( -z "$3" ) ]]
then
        echo -e "\\nUsage: $0  [n-max-core]  [n-file-tot]  [dir]\\n"
        exit 1
fi

maxcore=$1
nfile=$2
dir=${3%/}

#list_file=($(ls -1 ${dir} | head -n ${nfile}))

#TODO assicurarsi che list_file sia lungo proprio nfile

nth=1
while [ "${nth}" -le "${maxcore}" ]; do
        for num in $(seq $nth); do
                file_in=$((nfile / nth))
                dest="/mnt/avoton/fog/data/dataset4test/th${nth}/${dir}_${num}_${file_in}/"
                cp -rL ${dest} ../
        done

        ./run_deepnano_kraken.sh  ${nth}  ${dir}  ${file_in}

        nth=$((nth * 2))
done

rm -rf  ../${dir}_*
