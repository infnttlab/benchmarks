#!/bin/bash

if [[ ( $1 == "--help" ) || ( $1 == "-h" ) ]]
then
        echo -e "\\nUsage: $0  [n-max-core]  [n-file-tot]  [db-name]  [db-path]"
        echo "Where [db-path] is the full path before 'splitDB/' dir."
        echo -e "i.e: /mnt/avoton/fog-computing/data/test/splitDB/ ---> [db-path] = /mnt/avoton/fog-computing/data/test/\\n"
        exit 0
fi

if [[ ( -z "$3" ) ]]
then
        echo -e "\\nUsage: $0  [n-max-core]  [n-file-tot]  [db-name]  [db-path]"
        echo "Where [db-path] is the path before 'splitDB/' dir."
        echo -e "i.e: /mnt/avoton/fog-computing/data/test/splitDB/ ---> [db-path] = /mnt/avoton/fog-computing/data/test/\\n"
        exit 1
fi

maxcore=$1
nfile=$2
dir=${3%/}
pathDB=${4%/}

nth=1
while [ "${nth}" -le "${maxcore}" ]; do
        for num in $(seq $nth); do
                file_in=$((nfile / nth))
                dest="${pathDB}/splitDB/split${nth}/${dir}_${num}_${file_in}/"
                cp -rL ${dest} ../
        done

        ./run_deepnano_kraken.sh  ${nth}  ${dir}  ${file_in}

        nth=$((nth * 2))
done

rm -rf  ../${dir}_*
