#!/bin/bash

if [[ ( $1 == "--help" ) || ( $1 == "-h" ) ]]
then
        echo -e "\\nUsage: $0  [dest-path]   [n-max-core]  [n-file]  [db-name]\\n"
        exit 0
fi

if [[ ( -z "$4" ) ]]
then
        echo -e "\\nUsage: $0  [dest-path]   [n-max-core]  [n-file]  [db-name]\\n"
        exit 1
fi

path=${1%/}
maxcore=$2
nfile=$3
dir=${4%/}

list_file=($(ls -1 ${dir} | head -n ${nfile}))

#TODO assicurarsi che list_file sia lungo proprio nfile

nth=1
while [ "${nth}" -le "${maxcore}" ]; do
        for num in $(seq $nth); do
                file_in=$((nfile / nth))
                dest="${path}/splitDB/split${nth}/${dir}_${num}_${file_in}/"
                mkdir -p "${dest}"

                low=$(((num - 1) * file_in))
                high=$((low + file_in - 1))

                for i in $(seq $low $high); do
                        if [ -f "$PWD/${dir}/${list_file[$i]}" ]; then
                                cp  "$PWD/${dir}/${list_file[$i]}" "${dest}/${list_file[$i]}"
                        fi
                done
        done

        nth=$((nth * 2))
done
