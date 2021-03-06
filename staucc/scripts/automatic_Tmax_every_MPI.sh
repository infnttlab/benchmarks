#!/bin/bash

display_usage(){
        echo -e "\nUsage: bash $0 <COMPIL-VERSION>  <N-IT>  <N-MEM>
Es. bash $0 gcc4.8.2 10 81\n"
}


if [[ ( $1 == "--help" ) || $1 == "-h" || ( -z $1 ) || ( -z $2 ) || ( -z $3 ) ]]
then
        display_usage
else
        iterations=$2
        membrane=$3
        vgcc=$1

        for Tmax in 1 10 100
        do
                echo cat input/time_max
                cat input/time_max
                echo $Tmax > input/time_max
                echo cat input/time_max
                cat input/time_max

                for every in 10 100 1000 10000 100000
                do
                        echo cat input/every
                        cat input/every
                        echo $every 10000 > input/every
                        echo cat input/every
                        cat input/every

                        for proc in 2 4 6 8 10 12 23 24
                        #for proc in 23
                        do
                                echo bash run-parago.sh  $vgcc  $proc  $iterations  $membrane  $Tmax  $every
                                bash run-parago.sh  $vgcc  $proc  $iterations  $membrane  $Tmax  $every
                        done
                done
        done
        echo -e "\nFINE DI TUTTO!!!! :)\n"
fi
