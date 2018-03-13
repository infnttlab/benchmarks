#!/bin/sh

#db: Ecoli_2D Equal_v5_2D Equal_v6_2D Maeru_2D Pfluor_2D Rare_v6_2D Selong_2D Staggered_2D

for db in Ecoli_2D Equal_v5_2D Equal_v6_2D Maeru_2D Pfluor_2D Rare_v6_2D Selong_2D Staggered_2D
do
        ./run-single-db.sh 8 32 ${db}
done
