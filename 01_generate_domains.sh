#!/bin/bash

source 00_basic.sh

mkdir $domain_dir

for i in $(seq 1 $((${#infos[@]}/11))); do
    model=${infos[$((11*(i-1) + 0))]}
    example_file=${infos[$((11*(i-1) + 1))]}
    beg_year=${infos[$((11*(i-1) + 2))]}
    h_vertices=${infos[$((11*(i-1) + 3))]}
    lat_c=${infos[$((11*(i-1) + 4))]}
    lat_v=${infos[$((11*(i-1) + 5))]}
    lon_c=${infos[$((11*(i-1) + 6))]}
    lon_v=${infos[$((11*(i-1) + 7))]}
    lev_c=${infos[$((11*(i-1) + 8))]}
    lev_b=${infos[$((11*(i-1) + 9))]}
    lev_unit=${infos[$((11*(i-1) + 10))]}

    echo "Doing model: $model"


    julia tools/generate_domain_file.jl \
        --input-file "cmip_data/${model}/wo/${example_file}" \
        --output-file "${domain_dir}/domain.${model}.nc" \
        --vertices $h_vertices \
        --lat-c $lat_c \
        --lat-v $lat_v \
        --lon-c $lon_c \
        --lon-v $lon_v \
        --lev-c $lev_c \
        --lev-b $lev_b \
        --lev-unit $lev_unit

done
