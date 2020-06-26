#!/bin/bash

source 00_basic.sh

mkdir $domain_dir

for i in $(seq 1 $((${#infos[@]}/13))); do
    model=${infos[$((13*(i-1) + 0))]}
    beg_year=${infos[$((13*(i-1) + 1))]}
    h_vertices=${infos[$((13*(i-1) + 2))]}
    h_vertices_type=${infos[$((13*(i-1) + 3))]}
    lat_c=${infos[$((13*(i-1) + 4))]}
    lat_v=${infos[$((13*(i-1) + 5))]}
    lon_c=${infos[$((13*(i-1) + 6))]}
    lon_v=${infos[$((13*(i-1) + 7))]}
    lev_c=${infos[$((13*(i-1) + 8))]}
    lev_b=${infos[$((13*(i-1) + 9))]}
    lev_unit=${infos[$((13*(i-1) + 10))]}
    y_axis_type=${infos[$((13*(i-1) + 11))]}
    stress_type=${infos[$((13*(i-1) + 12))]}

    echo "Doing model: $model"

    p="cmip_data/${model}/wo"
    example_file=$( eval "ls -v $p/*wo*.nc" | head -n 1 )
        
    output_file="${domain_dir}/domain.${model}.nc"

    julia tools/generate_domain_file.jl  \
        --input-file "${example_file}"   \
        --output-file ${output_file}     \
        --vertices $h_vertices           \
        --vertices-type $h_vertices_type \
        --lat-c $lat_c \
        --lat-v $lat_v \
        --lon-c $lon_c \
        --lon-v $lon_v \
        --lev-c $lev_c \
        --lev-b $lev_b \
        --lev-unit $lev_unit 

    if [ "$y_axis_type" = "flipY" ] ; then
        echo "Flip y axis..."
        ncpdq -O -a '-Ny' $output_file $output_file
    fi
    

done
