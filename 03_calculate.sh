#!/bin/bash

source 00_basic.sh

mkdir $output_dir
mkdir $tmp_dir

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

    input_file="${data_dir}/input.${model}.nc"
    model_output_dir="${output_dir}/${model}"
    domain_file="${domain_dir}/domain.${model}.nc"

    mkdir $output_dir

    julia -p 6 tools/calculate_posterior.jl \
        --years       $years           \
        --input-file  "$input_file"    \
        --output-dir  "$model_output_dir"   \
        --domain-file "$domain_file"   

done
