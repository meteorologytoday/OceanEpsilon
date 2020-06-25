#!/bin/bash

source 00_basic.sh

mkdir $data_dir

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
 


    julia tools/gather_data.jl                                  \
        --models       "$model"                                 \
        --domain-file  "${domain_dir}/domain.${model}.nc"       \
        --vars         "tauuo,tauvo"                            \
        --time-range   "${beg_year},$(( beg_year + years - 1 ))"

  
    if [ ]; then 
    julia tools/interpolate_model.jl                            \
        --input-dir    "cmip_data/${model}/wo"                  \
        --domain-file  "${domain_dir}/domain.${model}.nc"       \
        --output-file  "${data_dir}/input.${model}.nc"          \
        --var          "wo"                                     \
        --lev-target   50                                       \
        --time-range   "${beg_year},$(( beg_year + years - 1 ))"
    
     fi

done

