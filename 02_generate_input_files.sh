#!/bin/bash

source 00_basic.sh

mkdir $data_dir

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


    file1="${data_dir}/input.${model}.tauo.nc"
    file2="${data_dir}/input.${model}.wo.nc"
    file_final="${data_dir}/input.${model}.nc"

    {
    julia tools/gather_data.jl                                  \
        --model        "$model"                                 \
        --domain-file  "${domain_dir}/domain.${model}.nc"       \
        --output-file  "$file1"                                 \
        --vars         "tauuo,tauvo"                            \
        --time-range   "${beg_year},$(( beg_year + years - 1 ))"

    if [ "$stress_type" = "flip_stress" ] ; then
        ncap2 -O \
            -s 'tauuo=-tauuo' \
            -s 'tauvo=-tauvo' \
            $file1 $file1
    fi
    } &

    {
    julia tools/interpolate_model.jl                            \
        --input-dir    "cmip_data/${model}/wo"                  \
        --domain-file  "${domain_dir}/domain.${model}.nc"       \
        --output-file  "$file2"                                 \
        --var          "wo"                                     \
        --lev-target   50                                       \
        --time-range   "${beg_year},$(( beg_year + years - 1 ))"
    }&

    wait

    ncks -A -v tauuo,tauvo $file1 $file2

    rm $file1

    if [ "$y_axis_type" = "flipY" ] ; then
        ncpdq -O -a '-Ny' $file2 $file2
    fi
    
    mv $file2 $file_final

done

