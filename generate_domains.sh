#!/bin/bash

data_root_dir=/glade/collections/cmip/CMIP6/CMIP

infos=(
    NCAR_CESM2 NCAR/CESM2/piControl/r1i1p1f1/Omon/wo/gn/latest/wo_Omon_CESM2_piControl_r1i1p1f1_gn_010001-019912.nc                              4 lat lat_bnds lon lon_bnds lev lev_bnds 
    EC_Earth   EC-Earth-Consortium/EC-Earth3/piControl/r1i1p1f1/Omon/wo/gn/v20190712/wo/wo_Omon_EC-Earth3_piControl_r1i1p1f1_gn_262601-262612.nc 4 latitude vertices_latitude longitude vertices_longitude lev lev_bnds
    E3SM       E3SM-Project/E3SM-1-0/piControl/r1i1p1f1/Omon/wo/gr/v20191007/wo/wo_Omon_E3SM-1-0_piControl_r1i1p1f1_gr_045601-046012.nc          2 lat lat_bnds lon lon_bnds lev lev_bnds         
)

domain_dir=domains

mkdir $domain_dir

for i in $(seq 1 $((${#infos[@]}/9))); do
    model=${infos[$((9*(i-1) + 0))]}
    example_file=${infos[$((9*(i-1) + 1))]}
    h_vertices=${infos[$((9*(i-1) + 2))]}
    lat_c=${infos[$((9*(i-1) + 3))]}
    lat_v=${infos[$((9*(i-1) + 4))]}
    lon_c=${infos[$((9*(i-1) + 5))]}
    lon_v=${infos[$((9*(i-1) + 6))]}
    lev_c=${infos[$((9*(i-1) + 7))]}
    lev_b=${infos[$((9*(i-1) + 8))]}

    echo "Doing model: $model"


    julia tools/generate_domain_file.jl \
        --input-file "${data_root_dir}/${example_file}" \
        --output-file "${domain_dir}/domain.${model}.nc" \
        --vertices $h_vertices \
        --lat-c $lat_c \
        --lat-v $lat_v \
        --lon-c $lon_c \
        --lon-v $lon_v \
        --lev-c $lev_c \
        --lev-b $lev_b

done
