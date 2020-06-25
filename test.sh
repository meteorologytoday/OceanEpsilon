#!/bin/bash


data_root_dir=/glade/collections/cmip/CMIP6/CMIP

julia interpolate_model.jl \
    --input-file   "$data_root_dir/NCAR/CESM2/piControl/r1i1p1f1/Omon/wo/gn/latest/wo_Omon_CESM2_piControl_r1i1p1f1_gn_010001-019912.nc" \
    --output-file  "NCAR_wo.nc" \
    --var          "wo"         \
    --lev          "lev"        \
    --lev-unit     "cm"         \
    --lev-target   50           \
    --time-range   1,120



