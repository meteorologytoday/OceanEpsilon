#!/bin/bash

data_root_dir=/glade/collections/cmip/CMIP6/CMIP

infos=(
    NCAR_CESM2 wo_Omon_CESM2_piControl_r1i1p1f1_gn_010001-019912.nc      100  4 lat lat_bnds lon lon_bnds lev lev_bnds cm
    E3SM       wo_Omon_E3SM-1-0_piControl_r1i1p1f1_gr_045601-046012.nc   101  2 lat lat_bnds lon lon_bnds lev lev_bnds m
)
#    EC_Earth   wo_Omon_EC-Earth3_piControl_r1i1p1f1_gn_262601-262612.nc 2400  4 latitude vertices_latitude longitude vertices_longitude lev lev_bnds m

years=5

domain_dir=domains
data_dir=data
output_dir=result
tmp_dir=tmp
