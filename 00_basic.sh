#!/bin/bash

data_root_dir=/glade/collections/cmip/CMIP6/CMIP

infos=(
    MOHC.HadGEM3-GC31-LL            1850  4 A latitude vertices_latitude longitude vertices_longitude lev lev_bnds m normY stress_normal
    CNRM-CERFACS.CNRM-ESM2-1        1850  4 A lat bounds_lat lon bounds_lon lev lev_bounds m normY stress_normal
    MPI.ESM1-2-LR                   1900  4 C latitude vertices_latitude longitude vertices_longitude lev lev_bnds m flipY stress_normal
    EC-Earth-Consortium.EC-Earth3   2260  4 A latitude vertices_latitude longitude vertices_longitude lev lev_bnds m normY stress_normal
    NCAR.CESM2                         1  4 A lat lat_bnds lon lon_bnds lev lev_bnds cm normY stress_normal
    E3SM.1-1                        1850  2 A lat lat_bnds lon lon_bnds lev lev_bnds m normY stress_normal
    CanESM5                         5201  4 B latitude vertices_latitude longitude vertices_longitude lev lev_bnds m normY stress_normal
    SNU.SAM0-UNICORN                   1  4 A latitude vertices_latitude longitude vertices_longitude lev lev_bnds m normY flip_stress
    UA.MCM-UA-1-0                      1  2 A latitude lat_bnds longitude lon_bnds lev lev_bnds m normY stress_normal
)

# MRI does not have the same Ny in wo and tauuo/tauvo    
# MRI.ESM2-0       wo_Omon_MRI-ESM2-0_piControl_r1i1p1f1_gn_185001-189912.nc  1850  4 latitude vertices_latitude longitude vertices_longitude lev lev_bnds m

years=10

domain_dir=domains
data_dir=data
output_dir=result
tmp_dir=tmp
