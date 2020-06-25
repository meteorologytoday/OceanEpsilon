using Formatting
using NCDatasets
using ArgParse
using JSON

function adjustLev(
    lev  :: AbstractArray,
    unit :: String,
)
    factor = Dict(
        "m"  => 1.0,
        "cm" => .01,
    )[unit]

    # check if lev is all positive or all negative
    if ! ( all(lev .>= 0) || all(lev .<= 0) ) 
        println("Warning: Not all levs are positive or negative")
    end


    if all(lev .<= 0)
        factor *= -1
    end
    
    lev = lev * factor

    Δlev = lev[2:end] - lev[1:end-1]

    if any(Δlev .<= 0)
        throw(ErrorException("Lev is not monitically increasing with indices."))
    end

    return lev
end


#include("PolelikeCoordinate.jl")

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Input file."
            arg_type = String
            required = true
 
        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
        "--vertices"
            help = "Number of vertices of each horizontal grid"
            arg_type = Int64
            required = true
            
        "--lat-v"
            help = "Varname of latitude on grid vertices."
            arg_type = String
            required = true

        "--lon-v"
            help = "Varname of longitude on grid vertices."
            arg_type = String
            required = true

        "--lev-b"
            help = "Varname of longitude on grid vertices."
            arg_type = String
            required = true


        "--lat-c"
            help = "Varname of latitude on grid center."
            arg_type = String
            required = true

        "--lon-c"
            help = "Varname of longitude on grid center."
            arg_type = String
            required = true

        "--lev-c"
            help = "Varname of lev on grid center."
            arg_type = String
            required = true

        "--lev-unit"
            help = "Unit of lev. Can be `m` or `cm`."
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))



Dataset(parsed["input-file"], "r") do ds

    global lev_b = nomissing(ds[parsed["lev-b"]][:])
    global lev_c = nomissing(ds[parsed["lev-c"]][:])
  

    lev_c = adjustLev(lev_c, parsed["lev-unit"])
    lev_b[1, :] = adjustLev(lev_b[1, :], parsed["lev-unit"])
    lev_b[2, :] = adjustLev(lev_b[2, :], parsed["lev-unit"])


 
    global Nz = length(lev_c)
    global Nb = size(lev_b)[1]

    if Nb != 2
        throw(ErrorException("Unexpected value of Nb"))
    end
 
    if parsed["vertices"] == 4

        global lon_v = nomissing(ds[parsed["lon-v"]][:])
        global lat_v = nomissing(ds[parsed["lat-v"]][:])

        global lon_c = nomissing(ds[parsed["lon-c"]][:])
        global lat_c = nomissing(ds[parsed["lat-c"]][:])

    
    elseif parsed["vertices"] == 2

        _lon_bnd = nomissing(ds[parsed["lon-v"]][:])
        _lat_bnd = nomissing(ds[parsed["lat-v"]][:])
        
        _lon_c = nomissing(ds[parsed["lon-c"]][:])
        _lat_c = nomissing(ds[parsed["lat-c"]][:])
        
        Nx = length(_lon_c)
        Ny = length(_lat_c)

        global lon_v = zeros(Float64, 4, Nx, Ny)
        global lat_v = similar(lon_v)

        if any( ( _lat_bnd[1, 2:end] - _lat_bnd[1, 1:end-1] ) .< 0 )

            throw(ErrorException("Error: latitude not increasing monotonically"))

        elseif any( ( _lon_bnd[1, 2:end] - _lon_bnd[1, 1:end-1] ) .< 0 )

            throw(ErrorException("Error: longitude not increasing monotonically"))

        end

        for i=1:Nx, j=1:Ny

            lon_v[1, i, j] = _lon_bnd[1, i]
            lon_v[2, i, j] = _lon_bnd[2, i]
            lon_v[3, i, j] = _lon_bnd[2, i]
            lon_v[4, i, j] = _lon_bnd[1, i]

            lat_v[1, i, j] = _lat_bnd[1, j]
            lat_v[2, i, j] = _lat_bnd[1, j]
            lat_v[3, i, j] = _lat_bnd[2, j]
            lat_v[4, i, j] = _lat_bnd[2, j]

        end



        global lon_c = repeat(reshape(_lon_c, :, 1), outer=(1, Ny))
        global lat_c = repeat(reshape(_lat_c, 1, :), outer=(Nx, 1))

    else
        throw(ErrorException("Unexpected vertices: " * string(Nv)))
    end



    global Nv, Nx, Ny = size(lon_v)
    
end 

#=
gi = CurvilinearSphericalGridInfo(;
    R = 6371e3,
    Ω = 2π/(86400.0 * 364 / 365),
    Nx=Nx,
    Ny=Ny,
    c_lon=lon_c,
    c_lat=lat_c,
    vs_lon=lon_v,
    vs_lat=lat_v,
    angle_unit=:deg,
)
=#         

Dataset(parsed["output-file"], "c") do ods

    defDim(ods, "Nx", Nx)
    defDim(ods, "Ny", Ny)
    defDim(ods, "Nv", Nv)
    defDim(ods, "Nz", Nz)
    defDim(ods, "Nb", Nb)
    defDim(ods, "time", Inf)

    for (varname, vardata, vardim, attrib) in [
        ("lon_v", lon_v, ("Nv", "Nx", "Ny"), Dict()),
        ("lat_v", lat_v, ("Nv", "Nx", "Ny"), Dict()),
        ("lon_c", lon_c, ("Nx", "Ny"), Dict()),
        ("lat_c", lat_c, ("Nx", "Ny"), Dict()),
        ("lev_b", lev_b, ("Nb", "Nz"), Dict()),
        ("lev_c", lev_c, ("Nz", ), Dict()),
    ]
        println("Doing varname:", varname)
        var = defVar(ods, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20

        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata
    end
end


