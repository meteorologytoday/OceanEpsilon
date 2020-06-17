using NCDatasets
using ArgParse
using JSON

function adjustLev(
    lev  :: AbstractArray{Float64, 1},
    unit :: String,
)
    factor = Dict(
        "m"  => 1.0,
        "cm" => .01,
    )[unit]

    # check if lev is all positive or all negative
    if ! ( all(lev .>= 0) || all(lev .<= 0) ) 
        throw(ErrorException("Not all levs are positive or negative"))
    end


    if all(lev .<= 0)
        factor *= -1
    end
    
    lev = lev * factor

    Δlev = lev[2:end] - lev[1:end-1]

    if any(Δlev .<= 0)
        throw(ErrorException("Lev is not monitically increasing with indices."))
    end

    return lev, Δlev
end

function interpInfo(
    lev :: AbstractArray{Float64, 1},
    target_lev :: Float64,
)

    target_i = -1
    for i=1:length(lev)-1
        if lev[i] <= target_lev <= lev[i+1]
            target_i = i
            break
        end
    end

    if target_i == -1
        throw(ErrorException("Cannot find target_lev: " * string(target_lev)))
    end

    wgt_of_i = (lev[i+1] - target_lev) / (lev[i+1] - lev[i])

    return target_i, wgt_of_i

end

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
 
        "--var"
            help = "Name of target variable"
            arg_type = String
            required = true
 
        "--lev"
            help = "Level coordinate varname."
            arg_type = String
            default = "lev"
 
        "--lat"
            help = "Latitude coordinate varname."
            arg_type = String
            default = "lat"
 
        "--lon"
            help = "Longitude coordinate varname."
            arg_type = String
            default = "lon"
 
        "--mask"
            help = "Mask varname."
            arg_type = String
            default = "mask"
 
        "--lev-unit"
            help = "Level coordinate unit."
            arg_type = String
            default = "m"
 
        "--lev-target"
            help = "Level interpolated (in meter, positive number)."
            arg_type = Float64
            default = 50.0

        "--time-range"
            help = "A range specifying the time. In format of [beg_time],[end_time]. Ex: 10 years of monthly data 1,120. `:` implies all data."
            arg_type = String
            default = ":"
     
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

lev_targets = parse(Float64, split(parsed["lev-target"], ","))
wgt = zeros(Float64, length(lev_targets), 2)

Dataset(parsed["input-file"], "r") do ds
    global lev, Δlev = adjustLev(ds[parsed["lev"]], parsed["lev-unit"])
    global Nx, Ny, Nt = size(ds[parsed["var"]])
end

if parsed["time-range"] == ":"
    time_range = 1:Nt
else
    time_range = parse(Int64, split(parsed["time-range"], ","))
    time_range = time_range[1]:time_range[2]
end

# Generate interpolation information
for i = 1:length(lev_targets)
    wgt[i, 1], wgt[i, 2] = interpInfo(lev, lev_targets[i])
end

result = zeros(Float64, Nx, Ny, length(lev_targets), Nt)
Dataset(parsed["input-file"], "r") do ds
   
    for l=1:length(lev_targets) 

        upper_layer_k   = wgt[l, 1]
        upper_layer_wgt = wgt[l, 2]

        lower_layer_k   = upper_layer_k + 1
        lower_layer_wgt = 1.0 - lower_layer_wgt

        _data = ds[parsed["var"]][time_range, upper_layer_k:lower_layer_k, :, :] |> nomissing
        for i=1:Nx, j=1:Ny, t=1:length(time_range)
            result[i, j, l, t] = _data[i, j, upper_layer_k, t] * upper_layer_wgt
                               + _data[i, j, lower_layer_k, t] * lower_layer_wgt
        end
    end

    Dataset(parsed["output-file"], "c") do ods

        defDim(ods, "Nx", Nx)
        defDim(ods, "Ny", Ny)
        defDim(ods, "target_lev", length(target_lev))
        defDim(ods, "time", Inf)

        for (varname, vardata, vardim, attrib) in [
            ("lat",  ods[parsed["lat"]][:]  |> nomissing, ("Nx", "Ny"), Dict()),
            ("lon",  ods[parsed["lon"]][:]  |> nomissing, ("Nx", "Ny"), Dict()),
            ("mask", ods[parsed["mask"]][:] |> nomissing, ("Nx", "Ny"), Dict()),
            ("target_lev", target_lev, ("target_lev",), Dict()),
            (parsed["var"], result, ("Nx", "Ny", "target_lev", length(time_range)), Dict()),
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
end



