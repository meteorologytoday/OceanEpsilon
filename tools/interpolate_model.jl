using NCDatasets
using ArgParse
using JSON
using Formatting

include("DataReader.jl")
using .DataReader

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

    wgt_of_i = (lev[target_i+1] - target_lev) / (lev[target_i+1] - lev[target_i])

    return target_i, wgt_of_i

end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-dir"
            help = "Directory containing Input files."
            arg_type = String
            required = true
 
        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--var"
            help = "Name of target variable"
            arg_type = String
            required = true
 
        "--lev-target"
            help = "Level interpolated (in meter, positive number). If more than one target then use comma ',' to separate them"
            arg_type = String
            default = "50.0"

        "--time-range"
            help = "A range specifying the time (year). In format of [beg_time],[end_time]. Ex: 10 years of monthly data 1,10."
            arg_type = String
            required = true
     
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

time_range = parse.(Int64, split(parsed["time-range"], ","))

lev_targets = parse.(Float64, split(parsed["lev-target"], ","))
wgt_i = zeros(Int64, length(lev_targets))
wgt = zeros(Float64, length(lev_targets))




Dataset(parsed["domain-file"], "r") do ds
    global lev = nomissing( ds["lev_c"][:])
    global Δlev = lev[2:end] - lev[1:end-1]
    global Nx = ds.dim["Nx"]
    global Ny = ds.dim["Ny"]
end
    
Nt = 12*(time_range[2] - time_range[1] + 1)

# Generate interpolation information
for i = 1:length(lev_targets)
    wgt_i[i], wgt[i] = interpInfo(lev, lev_targets[i])
end

result = zeros(Float64, Nx, Ny, length(lev_targets), Nt)

   
println("Interpolating...")
for l=1:length(lev_targets) 

    upper_layer_k   = Int64(wgt_i[l])
    upper_layer_wgt = wgt[l]

    lower_layer_k   = upper_layer_k + 1
    lower_layer_wgt = 1.0 - upper_layer_wgt
    
    println(format("Doing level: {:f}. Between layer {:d} and {:d}. Upper layer wgt: {:f}. Lower layer wgt: {:f}.", lev_targets[l], upper_layer_k, lower_layer_k, upper_layer_wgt, lower_layer_wgt))

    _data = DataReader.getData(
        parsed["input-dir"],
        parsed["var"],
        time_range,
        (:, :, upper_layer_k:lower_layer_k);
    )

    for i=1:Nx, j=1:Ny, t=1:Nt
        result[i, j, l, t] = _data[i, j, 1, t] * upper_layer_wgt + _data[i, j, 2, t] * lower_layer_wgt
    end
end

println("Output result... ")

Dataset(parsed["output-file"], "c") do ods

    defDim(ods, "Nx", Nx)
    defDim(ods, "Ny", Ny)
    defDim(ods, "lev_targets", length(lev_targets))
    defDim(ods, "time", Inf)

    for (varname, vardata, vardim, attrib) in [
#        ("time",  ds["time"][time_range]  |> nomissing, ("time",), ds["time"].attrib),
#        ("lat",  ds[parsed["lat"]][:]  |> nomissing, ("Nx", "Ny"), Dict()),
#        ("lon",  ds[parsed["lon"]][:]  |> nomissing, ("Nx", "Ny"), Dict()),
        ("lev_targets", lev_targets, ("lev_targets",), Dict()),
        (parsed["var"], result, ("Nx", "Ny", "lev_targets", "time"), Dict()),
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



