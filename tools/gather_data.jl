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

        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--models"
            help = "Name of target variables, split by comma `,`."
            arg_type = String
            required = true
 
        "--vars"
            help = "Name of target variables, split by comma `,`."
            arg_type = String
            required = true
 
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
varnames   = String.(split(parsed["vars"], ","))
models     = String.(split(parsed["models"], ","))

Dataset(parsed["domain-file"], "r") do ds
    global Nx = ds.dim["Nx"]
    global Ny = ds.dim["Ny"]
end

Nt = 12*(time_range[2] - time_range[1] + 1)

data = Dict()

println("Output result... ")
for model in models

    output_file = format("data/input.{:s}.nc", model)

    Dataset(output_file, "c") do ds
        
        defDim(ds, "Nx", Nx)
        defDim(ds, "Ny", Ny)
        defDim(ds, "time", Inf)
        
        for varname in varnames
             
            var = defVar(ds, varname, Float64, ("Nx", "Ny", "time",))
            var.attrib["_FillValue"] = 1e20
             
            var[:, :, 1:Nt] = DataReader.getData(
                 format("cmip_data/{:s}/{:s}", model, varname),
                 varname,
                 time_range,
                 (:, :),
            )
            
        end
        
    end
    
end

