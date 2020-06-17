using Statistics
using NCDatasets

using ArgParse
using JSON

include("nanop.jl")
include("PCA.jl")

using .PCA

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--output-file"
            help = "Output file."
            arg_type = String
            required = true
 
        "--data-file"
            help = "PSLA data file."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

        "--sparsity"
            help = "Because solving PCA matrix is time consuming. This parameter let user to use coarse data to derive PCA. Sparseity `n` means to skip every `n` other grid point in lon / lat. So sparsity 0 (default) means do not skip any data point. Sparsity 1 means means data density will drop to 1/4 of the original (1/2 * 1/2)."
            arg_type = Int64
            default = 0


     
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

sparsity = parsed["sparsity"]
inc = sparsity + 1

Dataset(parsed["domain-file"], "r") do ds
    global lon = replace(ds["xc"][:, 1], missing=>NaN)
    global lat = replace(ds["yc"][1, :], missing=>NaN)

    global Nx = length(lon)
    global Ny = length(lat)

    global valid_grids = []
    global share_grids = []

    # detect the first lat
    global southest_j = 1
    for j = 1:Ny
        if 20.0 <= lat[j]
            global southest_j = j
            break
        end
    end


    for i = 1:inc:Nx, j = southest_j:inc:Ny
        push!(valid_grids, (i, j))
    end
end

Dataset(parsed["data-file"], "r") do ds
    
    global PSLA  = replace(ds["PSL_MA"][:, :, 1, :], missing=>NaN)
    global months  = ds.dim["time"]
    
    global EOF_input = zeros(Float64, length(valid_grids), months)
    for i = 1:length(valid_grids)
        EOF_input[i, :] = PSLA[valid_grids[i][1], valid_grids[i][2], :]
    end
end

modes = 2

println("Solving for PCA...")
eigen_vectors = PCA.findPCAs(EOF_input, num=modes)
println("done.")

PCAs = zeros(Float64, Nx, Ny, modes)
PCAs .= NaN

for i = 1:length(valid_grids)
    PCAs[valid_grids[i][1], valid_grids[i][2], :] = eigen_vectors[i, :]
end

# Map grid back because of sparsity
if sparsity > 0
    for m = 1:modes, i = 1:inc:Nx, j=southest_j:inc:Ny
        PCAs[i:min(Nx, i+sparsity), j:min(Ny, j+sparsity), m] .= PCAs[i, j, m]
    end
end



AO = zeros(Float64, months)

AO_mode = view(PCAs, :, :, 1)
for i = 1:length(AO)
    AO[i] = nansum(AO_mode .* view(PSLA, :, :, i))
end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "Nx", Nx)
    defDim(ds, "Ny", Ny)
    defDim(ds, "time", Inf)
    defDim(ds, "modes", modes)

    for (varname, vardata, vardim, attrib) in [
        ("PCAs",  PCAs, ("Nx", "Ny", "modes",), Dict()),
        ("AO",   AO, ("time",), Dict()),
    ]
        println("Doing varname:", varname)
        var = defVar(ds, varname, Float64, vardim)
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






