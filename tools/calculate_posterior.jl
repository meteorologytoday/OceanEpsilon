using Formatting
using NCDatasets
using ArgParse
using JSON
using Distributed
using SharedArrays

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Input directory."
            arg_type = String
            required = true
 
        "--output-dir"
            help = "Output file of posterior."
            arg_type = String
            required = true

        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)
end

local_parsed = parse_commandline()
print(json(local_parsed, 4))

wkrs = workers()
nwkrs = length(wkrs)

for p in wkrs
    @spawnat p let
        global parsed = local_parsed 
    end
end

@everywhere include(joinpath( @__DIR__ , "worker_include.jl"))

y_rng = 1:5

mkpath(parsed["output-dir"])

# doing Bayesian
N   = 201
ϵ   = range(0.5, 2.5, length=N) / 86400.0
Δϵ  = ϵ[2] - ϵ[1]

#lat_flag = ( lat .>=  -20.0 ) .& ( lat .<= 20.0 )

σ_w = 1e-4  # 1cm/day tolerance

log_post = SharedArray{Float64}(N, length(workers()))
log_post .= 0.0

total_years = y_rng[end] - y_rng[1] + 1

# Distribute jobs
njob = Int64(total_years)
ptr = 1

println("Initializing workers: ")
for p in wkrs
    @spawnat p let
        global data_file = parsed["input-file"]
        global gi = loadDomain(parsed["domain-file"])
        global output_dir = parsed["output-dir"]
        global w  = zeros(Float64, gi.Nx, gi.Ny)
        global ϵ = copy(ϵ)
        global Nϵ = length(ϵ)

        global mask = ( -10 .<= gi.c_lat .<= 10 )
    end
end



println("workers: ", wkrs)
@sync for job=1:njob

    global ptr

    println("Spawn job at worker: ", wkrs[ptr]) 
    
    y = y_rng[1] + job - 1
    @spawnat wkrs[ptr] let 
        @elapsed work(
            year=y,
        )
    end

    ptr = mod(ptr, nwkrs) + 1
end

println("calculation done.")

