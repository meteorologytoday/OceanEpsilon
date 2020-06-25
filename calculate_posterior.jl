using Formatting
using NCDatasets
using ArgParse
using JSON
using Distributed

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-dir"
            help = "Input directory."
            arg_type = String
            required = true
 
        "--output-file"
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

parsed = parse_commandline()
print(json(parsed))

@everywhere parsed = parsed 
@everywhere include("worker_include.jl")

y_rng = 1:100

mkpath("output")

# doing Bayesian
N   = 201
ϵ   = range(0.5, 2.5, length=N) / 86400.0
Δϵ  = ϵ[2] - ϵ[1]

lat_flag = ( lat .>=  -20.0 ) .& ( lat .<= 20.0 )

σ_w = 1e-4  # 1cm/day tolerance

log_post = SharedArray{Float64}(N, length(workers()))
log_post .= 0.0

total_years = y_rng[end] - y_rng[1] + 1

# Distribute jobs
njob = Int64(total_years)
wkrs = workers()
nwkrs = length(wkrs)
ptr = 1

println("Initializing workers: ")
@sync for p in wkrs 
    @spawnat p let
        worker_init(
            data_file   = format("data/input.{:s}.nc", parse["model"]),
            domain_file = format("domains/domain.{:s}.nc", parse["model"])
            output_prefix = parse["model"],
            ϵ = ϵ,
        ) 
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

