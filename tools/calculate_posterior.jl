using Formatting
using NCDatasets
using ArgParse
using JSON
using Distributed
using SharedArrays

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--years"
            help = "Years"
            arg_type = Int64
            required = true
 
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

@everywhere include(joinpath( @__DIR__ , "worker_include.jl"))

for p in procs()
    @spawnat p let

        global parsed = local_parsed 
        global y_rng = 1:parsed["years"]
        global N   = 201
        global ϵ   = range(0.5, 2.5, length=N) / 86400.0
        global Δϵ  = ϵ[2] - ϵ[1]
        global σ_w = 1e-4  # 1cm/day tolerance
        global total_years = y_rng[end] - y_rng[1] + 1

        global Nϵ = length(ϵ)

        global data_file = parsed["input-file"]
        global gi = loadDomain(parsed["domain-file"])
        global output_dir = parsed["output-dir"]
        global w  = zeros(Float64, gi.Nx, gi.Ny)

        _lon = mod.(gi.c_lon, 360.0) # need 160E ~ 90W
        global mask = (( -5 .<= gi.c_lat .<= 5 ) .& ( 160 .<= _lon .<= (360-90) ))

    end
end

mkpath(local_parsed["output-dir"])

# Distribute jobs
njob = Int64(total_years)


println("workers: ", wkrs)
ptr = 1
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

