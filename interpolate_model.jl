using NCDatasets
using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "Input file."
            arg_type = String
            required = true
 
        "--lev"
            help = "Level coordinate varname."
            arg_type = String
            default = "lev"
 
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


function adjustLev(
    lev :: AbstractArray{Float64, 1},
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
    

    return lev * factor
end


wgt = zeros(Float64, parsed

Dataset(parsed["input-file"], "r") do ds

    global lev = adjustLev(ds[parsed["lev"]], parsed["lev-unit"])
    
         

end


# load z coordinate


# judging z coordinate
# 1. Increasing with index? positive or negative?
# 2. Determine which level to interpolate
# 3.

 
