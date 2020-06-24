using Formatting
using NCDatasets
include("DataReader.jl")

using .DataReader


cmip_root = "/glade/collections/cmip/CMIP6/CMIP"

tauuo = DataReader.getData(
    format("{:s}/NCAR/CESM2/piControl/r1i1p1f1/Omon/tauuo/gn/latest/", cmip_root),
    "tauuo",
    (795, 805),
    (:, :);
)

Nx, Ny, Nt = size(tauuo)

Dataset("chop_test.nc", "c") do ods

    defDim(ods, "Nx", Nx)
    defDim(ods, "Ny", Ny)
    defDim(ods, "time", Inf)

    for (varname, vardata, vardim, attrib) in [
        ("tauuo", tauuo, ("Nx", "Ny", "time"), Dict()),
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
