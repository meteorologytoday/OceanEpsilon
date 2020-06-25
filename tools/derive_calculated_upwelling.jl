using Formatting
using NCDatasets
include("DataReader.jl")

using .DataReader


cmip_root = "/glade/collections/cmip/CMIP6/CMIP"

tauuo, filenames = DataReader.getData(
    format("{:s}/NCAR/CESM2/piControl/r1i1p1f1/Omon/tauuo/gn/latest/", cmip_root),
    "tauuo",
    (795, 805),
    (:, :);
    return_filenames=true
)


Nx, Ny, Nt = size(tauuo)

@everywhere function calW!(w, τx, τy, ϵ)

    Mx  = zeros(Float64, Nx, Ny)
    My  = zeros(Float64, Nx, Ny)
    DIV = zeros(Float64, Nx, Ny)

    Mx .= NaN
    My .= NaN

    for i=1:Nx, j=1:Ny
        if mask[i, j] == 1
            continue
        end

        M̃ = (τx[i, j] + τy[i, j] * im) / (ϵ + f[j] * im)
        Mx[i, j] = real(M̃)
        My[i, j] = imag(M̃)
    end

    calDIV!(DIV, Mx, My)

    w .= DIV / ρ

end

@everywhere function calDIV!(
    DIV :: AbstractArray{Float64, 2},
    u :: AbstractArray{Float64, 2},
    v :: AbstractArray{Float64, 2},
)
    DIV[:, 1]  .= NaN
    DIV[:, Ny] .= NaN

    #println("Size of DIV", size(DIV))
    #println("Size of v", size(u))
    #println("Size of u", size(v))

    for i=1:Nx, j=2:Ny-1

        iw = (i==1)  ? Nx : i
        ie = (i==Nx) ? 1  : i

        DIV[i, j] = ( (u[ie, j] - u[iw, j]) * Δy
                + (v[i,j+1] + v[i,j])/2.0 * Δx_bnd[j] - (v[i,j-1] + v[i,j])/2.0 * Δx_bnd[j-1]
        ) / Δa[j]
    end

end







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
