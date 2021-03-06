include("PolelikeCoordinate.jl")

using Formatting
using NCDatasets
using JSON
using Distributed
using SharedArrays

const Re = 6.37122e6 
const Ωe = 2π * 366 / (365 * 86400.0)
const ρ  = 1e3

function loadDomain(filename)

    println("Loading domain: ", filename)

    ds = Dataset(filename, "r")

    lat_c = ds["lat_c"][:] |> nomissing
    lon_c = ds["lon_c"][:] |> nomissing

    lat_v = ds["lat_v"][:] |> nomissing
    lon_v = ds["lon_v"][:] |> nomissing

    close(ds)

    Nx, Ny = size(lat_c)

    return PolelikeCoordinate.CurvilinearSphericalGridInfo(;
        R = Re,
        Ω = Ωe,
        Nx = Nx,
        Ny = Ny,
        c_lon = lon_c,
        c_lat = lat_c,
        vs_lon = lon_v,
        vs_lat = lat_v,
        angle_unit=:deg,
    )

end


function loadData(filename, year)
    
    local τx, τy, wo

    rng = ((year-1)*12+1):year*12

    Dataset(filename, "r") do ds
        τx  = nomissing(ds["tauuo"][:, :,rng], NaN)
        τy  = nomissing(ds["tauvo"][:, :,rng], NaN)
        wo  = nomissing(ds["wo"][:,:,1,rng], NaN)
    end

    return τx, τy, wo

end


function calW!(
    gi :: PolelikeCoordinate.CurvilinearSphericalGridInfo,
    w  :: AbstractArray,
    τx :: AbstractArray,
    τy :: AbstractArray,
    ϵ  :: Float64,
)

    Nx = gi.Nx
    Ny = gi.Ny

    Mx  = zeros(Float64, Nx, Ny)
    My  = zeros(Float64, Nx, Ny)
    DIV = zeros(Float64, Nx, Ny)

    Mx .= NaN
    My .= NaN

    for i=1:Nx, j=1:Ny
        M̃ = (τx[i, j] + τy[i, j] * im) / (ϵ + gi.c_f[i, j] * im)
        Mx[i, j] = real(M̃) 
        My[i, j] = imag(M̃) 
    end

    calDIV!(gi, DIV, Mx, My)

    w .= DIV / ρ 

end

function calDIV!(
    gi  :: PolelikeCoordinate.CurvilinearSphericalGridInfo,
    DIV :: AbstractArray{Float64, 2},
    u :: AbstractArray{Float64, 2},
    v :: AbstractArray{Float64, 2},
)
    Nx = gi.Nx
    Ny = gi.Ny

    DIV[:,  1] .= NaN
    DIV[:, Ny] .= NaN

    for i=1:Nx, j=2:Ny-1

        iw = (i==1)  ? Nx : i
        ie = (i==Nx) ? 1  : i

        DIV[i, j] = (
                  ( (u[ie, j] + u[i, j]) / 2 * gi.DY[ie, j] - (u[iw, j] + u[i, j]) / 2.0 * gi.DY[i, j] )
                + ( (v[i,j+1] + v[i, j]) / 2 * gi.DX[i, j+1] - (v[i,j-1] + v[i, j]) / 2.0 * gi.DX[i, j] )
        ) / gi.dσ[i, j]

    end

end


function work(;
    year     :: Int64,
)
    global gi, w, ϵ, Nϵ, mask
    

 
    local τx, τy, w_obs, log_post
    τx, τy, w_obs = loadData(data_file, year)

    log_post = zeros(Float64, Nϵ)
   
    for m = 1:12
        println(format("[{:d}] Calculating {:d}-{:02d}. ", myid(), year, m))

        output_file = format("{:s}/epsilon_log_posterior_{:03d}-{:02d}.nc", output_dir, year, m)

        if isfile(output_file)
            println("[", output_file,"] already exists. Job skipped.")
            return
        end
     

        for k = 1:Nϵ

            calW!(gi, w, view(τx, :, :, m), view(τy,: ,:, m), ϵ[k])

            for i=1:gi.Nx, j=1:gi.Ny

                Δw = w[i, j] - w_obs[i, j, m]

                if isnan(Δw) || mask[i, j] == 0

                    continue
                    #println(i, "; ", j)
                    #throw(ErrorException("NaN generated..."))
                end

                log_post[k] += - ( Δw / σ_w )^2.0 / 2.0

            end

        end


        Dataset(output_file, "c") do ds
            
            defDim(ds, "N_eps", Nϵ)
            for (varname, vardata, vardim, attrib) in [
                ("log_post",  log_post, ("N_eps",), Dict()),
                ("eps",  ϵ, ("N_eps",), Dict()),
            ]

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

    end
end


