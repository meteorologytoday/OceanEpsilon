include("PolelikeCoordinate.jl")

using NCDatasets

const Re = 6.37122e6 
const Ωe = 2π * 366 / (365 * 86400.0)
const ρ  = 1e3

function loadDomain(filename)

    lat_c = ds["lat_c"][:] |> nomissing
    lon_c = ds["lon_c"][:] |> nomissing

    lat_v = ds["lat_v"][:] |> nomissing
    lon_v = ds["lon_v"][:] |> nomissing

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
    
    local τx, τy, wvel

    rng = ((year-1)*12):year*12

    Dataset(filename, "r") do ds
        τx  = - nomissing(ds["tauuo"][:, :,rng], NaN)
        τy  = - nomissing(ds["tauvo"][:, :,rng], NaN)
        wve = nomissing(ds["WVEL_50m"][:,:,rng], NaN)
    end

    return τx, τy, wvel

end


function calW!(
    gi :: PolelikeCoordinate.CurvilinearSphericalGridInfo,
    w  :: AbstractArray{Float64},
    τx :: AbstractArray{Float64},
    τy :: AbstractArray{Float64},
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
    DIV[:,  1] .= NaN
    DIV[:, Ny] .= NaN

    for i=1:Nx, j=2:Ny-1

        iw = (i==1)  ? Nx : i
        ie = (i==Nx) ? 1  : i

        DIV[i, j] = (
                  ( (u[ie, j] + u[i, j]) / 2 * gi.DY[i+1, j] - (u[iw, j] + u[i, j]) / 2.0 * gi.DY[i, j] )
                + ( (v[i,j+1] + v[i, j]) / 2 * gi.DX[i, j+1] - (v[i,j-1] + v[i, j]) / 2.0 * gi.DX[i, j] )
        ) / dσ[i, j]

    end

end

function worker_init(;
    data_file :: String,
    domain_file :: String,
    output_prefix :: String,
    ϵ           :: AbstractArray{Float64},
)
    global data_file = data_file
    global output_prefix = output_prefix
    global gi = loadDomain(domain_file)
    global w  = zeros(Float64, Nx, Ny)
    global ϵ = copy(ϵ)
    global Nϵ = length(ϵ)
end


function work(;
    year     :: Int64,
)
    global gi, w, ϵ
    
    println(format("[{:d}] Calculating year {:d}. ", myid(), year))
 
    local τx, τy, w_obs, log_post
    τx, τy, w_obs = loadData(data_file, year)

    log_post = zeros(Float64, Nϵ)
   
    for m = 1:12

        output_file = format("./output/{:s}_epsilon_log_posterior_{:03d}-{:02d}.nc", output_prefix, year, m)

        if isfile(output_file)
            println("output already exists. Job skipped.")
            return
        end
     

        for k = 1:Nϵ

            calW!(w, view(τx, :, :, m), view(τy,: ,:, m), ϵ[k])

                Δw = w[i, j] - w_obs[i, j, t]

                if isnan(Δw)

                    continue
                    #println(i, "; ", j)
                    #throw(ErrorException("NaN generated..."))
                end

                log_post[k] += - ( Δw / σ_w )^2.0 / 2.0

            end
        end


        Dataset(output_file, "c") do ds
            
            defDim(ds, "N_epsilon", Nϵ)
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


