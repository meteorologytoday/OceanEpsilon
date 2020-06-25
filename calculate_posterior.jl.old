@everywhere using NCDatasets
@everywhere using Statistics
@everywhere using Formatting
@everywhere using SharedArrays
@everywhere using Distributed

y_rng = 1:100

mkpath("output")

@everywhere const Re = 6.37122e6 
@everywhere const Ωe = 2π * 366 / (365 * 86400.0)
@everywhere const ρ  = 1e3

@everywhere Dataset("CESM_domains/domain.lnd.fv0.9x1.25_gx1v6.090309.nc", "r") do ds
    global mask = ds["mask"][:] |> nomissing
    global lat = ds["yc"][1, :] |> nomissing
    global lon = ds["xc"][:, 1] |> nomissing
end


@everywhere function loadData(rng)
    
    local τx, τy, wvel50

    #t = (365 * (y-1) + 1):(365*y)
    

    Dataset("TAUX.nc", "r") do ds
        τx = - (ds["TAUX"][:,:,rng] |> nomissing)    # TAUX is the stress force exerted to atmosphere
    end

    Dataset("TAUY.nc", "r") do ds
        τy = - (ds["TAUY"][:,:,rng] |> nomissing)
    end

    Dataset("WVEL_50m_rg.nc", "r") do ds
        wvel50 = ds["WVEL_50m"][:,:,rng] |> nomissing
    end

    return τx, τy, wvel50

end

@everywhere Nx = length(lon)
@everywhere Ny = length(lat)

println("Nx: ", Nx)
println("Ny: ", Ny)

@everywhere Δx = 2π * Re * cos.(deg2rad.(lat)) / Nx
@everywhere Δy =  π * Re / Ny
@everywhere Δa = Δx * Δy
@everywhere Δx_bnd = (Δx[1:end-1] + Δx[2:end]) / 2.0
@everywhere f  = 2Ωe * sin.(deg2rad.(lat))


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

@everywhere function work(
        d_rng;
        year::Integer,
        month::Integer,
    )

    w = zeros(Float64, Nx, Ny)

    println(format("[{:d}] Calculating {:03d}-{:02d}. ", myid(), year, month), d_rng)
    
    filename = format("./output/epsilon_log_posterior_{:03d}-{:02d}.nc", year, month)
    if isfile(filename)
        println("output already exists. Job skipped.")
        return
    end
 
    local τx, τy, w_obs, log_post

    τx, τy, w_obs = loadData(d_rng)
    log_post = zeros(Float64, N)

    for k=1:N, t=1:size(τx)[3]

        calW!(w, view(τx, :, :, t), view(τy,: ,:, t), ϵ[k])

        for i=1:Nx,j=1:Ny
            if (! lat_flag[j]) || mask[i, j] == 1
                continue
            end


            Δw = w[i, j] - w_obs[i, j, t]

            if isnan(Δw)

                continue
                #println(i, "; ", j)
                #throw(ErrorException("NaN generated..."))
            end

            log_post[k] += - ( Δw / σ_w )^2.0 / 2.0

        end
    end



    Dataset(filename, "c") do ds
        
        defDim(ds, "N", N)
        for (varname, vardata, vardim, attrib) in [
            ("log_post",  log_post, ("N",), Dict()),
            ("epsilon",  ϵ, ("N",), Dict()),
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


# doing Bayesian
@everywhere N   = 201
@everywhere ϵ   = range(0.5, 2.5, length=N) / 86400.0
@everywhere Δϵ  = ϵ[2] - ϵ[1]

@everywhere lat_flag = ( lat .>=  -20.0 ) .& ( lat .<= 20.0 )

@everywhere σ_w = 1e-4  # 1cm/day tolerance

log_post = SharedArray{Float64}(N, length(workers()))
log_post .= 0.0

# Years calculated
days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30 ,31]
beg_day_of_month = [ sum(days_of_month[1:i]) for i = 0:11 ] .+ 1

sum(days_of_month) != 365 && thorw(ErrorException("Sum is not right"))

total_years = y_rng[end] - y_rng[1] + 1

# Distribute jobs
njob = Int64(total_years * 12)
wkrs = workers()
nwkrs = length(wkrs)
ptr = 1

println("workers: ", wkrs)
@sync for job=1:njob

    global ptr

    println("Spawn job at worker: ", wkrs[ptr]) 

    @spawnat wkrs[ptr] let 
        y = floor(Int64, (job-1) / 12.0) + 1
        m = mod(job-1, 12) + 1

        beg_day = (y-1)*365 + beg_day_of_month[m]
        end_day = beg_day + days_of_month[m] - 1

        @elapsed work(
            beg_day:end_day;
            year=y,
            month=m,
        )
    end

    ptr = mod(ptr, nwkrs) + 1
end

println("calculation done.")


#=
println("make posterior")

log_post = sum(log_post, dims=(2,))[:,1]

println("Before substract max: ", log_post)
log_post .-= maximum(log_post)
println("Before exp: ", log_post)
posterior = exp.(log_post)
posterior /= sum(posterior) * Δϵ

println(posterior)

Dataset("epsilon_posterior.nc", "c") do ds
    
    defDim(ds, "N", N)
    for (varname, vardata, vardim, attrib) in [
        ("epsilon_posterior",  posterior, ("N",), Dict()),
        ("epsilon",  ϵ, ("N",), Dict()),
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




#=


 
Dataset("output.nc", "c") do ds

    defDim(ds, "Nx", Nx)
    defDim(ds, "Ny", Ny)
    defDim(ds, "time", Inf)


    for (varname, vardata, vardim, attrib) in [
        ("w",  w, ("Nx", "Ny"), Dict()),
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

=#



=#
