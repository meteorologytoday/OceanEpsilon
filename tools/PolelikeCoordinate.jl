
if ! isdefined(Main, :GridFiles)
    include("GridFiles.jl")
end


module PolelikeCoordinate

    using ..GridFiles
    using LinearAlgebra: ⋅, normalize, normalize!, norm

    """

    This function calculate the angle between the real north and the
    north defined by vertices. I call it T-north.

    Each T-grid is defined by four vertices. As shown below:

    [4] --- [3]
     |       | 
     |       | 
    [1] --- [2]

    T-north for this grid is defined by the average of [1]->[4] and
    [2]->[3]. Likewise, T-east is the east defined by vertices, the
    average of [1]->[2], [4]->[3]

    Positive angle α means the local T-north and T-east will align
    with true north and east if rotating true north counterclockwise
    by angle α. 

    """

    d2r = π / 180.0
    r2d = 180 / π


    function getCyclicNeighbors(Nx, Ny, i, j)
        return mod(i-1 - 1, Nx) + 1, mod(i+1 - 1, Nx) + 1, mod(j-1 - 1, Ny) + 1, mod(j+1 - 1, Ny) + 1
    end


    function extend(
        a         :: AbstractArray{Float64, 2};
        cyclic_x  :: Bool,
        cyclic_y  :: Bool,
    )
        Nx, Ny = size(a)

        aa = Array{Float64}(undef, Nx + 2, Ny + 2)
        fill!(aa, NaN)

        aa[2:end-1, 2:end-1] = a
        
        if cyclic_x
            aa[1,   2:end-1] = a[end, :]
            aa[end, 2:end-1] = a[1,   :]
        end

        if cyclic_y
            aa[2:end-1,   1] = a[:, end]
            aa[2:end-1, end] = a[:,   1]
        end

        return aa 

    end

    abstract type GridInfo end

    struct CurvilinearSphericalGridInfo <: GridInfo

        R     :: Float64
        Ω     :: Float64

        Nx    :: Integer
        Ny    :: Integer

        c_lon :: AbstractArray{Float64, 2}
        c_lat :: AbstractArray{Float64, 2}
        c_f   :: AbstractArray{Float64, 2}

        α     :: AbstractArray{Float64, 2}
        cosα  :: AbstractArray{Float64, 2}
        sinα  :: AbstractArray{Float64, 2}
        dx_w  :: AbstractArray{Float64, 2}
        dx_c  :: AbstractArray{Float64, 2}
        dx_e  :: AbstractArray{Float64, 2}
        dy_s  :: AbstractArray{Float64, 2}
        dy_c  :: AbstractArray{Float64, 2}
        dy_n  :: AbstractArray{Float64, 2}
        
        ds1   :: AbstractArray{Float64, 2}
        ds2   :: AbstractArray{Float64, 2}
        ds3   :: AbstractArray{Float64, 2}
        ds4   :: AbstractArray{Float64, 2}

        dσ    :: AbstractArray{Float64, 2}

        weight_e :: AbstractArray{Float64, 2}    # ( Nx , Ny   )
        weight_n :: AbstractArray{Float64, 2}    # ( Nx   , Ny+1 )
        DX    :: AbstractArray{Float64, 2}       # ( Nx   , Ny+1 )   The X-side grid size.
        DY    :: AbstractArray{Float64, 2}       # ( Nx , Ny   )   The Y-side grid size.

      
        function CurvilinearSphericalGridInfo(;
            R      :: Float64,
            Ω      :: Float64,
            Nx     :: Integer,
            Ny     :: Integer,
            c_lon  :: AbstractArray{Float64, 2},  # center longitude
            c_lat  :: AbstractArray{Float64, 2},  # center latitude
            vs_lon :: AbstractArray{Float64, 3},  # vertices longitude (4, Nx, Ny)
            vs_lat :: AbstractArray{Float64, 3},  # vertices latitude  (4, Nx, Ny)
            angle_unit :: Symbol = :deg,
            sub_yrng :: Union{Colon, UnitRange} = Colon(),
        )
       
            α    = zeros(Float64, Nx, Ny)
            cosα = zeros(Float64, Nx, Ny)    
            sinα = zeros(Float64, Nx, Ny)   
            dx_w = zeros(Float64, Nx, Ny)
            dx_c = zeros(Float64, Nx, Ny)
            dx_e = zeros(Float64, Nx, Ny)
            dy_s = zeros(Float64, Nx, Ny)
            dy_c = zeros(Float64, Nx, Ny)
            dy_n = zeros(Float64, Nx, Ny)
            dσ = zeros(Float64, Nx, Ny)
            
            ds1 = zeros(Float64, Nx, Ny)
            ds2 = zeros(Float64, Nx, Ny)
            ds3 = zeros(Float64, Nx, Ny)
            ds4 = zeros(Float64, Nx, Ny)
        
            ps = zeros(Float64, 3, 4)

            true_north  = zeros(Float64, 3)
            true_east   = zeros(Float64, 3)
            true_upward = zeros(Float64, 3)


            c_lon_rad = copy(c_lon)
            c_lat_rad = copy(c_lat)

            vs_lon_rad = copy(vs_lon)
            vs_lat_rad = copy(vs_lat)
            
            weight_e = zeros(Float64, Nx, Ny)
            weight_n = zeros(Float64, Nx, Ny+1)

            DX = zeros(Float64, Nx, Ny+1)
            DY = zeros(Float64, Nx, Ny)


            if angle_unit == :deg

                c_lon_rad .*= d2r
                c_lat_rad .*= d2r

                vs_lon_rad .*= d2r 
                vs_lat_rad .*= d2r
     
            elseif angle_unit == :rad
                
                # do nothing

            else
                throw(ErrorException("Unknown `angle_unit`: " * angle_unit))
            end


            for i = 1:Nx, j = 1:Ny

                for k = 1:4

                    λ = vs_lon_rad[k, i, j]
                    θ = vs_lat_rad[k, i, j]

                    ps[1, k] = cos(θ) * cos(λ)
                    ps[2, k] = cos(θ) * sin(λ)
                    ps[3, k] = sin(θ)

                end

                ps .*= R

                u1 = ps[:, 2] - ps[:, 1]
                u2 = ps[:, 3] - ps[:, 2]
                u3 = ps[:, 3] - ps[:, 4]
                u4 = ps[:, 4] - ps[:, 1]

                ds1[i, j] = norm(u1)
                ds2[i, j] = norm(u2)
                ds3[i, j] = norm(u3)
                ds4[i, j] = norm(u4)

                dx_c[i, j] = (norm(u1) + norm(u3)) / 2.0
                dy_c[i, j] = (norm(u2) + norm(u4)) / 2.0
                dσ[i, j] = dx_c[i, j] * dy_c[i, j]

                grid_east  = u1+u3
                grid_north = u2+u4

                λc = c_lon_rad[i, j]
                θc = c_lat_rad[i, j]

                true_north[:]  = [ - sin(θc) * cos(λc), - sin(θc) * sin(λc),   cos(θc) ]
                true_east[:]   = [ - sin(λc)          ,   cos(λc)          ,   0.0     ]
                true_upward[:] = [   cos(θc) * cos(λc),   cos(θc) * sin(λc),   sin(θc) ]

                #grid_north -= (grid_north ⋅ true_upward) * true_upward

                cos_α = (grid_north ⋅ true_north) / norm(grid_north) #/ norm(true_north)
                cos_β = (grid_north ⋅ true_east)  / norm(grid_north)

                if abs(cos_α) > 1.0 && abs(cos_α - 1 < 1e-3)
                    cos_α = (cos_α > 0) ? 1.0 : -1.0
                end

                α[i, j] = acos(cos_α)

                if cos_β < 0.0   # Flip angle
                    α[i, j] = 2*π - α[i, j]
                end

                cosα[i, j] = cos(α[i, j])
                sinα[i, j] = sin(α[i, j])

            end
           

            c_f = 2Ω * sin.(c_lat_rad)
     
            for i = 1:Nx, j = 1:Ny

                i_w, i_e, j_s, j_n = getCyclicNeighbors(Nx, Ny, i, j)

                dx_w[i, j] = (dx_c[i_w, j] + dx_c[i,   j]) / 2.0
                dx_e[i, j] = (dx_c[i,   j] + dx_c[i_e, j]) / 2.0

                dy_s[i, j] = (dy_c[i, j_s] + dy_c[i, j  ]) / 2.0
                dy_n[i, j] = (dy_c[i, j  ] + dy_c[i, j_n]) / 2.0

            end

            #
            # weight_e[i, ?] is the relative portion of grid weighting
            # to the west of the boundary of eastward vectors
            #
            #                     bnd
            #                     [i]
            #     |                |                  |
            #     |    grid[i-1]   |    grid[i]       |
            #     |                |                  |
            #     | <-- dx[i-1]--> | <--- dx[i] --->  |
            #     |                |                  |
            #
            #   weight_e[i, ?] = dx[i-1] / (dx[i-1] + dx[i])
            #
            #   weight_n would be the same idea with 
            #   portion represent the south grid
            #

            for j = 1:Ny  # i=1
                weight_e[1, j] = dx_c[Nx, j] / (dx_c[Nx, j] + dx_c[1, j])
            end
            for i = 2:Nx, j = 1:Ny
                weight_e[i, j] = dx_c[i-1, j] / (dx_c[i-1, j] + dx_c[i, j])
            end

            # Ignore the northest and the southest because information
            # is unknown
            for i = 1:Nx, j = 2:Ny
                weight_n[i, j] = dy_c[i, j-1] / ( dy_c[i, j-1] + dy_c[i, j] )
            end


            # Calculate DX
            for i = 1:Nx, j = 1:Ny
                DX[i, j] = ds1[i, j]
            end
            DX[:, Ny+1] = ds3[:, Ny]

            # Calculate DY
            for i = 1:Nx, j = 1:Ny
                DY[i, j] = ds4[i, j]
            end

            if sub_yrng == Colon()
                sub_yrng = 1:Ny
            end

            new_Ny = length(sub_yrng)
            sub_yrng_ext = sub_yrng[1]:sub_yrng[end]+1
           
            return new(
                R,
                Ω,
                Nx,
                new_Ny,
                c_lon_rad[:, sub_yrng],
                c_lat_rad[:, sub_yrng],
                c_f[:, sub_yrng],
                α[:, sub_yrng],
                cosα[:, sub_yrng],
                sinα[:, sub_yrng],
                dx_w[:, sub_yrng],
                dx_c[:, sub_yrng],
                dx_e[:, sub_yrng],
                dy_s[:, sub_yrng],
                dy_c[:, sub_yrng],
                dy_n[:, sub_yrng],
                ds1[:, sub_yrng],
                ds2[:, sub_yrng],
                ds3[:, sub_yrng],
                ds4[:, sub_yrng],
                dσ[:, sub_yrng],
                weight_e[:, sub_yrng],
                weight_n[:, sub_yrng_ext],
                DX[:, sub_yrng_ext],
                DY[:, sub_yrng],
            )
     
        end
    end


    function projectSpherical!(
        gi    :: CurvilinearSphericalGridInfo,
        ivf_e :: AbstractArray{Float64, 2},     # input vector field east
        ivf_n :: AbstractArray{Float64, 2},     # input vector field north
        ovf_e :: AbstractArray{Float64, 2},    # output vector field east
        ovf_n :: AbstractArray{Float64, 2};    # output vector field north
        direction = :Forward,
    )

        if direction == :Forward   # from outside world onto dispalced pole grid

            for i=1:gi.Nx, j=1:gi.Ny
                ovf_e[i, j] =   ivf_e[i, j] * gi.cosα[i, j] - ivf_n[i, j] * gi.sinα[i, j]
                ovf_n[i, j] =   ivf_e[i, j] * gi.sinα[i, j] + ivf_n[i, j] * gi.cosα[i, j]
            end
     
        else                       # from displaced pole grid onto outside world

            for i=1:gi.Nx, j=1:gi.Ny
                ovf_e[i, j] =   ivf_e[i, j] * gi.cosα[i, j] + ivf_n[i, j] * gi.sinα[i, j]
                ovf_n[i, j] = - ivf_e[i, j] * gi.sinα[i, j] + ivf_n[i, j] * gi.cosα[i, j]
            end

        end

    end

    mutable struct CylindricalGridInfo <: GridInfo

        R     :: Float64
        Ω     :: Float64
        Ly    :: Float64
        β      :: Float64

        Nx    :: Integer
        Ny    :: Integer

        c_lat :: AbstractArray{Float64, 2}
        c_lon :: AbstractArray{Float64, 2}
        c_y   :: AbstractArray{Float64, 2}
        c_f   :: AbstractArray{Float64, 2}

        dx_w  :: AbstractArray{Float64, 2}
        dx_c  :: AbstractArray{Float64, 2}
        dx_e  :: AbstractArray{Float64, 2}
        dy_s  :: AbstractArray{Float64, 2}
        dy_c  :: AbstractArray{Float64, 2}
        dy_n  :: AbstractArray{Float64, 2}
        
        dσ    :: AbstractArray{Float64, 2}

        weight_e :: AbstractArray{Float64, 2}    # ( Nx , Ny   )
        weight_n :: AbstractArray{Float64, 2}    # ( Nx   , Ny+1 )
        DX    :: AbstractArray{Float64, 2}       # ( Nx   , Ny+1 )   The X-side grid size.
        DY    :: AbstractArray{Float64, 2}       # ( Nx , Ny   )   The Y-side grid size.
      
        function CylindricalGridInfo(;
            R       :: Float64,
            Ω       :: Float64,
            Nx      :: Integer,
            Ny      :: Integer,
            Ly      :: Float64,
            lat0    :: Float64 = 0,  # in rad
            β       :: Float64 = 0,  
            sub_yrng :: Union{Colon, UnitRange} = Colon(),
        )

            Δx = 2π * R / Nx
            Δy = Ly     / Ny

            Δλ = 2π / Nx
            c_lon = collect(Float64, (1:Nx) .- 0.5) * Δλ

            if mod(Ny, 2) == 0
                c_y = collect(Float64, (-Ny/2+0.5:1:Ny/2-0.5)) * Δy
            else
                o = (Ny-1) / 2
                c_y = collect(Float64, (-o:1:o)) * Δy
            end

    #        if lat0 + Ly/2*β >= 90.0 || lat0 - Ly/2*β <= -90.0
    #            throw(ErrorException("Bad choice of Ly, lat0 and β"))
    #        end



            c_lat = zeros(Float64, Nx, Ny)  # meaningless
            c_lat .= lat0

            c_lon = repeat(reshape(c_lon, :, 1), outer=(1, Ny))
            c_y   = repeat(reshape(c_y, 1, :), outer=(Nx, 1))

            c_f   = 2Ω * sin(lat0) .+ β * c_y

            dx_w = zeros(Float64, Nx, Ny)
            dx_w .= Δx
            
            dx_c = copy(dx_w)
            dx_e = copy(dx_w)

            dy_s = zeros(Float64, Nx, Ny)
            dy_s .= Δy

            dy_c = copy(dy_s)
            dy_n = copy(dy_s)
            
            dσ = dx_c .* dy_c
            
            weight_e = zeros(Float64, Nx, Ny)
            weight_n = zeros(Float64, Nx, Ny+1)

            weight_e .= .5
            weight_n .= .5

            DX = zeros(Float64, Nx, Ny+1)
            DY = zeros(Float64, Nx, Ny)

            DX .= Δx
            DY .= Δy

            if sub_yrng == Colon()
                sub_yrng = 1:Ny
            end

            new_Ny = length(sub_yrng)
            sub_yrng_ext = sub_yrng[1]:sub_yrng[end]+1
           
            return new(
                R,
                Ω,
                Ly,
                β,
                Nx,
                new_Ny,
                c_lat[:, sub_yrng],
                c_lon[:, sub_yrng],
                c_y[:, sub_yrng],
                c_f[:, sub_yrng],
                dx_w[:, sub_yrng],
                dx_c[:, sub_yrng],
                dx_e[:, sub_yrng],
                dy_s[:, sub_yrng],
                dy_c[:, sub_yrng],
                dy_n[:, sub_yrng],
                dσ[:, sub_yrng],
                weight_e[:, sub_yrng],
                weight_n[:, sub_yrng_ext],
                DX[:, sub_yrng_ext],
                DY[:, sub_yrng],
            )
     
        end
    end

    function genGridInfo(
        gf       :: GridFiles.GridFile;
        sub_yrng :: Union{Colon, UnitRange} = Colon(),
    )

        local gi

        if typeof(gf) <: GridFiles.CurvilinearSphericalGridFile
            
            gi = CurvilinearSphericalGridInfo(;
                R=gf.R,
                Ω=gf.Ω,
                Nx=gf.Nx,
                Ny=gf.Ny,
                c_lon=gf.xc,
                c_lat=gf.yc,
                vs_lon=gf.xv,
                vs_lat=gf.yv,
                area=gf.area,
                angle_unit=:deg,
                sub_yrng=sub_yrng,
            )
           

        elseif typeof(gf) <: GridFiles.CylindricalGridFile
 
            gi = CylindricalGridInfo(;
                R=gf.R,
                Ω=gf.Ω,
                Nx=gf.Nx,
                Ny=gf.Ny,
                Ly=gf.Ly,
                lat0=gf.lat0,
                β=gf.β,
                sub_yrng=sub_yrng,
            )
            
        end


        return gi
    end

end
