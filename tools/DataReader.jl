module DataReader

    using NCDatasets
    using Formatting
    export getData  

#    include("nanop.jl")

    function getData(
        dir      :: String,
        varnames :: Union{String, Array},
        year_rng :: Union{Tuple, Array},
        idxes...;
        
        verbose=true,
        return_filenames=false,
    )
        
        beg_year, end_year = year_rng

        local only_one_variable = false

        if typeof(varnames) <: AbstractString
            only_one_variable = true
            varnames = [varnames]
        end

        if ! only_one_variable && length(idxes) == 1
            idxes = [ idxes[1] for _ in 1:length(varnames) ]  
        end

        if ( only_one_variable && length(idxes) > 1 ) || length(idxes) != length(varnames)
            println(varnames)
            println(idxes)
            throw(ErrorException("Number of idxes should match number of varnames or only 1 (applied to all varnames)."))
        end

        for k = 1:length(idxes)
            if length(idxes[k]) == 0
                throw(ErrorException("Spatial range not specified."))
            end
        end

        local data, new_idxes

        datas = Array{Any}(undef, length(varnames))
        new_idxes = Array{Any}(undef, length(varnames))
        flag_1s = [true for _ in 1:length(varnames)]  # flags to indicate reading first batch of data

        casename, filenames, file_year_rngs = statDirectory(dir)
        
        for y in beg_year:end_year
            index = searchIndex(y, file_year_rngs)

            filename = filenames[index]
            file_year_rng = file_year_rngs[index]

            filename = format("{:s}/{:s}", dir, filename)

            verbose && println("Loading file: ", filename)
            ds = Dataset(filename, "r")

            year_index_in_the_file = y - file_year_rng[1] + 1

            for (k, varname) in enumerate(varnames)

                rng = (idxes[k]..., ((year_index_in_the_file - 1) * 12 + 1):(year_index_in_the_file * 12) )

                partial_data = nomissing(ds[varname][rng...], NaN)

                if flag_1s[k]
                    new_size = size(partial_data)[1:end-1]
                    new_idxes[k]  = [Colon() for _ in 1:length(new_size)]
                    datas[k] = zeros(Float64, new_size..., (end_year - beg_year + 1) * 12)
                    datas[k] .= NaN
                    flag_1s[k] = false
                end

                y_relative_to_beg_year = y - beg_year
                datas[k][new_idxes[k]..., (y_relative_to_beg_year * 12 + 1):(y_relative_to_beg_year * 12 + 12)] = partial_data
              
            end

            close(ds)

        end 

        if return_filenames 
            return ((only_one_variable) ? datas[1] : datas), filenames
        else    
            return (only_one_variable) ? datas[1] : datas
        end
    end

    function searchIndex(
        y,
        year_rngs,
    )
        for (i, (first_year, last_year)) in enumerate(year_rngs)
            if first_year <= y <= last_year
                return i
            end
        end
    end


    function statDirectory(
        dir :: String
    )

        pattern_first = r"([A-Za-z0-9\-\_]+)_([0-9]{4})01-([0-9]{4})12\.nc"
        pattern_general = nothing
        casename = nothing
        
        beg_year = -1
        end_year = -1

        all_files = sort(readdir(dir))
        
        valid_filenames = []
        year_rngs = []
        
        for (i, filename) in enumerate(all_files)

            filename = basename(filename)

            if casename == nothing
                m = match(pattern_first, filename)
                if m == nothing
                    continue
                else
                    casename = m.captures[1]
                    pattern_general = Regex(format("{:s}_([0-9]{{4}})01-([0-9]{{4}})12\\.nc", casename))
                    first_year = parse(Int64, m.captures[2])
                    last_year = parse(Int64, m.captures[3])
                    beg_year = first_year
                    end_year = last_year
#                    println(format("First_year {}, last_year {}, years_per_file: {}", first_year, last_year, years_per_file))
                end
            end
            
            m = match(pattern_general, filename)
            
            if m == nothing
                
                continue
                
            else

                first_year = parse(Int64, m.captures[1])
                last_year = parse(Int64, m.captures[2])

                if first_year < beg_year
                    beg_year = first_year
                end

                if last_year > end_year
                    end_year = last_year
                end
                
                push!(valid_filenames, filename)
                push!(year_rngs, (first_year, last_year))
            end
        end
        
        if length(valid_filenames) == 0
            throw(ErrorException("Error: No valid file!"))
        end
         
        return casename,  valid_filenames, year_rngs
    end
end
