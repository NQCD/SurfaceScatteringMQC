
using DrWatson
using DataFrames
using HDF5

get_check(::Missing) = x->ismissing.(x)
get_check(v::Any) = x->x .== v

function select_data(results; kwargs...)
    df = copy(results)
    for (k, v) in kwargs
        subset!(df, k => get_check(v); skipmissing=true)
    end
    if nrow(df) == 1
        return df[1,:]
    else
        display(kwargs)
        throw(error("Data incorrect? Found $(nrow(df)) rows instead of 1."))
    end
end

function read_data_file(directory, params, quantity, dtype::Type{T})::Vector{T} where {T}
    output = Dict{String,Any}()
    filename = savename(params, "h5")

    path = joinpath(directory, filename)
    if !isfile(path)
        @error "Cannot find" path
        error("File does not exist!")
    end

    fid = h5open(joinpath(directory, filename), "r")
    if !haskey(first(fid), quantity)
        close(fid)
        error("Quantity does not exist!")
    end

    output = dtype[]
    sizehint!(output, params["trajectories"])
    for traj in keys(fid)
        dst = open_dataset(fid, joinpath(traj, quantity))
        if dtype <: Integer
            data = HDF5.read(dst, dtype)
        else
            data = HDF5.read(dst)
        end
        push!(output, data)
        close(dst)
    end

    close(fid)

    return output
end

select_data_entry(results; kwargs...) = results[savename(kwargs)]
