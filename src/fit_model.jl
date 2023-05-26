using Optim
using NQCModels
using LinearAlgebra: eigvals, LinearAlgebra, diagind
using Unitful, UnitfulAtomic
using CSV, DataFrames

struct FittingPoint{T}
    x::T
    z::T
    target_energy::T
end

function load_fitting_data(directory)
    files = readdir(directory)
    fitting_points = FittingPoint{Float64}[]
    for f in files
        df = CSV.read(joinpath(directory, f), DataFrame; header=false)
        energy = parse(Float64, split(f, "eV.csv"; keepempty=false)[1])
        for row in eachrow(df)
            push!(fitting_points, FittingPoint(row[1], row[2], energy))
        end
    end
    return fitting_points
end

function calculate_model_energy(model, x, z)
    x = austrip(x * u"Å")
    z = austrip(z * u"Å")

    position = [x z]
    eigs = eigvals(potential(model, position))
    a = sum(eigs[i] for i in 1:div(length(eigs),2))
    b = NQCModels.state_independent_potential(model, position)
    return ustrip(auconvert(u"eV",  a + b))
end

function ground_state_potential(model, x, z, reference_point::FittingPoint)
    shift = calculate_model_energy(model, reference_point.x, reference_point.z)
    return calculate_model_energy(model, x, z) - shift + reference_point.target_energy
end

function ground_state_potential(model, x, z, reference_point)
    shift = calculate_model_energy(model, reference_point.x, reference_point.z)
    return calculate_model_energy(model, x, z) - shift
end
