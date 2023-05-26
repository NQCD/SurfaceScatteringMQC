
using StaticArrays: SMatrix, SVector
using LinearAlgebra: Hermitian
using DelimitedFiles: readdlm
using DrWatson: datadir


struct Morse{T}
    a::T
    x0::T
    D::T
end

function (morse::Morse)(x)
    (;a, x0, D) = morse
    return D*(exp(-2a*(x-x0)) - 2exp(-a*(x-x0)))
end

function ∂(morse::Morse, x)
    (;a, x0, D) = morse
    return D*(-2a*exp(-2a*(x-x0)) + 2a*exp(-a*(x-x0)))
end

MorseNO() = return Morse(1.48, 2.175, austrip(6.610u"eV"))
MorseNO⁻() = return Morse(1.20, 2.393, austrip(5.161u"eV"))

struct Repel{T}
    b::T
    x0::T
end

function (repel::Repel)(x)
    (;b, x0) = repel
    return exp(-b*(x-x0))
end

function ∂(repel::Repel, x)
    (;b, x0) = repel
    return -b*exp(-b*(x-x0))
end

struct TwoDimensionalPotential{X,Y,T}
    x_component::X
    y_component::Y
    shift::T
end

function (two_dimensional_potential::TwoDimensionalPotential)(x, y)
    (;x_component, y_component, shift) = two_dimensional_potential
    return x_component(x) + y_component(y) + shift
end

function ∂(two_dimensional_potential::TwoDimensionalPotential, x, y)
    (;x_component, y_component) = two_dimensional_potential
    return SVector{2}(∂(x_component, x), ∂(y_component, y))
end

struct SmoothStep{T}
    V̄ₖ::T
    x̃::T
    ã::T
end

function (smooth_step::SmoothStep)(x)
    (;V̄ₖ, x̃, ã) = smooth_step
    return V̄ₖ * (1 - tanh((x-x̃)/ã))
end

function ∂(smooth_step::SmoothStep, x)
    (;V̄ₖ, x̃, ã) = smooth_step
    return -V̄ₖ * sech((x-x̃)/ã)^2 / ã
end

struct PogoModel{T} <: NQCModels.DiabaticModels.DiabaticModel
    U0::TwoDimensionalPotential{Morse{T},Repel{T}}
    U1::TwoDimensionalPotential{Morse{T},Morse{T}}
    Vk::SmoothStep{T}
end

function PogoModel(modelid::String; Γ)
    parameters = readdlm(datadir("fitting", modelid, "parameters.txt"))
    return PogoModel(parameters; Γ)
end

function PogoModel(parameters::AbstractArray; Γ=0.0)
    U0 = TwoDimensionalPotential(
        MorseNO(),
        Repel(
            austrip(parameters[1] * u"Å^-1"),
            austrip(parameters[2] * u"Å")
        ),
        austrip(parameters[3] * u"eV")
    )

    U1 = TwoDimensionalPotential(
        Morse(
            austrip(parameters[4] * u"Å^-1"),
            austrip(parameters[5] * u"Å"),
            austrip(parameters[6] * u"eV")
        ),
        Morse(
            austrip(parameters[7] * u"Å^-1"),
            austrip(parameters[8] * u"Å"),
            austrip(parameters[9] * u"eV")
        ),
        austrip(parameters[10] * u"eV")
    )

    V̄ₖ = sqrt(austrip(Γ) / 2π)
    x̃ = austrip(0.0u"Å")
    ã = austrip(10.0u"Å")
    Vk = SmoothStep(V̄ₖ, x̃, ã)

    return PogoModel(U0, U1, Vk)
end

NQCModels.nstates(::PogoModel) = 2
NQCModels.ndofs(::PogoModel) = 1

function NQCModels.potential(model::PogoModel, r::Real)
    return NQCModels.potential(model, SVector{2}(r, 20.0))
end

function NQCModels.potential(model::PogoModel, r::AbstractVector)
    (;U0, U1, Vk) = model
    x, y = r
    V11 = U0(x, y)
    V22 = U1(x, y) 
    V12 = Vk(y)
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NQCModels.derivative!(model::PogoModel, D, r::AbstractVector)
    (;U0, U1, Vk) = model
    x, y = r
    ∂U0 = ∂(U0, x, y)
    ∂U1 = ∂(U1, x, y)
    ∂Vk = ∂(Vk, y)
    D[1] = Hermitian(SMatrix{2,2}(∂U0[1], 0.0, 0.0, ∂U1[1]))
    D[2] = Hermitian(SMatrix{2,2}(∂U0[2], ∂Vk, ∂Vk, ∂U1[2]))
    return D
end