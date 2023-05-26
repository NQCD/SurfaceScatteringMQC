using DrWatson
@quickactivate "SurfaceScatteringMQC"
using NQCDynamics
using SurfaceScatteringMQC
using Unitful, UnitfulAtomic
using LinearAlgebra: BLAS
BLAS.set_num_threads(1)

include(srcdir("outputs.jl"))
include(srcdir("generate_distribution.jl"))
include(srcdir("terminate.jl"))

function run_simulation(params)

    @unpack trajectories, dt, vibrational_number, incident_energy, M, temperature, method, model_id, Γ, bandwidth = params
    dt = dt*u"fs"
    temperature = temperature*u"K"
    bandwidth = austrip(bandwidth * u"eV")
    Γ = austrip(Γ * u"eV")

    Nmass = austrip(14u"u")
    Omass = austrip(16u"u")
    totalmass = Nmass + Omass
    μ = (Nmass * Omass) / totalmass
    atoms = Atoms([μ, totalmass])

    metalmodel = SurfaceScatteringMQC.get_model(model_id, M, Γ, bandwidth)
    adiabatic_model = SurfaceScatteringMQC.get_adiabatic_model(model_id, metalmodel)

    d = generate_distribution(vibrational_number, incident_energy, atoms.masses[1], atoms.masses[2], 10trajectories, adiabatic_model)

    if method === :DiabaticMDEF
        frictionmodel = AndersonHolstein(metalmodel.model, TrapezoidalRule(M, -bandwidth/2, bandwidth/2))
        sim = Simulation{eval(method)}(atoms, frictionmodel;
            friction_method=ClassicalMethods.WideBandExact((M+1)/bandwidth, 1/austrip(temperature)), temperature
        )
    elseif method === :Classical
        sim = Simulation{eval(method)}(atoms, metalmodel;
            temperature
        )
    elseif method === :AdiabaticIESH
        d *= FermiDiracState(0.0, temperature)
        @unpack decoherence = params
        if decoherence == :EDC
            sim = Simulation{eval(method)}(atoms, metalmodel;
                temperature, decoherence=SurfaceHoppingMethods.DecoherenceCorrectionEDC()
            )
        else
            sim = Simulation{eval(method)}(atoms, metalmodel;
                temperature,
            )
        end
    elseif method === :EhrenfestNA
        d *= FermiDiracState(0.0, temperature)
        sim = Simulation{eval(method)}(atoms, metalmodel;
            temperature
        )
    elseif method === :BCME
        d *= PureState(1, Diabatic())
        sim = Simulation{eval(method)}(atoms, metalmodel.model;
            temperature, bandwidth
        )
    end

    output = (OutputReactionOutcome, OutputVibrationalState(adiabatic_model, μ))

    result = run_dynamics(sim, (0.0, 1000u"fs"), d;
        output,
        selection=1:trajectories,
        dt,
        trajectories,
        callback=terminate,
        saveat=5u"fs",
        reduction=FileReduction(datadir(savename(params, "h5"))),
        abstol=1e-10,
        reltol=1e-10,
    )

    return result
end

all_params = Dict{String,Any}(
    "trajectories"=>[1],
    "dt"=>[0.25],
    "M"=>[200],
    "Γ"=>[
        @onlyif(("model_id"=="NOAu") && ("vibrational_number"==16), 0.15),
        1.5,
        @onlyif(("model_id"=="NOAu") && ("vibrational_number"==16), 15.0),
    ],
    "bandwidth"=>[100],
    "vibrational_number"=>[
        @onlyif(("model_id"=="NOAg") || ("model_id"=="NOAu"), 3),
        @onlyif(("model_id"=="NOAg") || ("model_id"=="NOAu"), 16),
    ],
    "incident_energy"=>collect(0.2:0.1:1.0),
    "temperature"=>[300],
    "method"=>[:Classical, :DiabaticMDEF, :AdiabaticIESH, :EhrenfestNA, :BCME],
    "model_id"=>["NOAu", "NOAg"],
    "decoherence"=>[
        @onlyif(("method"==:AdiabaticIESH) && (("incident_energy"==0.2) || ("incident_energy"==1.0)), :none),
        @onlyif("method"==:AdiabaticIESH, :EDC),
    ]
)

params = dict_list(all_params)[1] # Change index for other parameters

run_simulation(params)
