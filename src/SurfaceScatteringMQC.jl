module SurfaceScatteringMQC

    include("fit_model.jl")
    include("pogo_model.jl")

    function get_model(modelid, M, Γ, bandwidth)
        AndersonHolstein(
            PogoModel(modelid; Γ),
            ShenviGaussLegendre(M, -bandwidth/2, bandwidth/2)
        )
    end

    get_adiabatic_model(metalmodel) = AdiabaticStateSelector(metalmodel.model, 1)

end # module
