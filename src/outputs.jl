
function OutputReactionOutcome(sol, i)::Int
    final_position = DynamicsUtils.get_positions(last(sol.u))
    x = final_position[1]
    z = final_position[2]
    z > austrip(5u"Å") && return 1
    x > austrip(2.25u"Å") && return -1

    return 0
end

struct OutputVibrationalState{M}
    adiabatic_model::M
    μ::Float64
end

function (output::OutputVibrationalState)(sol, i)::Int
    ufinal = last(sol.u)
    positions = DynamicsUtils.get_positions(ufinal)
    z = positions[2]

    bond_lengths = 0.5:0.01:10.0

    if z >= austrip(5u"Å")
        r = positions[1]
        v = DynamicsUtils.get_velocities(ufinal)[1]
        try
            return InitialConditions.QuantisedDiatomic.quantise_1D_vibration(
                output.adiabatic_model, output.μ, r, v;
                bond_lengths
            )
        catch e
            if e isa ArgumentError
                return -100
            else
                rethrow(e)
            end
        end
    else
        return -10
    end
end
