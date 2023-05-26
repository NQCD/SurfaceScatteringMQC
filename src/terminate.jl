using NQCDynamics: DynamicsUtils
using Unitful, UnitfulAtomic

function termination_condition(u, t, integrator)::Bool
    x = DynamicsUtils.get_positions(u)[1]
    z = DynamicsUtils.get_positions(u)[2]

    x > austrip(2.25u"Å") && return true
    z > austrip(5u"Å") && return true

    return false
end
terminate = DynamicsUtils.TerminatingCallback(termination_condition)
