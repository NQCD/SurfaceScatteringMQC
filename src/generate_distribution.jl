using Unitful, UnitfulAtomic
using NQCDynamics: InitialConditions, DynamicalDistribution
using Distributions: Normal

function generate_distribution(vibrational_state, incident_energy, μ, total_mass, trajectories, adiabatic_model)

    bond_lengths = 1.5:0.01:7.0

    bond_lengths, bond_velocities = InitialConditions.QuantisedDiatomic.generate_1D_vibrations(
        adiabatic_model, μ, vibrational_state;
        samples=trajectories, bond_lengths
    )

    incident_energy = austrip(incident_energy * u"eV")
    vertical_velocity = sqrt(2incident_energy / total_mass)
    height_velocity_distribution = Normal(-vertical_velocity, 0.0)
    height_position_distribution = Normal(austrip(5u"Å"), 0.0)

    velocity = [[v rand(height_velocity_distribution)] for v in bond_velocities]
    position = [[r rand(height_position_distribution)] for r in bond_lengths]

    d = DynamicalDistribution(velocity, position, (1,2))
    return d
end
