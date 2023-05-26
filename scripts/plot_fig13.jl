using DrWatson
@quickactivate "SurfaceScatteringMQC"
using NQCDynamics
using SurfaceScatteringMQC
using Unitful, UnitfulAtomic
using CairoMakie
using JamesPlots

include(srcdir("plotting_variables.jl"))

M = 200
bandwidth = austrip(100u"eV")
Γ = austrip(1.5u"eV")

metalmodel = SurfaceScatteringMQC.get_model("NOAu", M, Γ, bandwidth)

Nmass = austrip(14u"u")
Omass = austrip(16u"u")
totalmass = Nmass + Omass
μ = (Nmass * Omass) / totalmass
m = μ
model = metalmodel.model.U0.x_component
(;D, a) = model
ω = sqrt(2D * a^2 / m)
Ek = austrip(1.0u"eV")

function gamma_function(gamma, x)
    metalmodel = SurfaceScatteringMQC.get_model("NOAu", M, gamma, bandwidth)
    Vkfunc = metalmodel.model.Vk
    Γ = 2π * Vkfunc(x)^2
    return ustrip(auconvert(u"eV", Γ))
end

function nuclear_frequency(Ek, n)
    return ustrip(auconvert(u"eV", 1 / ( (1/Ek) + (1/(n*ω)) )))
end

Γ = austrip(1.5u"eV")

function shade_region!(ax, Ek, gamma)
    heights = [0.0, 3.5]
    gamma_limits = gamma_function.(gamma, austrip.(heights .* u"Å"))
    band!(ax, Ek, gamma_limits..., color=:gray80)
    hlines!(ax, gamma_limits, linestyle=dash, color=:black)
end

function plot_fig()
    fig = Figure(figure_padding=(2, 1, 2, 1))
    ax = MyAxis(fig[1,1], ylabel="Energy /eV", xlabel="Incidence energy /eV",
        limits=(0.15, 1.05, nothing, nothing), yscale=log10)

    Ek = range(0.15, 1.05, length=100)
    Ek_hartree = austrip.(Ek .* u"eV")

    shade_region!(ax, Ek, austrip(1.5u"eV"))
    shade_region!(ax, Ek, austrip(15u"eV"))
    shade_region!(ax, Ek, austrip(0.15u"eV"))

    lines!(ax, Ek, nuclear_frequency.(Ek_hartree, 3), label="νᵢ = 3")
    lines!(ax, Ek, nuclear_frequency.(Ek_hartree, 16), label="νᵢ = 16")

    Label(fig[1,1], "Γ = 15.0 eV"; halign=:left, valign=1.0, tellwidth=false, tellheight=false, padding=(10, 10, 10, 10))
    Label(fig[1,1], "Γ = 1.50 eV"; halign=:left, valign=0.49, tellwidth=false, tellheight=false, padding=(10, 10, 10, 10))
    Label(fig[1,1], "Γ = 0.15 eV"; halign=:left, valign=-0.03, tellwidth=false, tellheight=false, padding=(10, 10, 10, 10))

    Legend(fig[1,1], ax; valign=:top, halign=:right, tellwidth=false, tellheight=false, margin=(10, 10, 10, 10))

    return fig
end

save_figure(plotsdir("fig13_regimes.png"), plot_fig, fontsize=9)
