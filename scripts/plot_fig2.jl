using DrWatson
@quickactivate "SurfaceScatteringMQC"
using NQCModels
using CairoMakie
using JamesPlots
using Unitful, UnitfulAtomic
using ColorSchemes
using SurfaceScatteringMQC

colormap = reverse(ColorSchemes.linear_bmy_10_95_c71_n256)

function plot_single_model!(ax, modelid, energy_limits, energy_levels, gamma)

    bandwidth = austrip(100u"eV")
    gamma = austrip(gamma*u"eV")
    metalmodel = SurfaceScatteringMQC.get_model(modelid, 100, gamma, bandwidth)

    xgrid = range(0.9, 2.9, length=60)
    ygrid = range(0.5, 6, length=60)

    contourf!(ax, xgrid, ygrid, (x, y)->SurfaceScatteringMQC.ground_state_potential(metalmodel, x, y, (;x=1.15, z=6.0)),
        levels=energy_levels,
        colormap=ColorSchemes.get(colormap, energy_levels, energy_limits), extendlow=:auto
    )
    contour!(ax, xgrid, ygrid, (x, y)->SurfaceScatteringMQC.ground_state_potential(metalmodel, x, y, (;x=1.15, z=6.0)),
        levels=energy_levels,
        color=:black,
    )
end

function plot_both_models()
    fig = Figure(resolution=(JamesPlots.COLUMN_WIDTH, 1.0JamesPlots.RESOLUTION[2]), figure_padding=(3,2,1,1))
    axes = [MyAxis(fig[1,j], limits=(nothing, 2.9, 0.5, 6), xticks=1.0:0.5:2.5, xminorticks=1.25:0.5:2.75, ylabel="Height above the surface /Å") for i in 1:1, j in 1:2]

    hideydecorations!(axes[2]; ticks=false, minorticks=false)

    energy_limits = (0.0, 5.0)
    energy_levels = range(energy_limits...; length=6)
    plot_single_model!(axes[1,1], "NOAu", energy_limits, energy_levels, 1.5)
    plot_single_model!(axes[1,2], "NOAg", energy_limits, energy_levels, 1.5)

    Colorbar(fig[1,3]; limits=(-1.0, 5.0), colormap=cgrad(colormap, length(energy_levels), categorical=true),
        ticks=energy_levels, ticksize=19, tickalign=0.85, ticklabelpad=1, label="Potential energy /eV", flip_vertical_label=true, labelpadding=2
    )
    Label(fig[end+1,:], "Bond length /Å", padding=(0, 0, 0, 3))

    padding = (1,1,1,1)
    Label(fig[0,1], "NO/Au"; padding, tellwidth=false, tellheight=true)
    Label(fig[0,2], "NO/Ag"; padding, tellwidth=false, tellheight=true)

    return fig
end

save_figure(plotsdir("fig2_adiabatic_surfaces.png"), plot_both_models, fontsize=9)
