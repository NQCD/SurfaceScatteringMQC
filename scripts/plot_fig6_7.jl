
using DrWatson
@quickactivate "ElbowDynamics"
using DataFrames
using CairoMakie
using JamesPlots
using Unitful, UnitfulAtomic
using CSV, DataFrames
using ColorSchemes

include(srcdir("plotting_variables.jl"))

M = 200
bandwidth = 100
gamma = 1.5
dt = 0.25

function plot_entry!(ax, input, method)
    df = CSV.read(input, DataFrame)

    map1 = range(colorant"white", color_dict[method]; length=6)[end-2:end-1]
    map2 = range(color_dict[method], colorant"black"; length=10)[begin:begin+1]
    colormap = ColorScheme(vcat(map1, map2))

    energies = 0.2:0.1:1.0
    for i=0:3
        selection = Symbol(:state, i)
        y = getproperty(df, selection)
        color = colormap[i+1]
        if i == 3; linestyle = :solid; marker=:circle end
        if i == 2; linestyle = dash; marker=:utriangle end
        if i == 1; linestyle = :dash; marker=:rect end
        if i == 0; linestyle = :dot; marker=:dtriangle end
        lines!(ax, energies, y ./ df.total; label="ν = $(i)", color, linestyle)
        scatter!(ax, energies, y ./ df.total; label="ν = $(i)", color, marker)
    end
end

function plot_fig(model_id)

    fig = Figure(resolution=(JamesPlots.COLUMN_WIDTH, 1.2JamesPlots.RESOLUTION[2]), figure_padding=(1,1,1,1))
    axes = [MyAxis(fig[i,j], limits=(0.10, 1.10, -0.1, 1.1)) for i=1:2, j=1:2]

    hidexdecorations!.(axes[1:end-1,:]; ticks=false, minorticks=false)
    hideydecorations!.(axes[:,2:end]; ticks=false, minorticks=true)
    linkyaxes!(axes...)

    initial_state = 3

    methods = ["MDEF" "IESH"; "Ehrenfest" "BCME"]

    vibrational_number = initial_state
    for i in CartesianIndices(methods)
        method = methods[i]
        Label(fig[i[1], i[2]], method; tellwidth=false, tellheight=false, valign=:top, halign=:left, padding=(10,10,10,15))
        if method == "IESH"
            input = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt-decoherence=EDC.csv")
            plot_entry!(axes[i], input, method)
        else
            input = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt.csv")
            plot_entry!(axes[i], input, method)
        end
    end

    Label(fig[:,0], "Probability", rotation=π/2, padding=(0,3,0,0))
    Label(fig[3,1:2], "Incidence energy /eV")

    labels = ["ν = $i" for i=0:3]
    elements = [
        [LineElement(linestyle=:dot, color=:black), MarkerElement(marker=:dtriangle, color=:white)],
        [LineElement(linestyle=:dash, color=:black), MarkerElement(marker=:rect, color=:white)],
        [LineElement(linestyle=dash, color=:black), MarkerElement(marker=:utriangle, color=:white)],
        [LineElement(linestyle=:solid, color=:black), MarkerElement(marker=:circle, color=:white)],
    ]
    Legend(fig[0,1:2], elements, labels, orientation=:horizontal)

    return fig
end

save_figure(plotsdir("fig6_vibrational_deexcitation_NOAg_nu=3.png"), ()->plot_fig("NOAg"); ncolors=4, fontsize=9, markersize=5)
save_figure(plotsdir("fig7_vibrational_deexcitation_NOAu_nu=3.png"), ()->plot_fig("NOAu"); ncolors=4, fontsize=9, markersize=5)
