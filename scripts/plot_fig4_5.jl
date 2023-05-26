using DrWatson
@quickactivate "SurfaceScatteringMQC"
using DataFrames, CSV
using CairoMakie
using JamesPlots
using Unitful, UnitfulAtomic

include(srcdir("plotting_variables.jl"))

M = 200
gamma = 1.5

function plot_fig(incidence_energy)

    fig = Figure(resolution=(JamesPlots.COLUMN_WIDTH, 1.5JamesPlots.RESOLUTION[2]), figure_padding=(1,1,1,1))
    toprow = [MyAxis(fig[i,j], limits=(-0.5, 5.5, nothing, nothing), xminorticksize=0, xticks=0:5) for i=1:1, j=1:2]
    if incidence_energy == 1.0
        bottomrow = [MyAxis(fig[i,j], limits=(-2, 22, 0, nothing), xminorticksize=0, xticks=0:4:20) for i=2:2, j=1:2]
    else
        bottomrow = [MyAxis(fig[i,j], limits=(-2, 22, 0, nothing), xminorticksize=0, xticks=0:4:20) for i=2:2, j=1:2]
    end
    axes = [toprow; bottomrow]

    linkyaxes!(axes[1,:]...)
    linkyaxes!(axes[2,:]...)
    hideydecorations!.(axes[:,2:end]; ticks=false, minorticks=false)

    for (col, model_id) in enumerate(["NOAu", "NOAg"])
    for (row_i, vibrational_number) in enumerate([3, 16])

        vlines!(axes[row_i,col], vibrational_number; color=:black, linestyle=dash, linewidth=1)
        if (vibrational_number == 3)
            extra_states = 2
        else
            extra_states = 4
        end

        input = projectdir("figure_data",
            "model=$model_id-method=IESH-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=100-dt=0.25-decoherence=none.csv")
        df = CSV.read(input, DataFrame)
        display(df)
        row = subset(df, :energies => a -> a .== incidence_energy)[1,:]
        probabilities = collect(row[2:22])
        errs = collect(row[23:end])
        probabilities ./= sum(probabilities)

        avg = sum((0:20) .* probabilities[1:21])
        errorbars!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], errs[1:vibrational_number+extra_states+1]; whiskerwidth=3, linewidth=1, color=COLORS[6])
        lines!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], color=COLORS[6], label="Standard IESH")
        scatter!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], color=COLORS[6], label="Standard IESH")

        input = projectdir("figure_data",
            "model=$model_id-method=IESH-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=100-dt=0.25-decoherence=EDC.csv")
        df = CSV.read(input, DataFrame)
        row = subset(df, :energies => a -> a .== incidence_energy)[1,:]
        probabilities = collect(row[2:22])
        probabilities ./= sum(probabilities)
        errs = collect(row[23:end])

        avg = sum((0:20) .* probabilities[1:21])
        errorbars!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], errs[1:vibrational_number+extra_states+1]; whiskerwidth=3, linewidth=1, color=COLORS[3])
        lines!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], color=COLORS[3], label="EDC IESH")
        scatter!(axes[row_i,col], 0:vibrational_number+extra_states, probabilities[1:vibrational_number+extra_states+1], color=COLORS[3], label="EDC IESH")
    end
    end
    Label(fig[:,0], "Probability", rotation=π/2, padding=(0, 2, 0, 0))
    Label(fig[end+1,1:2], "Final vibrational state")
    Label(fig[0,1], "NO/Au", tellwidth=false, tellheight=false, halign=:left, padding=(0, 0, -2, 0))
    Label(fig[0,2], "NO/Ag", tellwidth=false, tellheight=false, halign=:right, padding=(0, 0, -2, 0))

    Legend(fig[0,1:2], axes[1], orientation=:horizontal, merge=true, margin=(0, 0, -2, 0))
    rowgap!(fig.layout, 2)

    return fig
end

save_figure(plotsdir("fig4_decoherence_0.2eV.png"), () -> plot_fig(0.2); ncolors=6, fontsize=9, markersize=5)
save_figure(plotsdir("fig5_decoherence_1.0eV.png"), () -> plot_fig(1.0); ncolors=6, fontsize=9, markersize=5)
