using DrWatson
@quickactivate "ElbowDynamics"
using DataFrames
using CairoMakie
using JamesPlots
using Unitful, UnitfulAtomic
using CSV, DataFrames

include(srcdir("plotting_variables.jl"))

M = 200
bandwidth = 100
dt = 0.25

function plot_entry!(ax, input, method)
    df = CSV.read(input, DataFrame)
    for row in Iterators.reverse(eachrow(df))
        x = Float64[]
        y = Int[]
        minstate = -1
        maxstate = -1
        for i=0:20
            occurrences = convert(Int, row[Symbol(:state, i)])
            if (minstate == -1) && occurrences != 0
                minstate = i
            end
            if occurrences != 0
                maxstate = i
            end
            for j=1:occurrences
                push!(x, row.energies)
                push!(y, i)
            end
        end
        counts = [convert(Int, row[Symbol(:state, i)]) for i=0:20]
        i_end = findfirst(x->x!=0, reverse(counts))
        if i_end === nothing
            i_end = 20
        end
        i_end = 1+length(counts) - i_end + 1
        i_end = min(21, i_end)
        counts /= sum(counts)
        if minstate == 0
            i_start = 1
        else
            i_start = minstate
        end
        offset = row.energies .* 0.75
        band!(ax, collect(0:20)[i_start:i_end], offset, counts[i_start:i_end] .+ offset; color=color_dict[method])
        lines!(ax, collect(0:20)[i_start:i_end], counts[i_start:i_end] .+ offset; color=:black, linewidth=1.0)
    end
end

function plot_fig(model_id, gamma)

    fig = Figure(resolution=(JamesPlots.COLUMN_WIDTH, 1.3JamesPlots.RESOLUTION[2]), figure_padding=(1,1,1,1))
    axes = [MyAxis(fig[i,j], xlabel="Final vibrational state", limits=(0, 19, 0.1, nothing), yticksvisible=false, yminorticksvisible=false) for i=1:2, j=1:2]

    hidexdecorations!.(axes[1:end-1,:]; ticks=false, minorticks=false)
    hideydecorations!.(axes[:,:]; ticks=true, minorticks=true)
    linkyaxes!(axes...)

    initial_state = 16

    methods = ["MDEF" "IESH"; "Ehrenfest" "BCME"]

    vibrational_number = initial_state
    for i in CartesianIndices(methods)
        method = methods[i]
        Label(fig[i[1], i[2]], method; tellwidth=false, tellheight=false, valign=:top, halign=0.6, padding=(5,5,5,5))
        vlines!(axes[i], initial_state; linestyle=dash, color=:black, linewidth=1)
        if method == "IESH"
            input = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt-decoherence=EDC.csv")
            plot_entry!(axes[i], input, method)
        else
            input = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt.csv")
            plot_entry!(axes[i], input, method)
        end
    end

    Label(fig[:,0], "Probability", rotation=π/2, padding=(0,3,0,0))

    if (model_id == "NOAg")
        Label(fig[1,1], L"$E_\text{kin}"; tellwidth=false, tellheight=false, halign=0.60, valign=0.20)
        arrows!(axes[1,1], [9], [0.25*0.75], [1.4], [0.6*0.75], arrowsize=10)
    end


    return fig
end

save_figure(plotsdir("fig8_vibrational_deexcitation_NOAg_1.5.png"), ()->plot_fig("NOAg", 1.5); ncolors=4, fontsize=9)
save_figure(plotsdir("fig9_vibrational_deexcitation_NOAu_1.5.png"), ()->plot_fig("NOAu", 1.5); ncolors=4, fontsize=9)
save_figure(plotsdir("fig11_vibrational_deexcitation_NOAu_15.0.png"), ()->plot_fig("NOAu", 15.0); ncolors=4, fontsize=9)
save_figure(plotsdir("fig12_vibrational_deexcitation_NOAu_0.15.png"), ()->plot_fig("NOAu", 0.15); ncolors=4, fontsize=9)
