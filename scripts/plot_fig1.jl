using DrWatson
@quickactivate "SurfaceScatteringMQC"
using NQCModels
using CairoMakie
using JamesPlots
using Unitful, UnitfulAtomic
using SurfaceScatteringMQC
using SurfaceScatteringMQC: FittingPoint

using CSV, DataFrames
using DelimitedFiles

function load_fitting_data(modelid)
    e0_points = FittingPoint[]
    h00_points = FittingPoint[]
    h11_points = FittingPoint[]

    dir = datadir("fitting", modelid, "r=1.6")
    df = CSV.read(joinpath(dir, "e0.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(e0_points, FittingPoint(1.6, row[1], row[2]))
    end
    df = CSV.read(joinpath(dir, "h00.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h00_points, FittingPoint(1.6, row[1], row[2]))
    end
    df = CSV.read(joinpath(dir, "h11.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h11_points, FittingPoint(1.6, row[1], row[2]))
    end

    dir = datadir("fitting", modelid, "r=1.17")
    df = CSV.read(joinpath(dir, "e0.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(e0_points, FittingPoint(1.17, row[1], row[2]))
    end
    df = CSV.read(joinpath(dir, "h00.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h00_points, FittingPoint(1.17, row[1], row[2]))
    end
    df = CSV.read(joinpath(dir, "h11.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h11_points, FittingPoint(1.17, row[1], row[2]))
    end

    dir = datadir("fitting", modelid, "z=1.6")
    df = CSV.read(joinpath(dir, "e0.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(e0_points, FittingPoint(row[1], 1.6, row[2]))
    end
    df = CSV.read(joinpath(dir, "h00.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h00_points, FittingPoint(row[1], 1.6, row[2]))
    end
    df = CSV.read(joinpath(dir, "h11.csv"), DataFrame; header=false)
    for row in eachrow(df)
        push!(h11_points, FittingPoint(row[1], 1.6, row[2]))
    end

    return e0_points, h00_points, h11_points
end

function U0(model, x, z)
    r = austrip.([x z] .* u"Å")
    return ustrip(auconvert(u"eV", NQCModels.potential(model, r)[1,1]))
end

function U1(model, x, z)
    r = austrip.([x z] .* u"Å")
    return ustrip(auconvert(u"eV", NQCModels.potential(model, r)[2,2]))
end

function plot_one_model!(axes, modelid)
    e0_points, h00_points, h11_points = load_fitting_data(modelid)

    bandwidth = austrip(100u"eV")
    gamma = austrip(1.5u"eV")
    metalmodel = SurfaceScatteringMQC.get_model(modelid, 100, gamma, bandwidth)
    model = metalmodel.model

    xgrid = range(0.9, 2.2, length=40)
    ygrid = range(0.5, 6, length=40)

    r0 = 1.17
    points = filter(point->point.x==r0, h00_points)
    scatter!(axes[1], [point.z for point in points], [point.target_energy for point in points])
    points = filter(point->point.x==r0, h11_points)
    scatter!(axes[1], [point.z for point in points], [point.target_energy for point in points])
    points = filter(point->point.x==r0, e0_points)
    scatter!(axes[1], [point.z for point in points], [point.target_energy for point in points])
    lines!(axes[1], ygrid, y->U0(model, r0, y), label="U₀")
    lines!(axes[1], ygrid, y->U1(model, r0, y), label="U₁")
    lines!(axes[1], ygrid, y->SurfaceScatteringMQC.ground_state_potential(metalmodel, r0, y, (;x=1.15, z=6.0)), label="E₀")

    r0 = 1.6
    points = filter(point->point.x==r0, h00_points)
    scatter!(axes[2], [point.z for point in points], [point.target_energy for point in points])
    points = filter(point->point.x==r0, h11_points)
    scatter!(axes[2], [point.z for point in points], [point.target_energy for point in points])
    points = filter(point->point.x==r0, e0_points)
    scatter!(axes[2], [point.z for point in points], [point.target_energy for point in points])
    lines!(axes[2], ygrid, y->U0(model, r0, y))
    lines!(axes[2], ygrid, y->U1(model, r0, y))
    lines!(axes[2], ygrid, y->SurfaceScatteringMQC.ground_state_potential(metalmodel, r0, y, (;x=1.15, z=6.0)))

    points = filter(point->point.z==1.6, h00_points)
    scatter!(axes[3], [point.x for point in points], [point.target_energy for point in points])
    points = filter(point->point.z==1.6, h11_points)
    scatter!(axes[3], [point.x for point in points], [point.target_energy for point in points])
    points = filter(point->point.z==1.6, e0_points)
    scatter!(axes[3], [point.x for point in points], [point.target_energy for point in points])
    lines!(axes[3], xgrid, x->U0(model, x, 1.6))
    lines!(axes[3], xgrid, x->U1(model, x, 1.6))
    lines!(axes[3], xgrid, x->SurfaceScatteringMQC.ground_state_potential(metalmodel, x, 1.6, (;x=1.15, z=6.0)))
end

function plot_both_models()
    fig = Figure(resolution=(JamesPlots.COLUMN_WIDTH, 1.7JamesPlots.RESOLUTION[2]), figure_padding=(1,1,1,1))
    leftcol = [MyAxis(fig[i,1], limits=(nothing, nothing, nothing, 6.5)) for i in 1:3]
    rightcol = [MyAxis(fig[i,2], limits=(nothing, nothing, nothing, 6.5)) for i in 1:3]
    axes = hcat(leftcol, rightcol)

    linkaxes!(axes[1,:]...)
    linkaxes!(axes[2,:]...)
    linkaxes!(axes[3,:]...)
    hidexdecorations!.(axes[1,:]; minorticks=false, ticks=false)
    hideydecorations!.(axes[:,2:end]; minorticks=false, ticks=false)

    Label(fig[2,:,Bottom()], "Height above the surface /Å", padding=(5,5,4,11))
    Label(fig[3,:,Bottom()], "Bond length /Å", padding=(5,5,0,11))
    Label(fig[1:2,0], "Energy /eV"; rotation=π/2, padding=(0, 3, 0, 0), tellheight=false)
    Label(fig[3,0], "Energy /eV"; rotation=π/2, padding=(0, 3, 0, 0), tellheight=false)

    plot_one_model!(axes[:,1], "NOAu")
    plot_one_model!(axes[:,2], "NOAg")

    Legend(fig[0,1:2], axes[1], orientation=:horizontal, margin=(0, 0, 0, 0), merge=true)
    Label(fig[0,1], "NO/Au", halign=0.2, tellwidth=false)
    Label(fig[0,2], "NO/Ag", halign=0.8, tellwidth=false)

    Label(fig[1,3], "r = 1.17 Å", rotation=-π/2, tellheight=false, padding=(0, 2, 0, 0))
    Label(fig[2,3], "r = 1.6 Å", rotation=-π/2, tellheight=false, padding=(0, 2, 0, 0))
    Label(fig[3,3], "z = 1.6 Å", rotation=-π/2, tellheight=false, padding=(0, 2, 0, 0))

    height = 12
    width = 12
    boxhalign = 0.92
    boxvalign = 0.1
    texthalign = 0.92
    textvalign = 0.1
    padding = (0, 3, 0, 3)
    Box(fig[1,1]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[1,1], "A"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)
    Box(fig[2,1]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[2,1], "B"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)
    Box(fig[3,1]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[3,1], "C"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)
    Box(fig[1,2]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[1,2], "D"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)
    Box(fig[2,2]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[2,2], "E"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)
    Box(fig[3,2]; tellwidth=false, tellheight=false, valign=boxvalign, halign=boxhalign, height, width, color=:white)
    Label(fig[3,2], "F"; tellwidth=false, tellheight=false, valign=textvalign, halign=texthalign, justification=:center, padding)

    return fig
end

save_figure(plotsdir("fig1_model_fitting.png"), plot_both_models; fontsize=9)
