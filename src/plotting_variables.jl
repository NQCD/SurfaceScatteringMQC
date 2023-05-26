
using JamesPlots

dash = [0.0, 3.2, 4.0]

markers = Dict(
    "Adiabatic"=>:circle,
    "MDEF"=>:utriangle,
    "IESH"=>:rect,
    "Ehrenfest"=>:dtriangle,
    "BCME"=>:diamond,
)

linestyles = Dict(
    "Adiabatic"=>:dash,
    "MDEF"=>dash,
    "IESH"=>:solid,
    "Ehrenfest"=>[0.0, 1.0, 1.8],
    "BCME"=>:dot,
)

method_symbols = Dict(
    "Adiabatic"=>:Classical,
    "MDEF"=>:DiabaticMDEF,
    "IESH"=>:AdiabaticIESH,
    "Ehrenfest"=>:EhrenfestNA,
    "BCME"=>:BCME,
)

color_dict = Dict(
    "MDEF"=>COLORS[1],
    "IESH"=>COLORS[3],
    "Ehrenfest"=>COLORS[4],
    "BCME"=>COLORS[5],
    "Adiabatic"=>COLORS[2],
)