
using DrWatson
@quickactivate "ElbowDynamics"
using DataFrames, CSV
using Statistics: mean, std
using StatsBase

include(srcdir("data.jl"))

function process_vibrational_transitions()
    energies = 0.2:0.1:1.0
    vibrational_states = [3, 16]
    gammas = [1.5]
    Ms = [200]
    dts = [0.25]

    models = ["NOAg", "NOAu"]
    methods = ["MDEF", "IESH", "Ehrenfest", "BCME"]
    method_symbols = Dict(
        "Adiabatic"=>:Classical,
        "MDEF"=>:DiabaticMDEF,
        "IESH"=>:AdiabaticIESH,
        "Ehrenfest"=>:EhrenfestNA,
        "BCME"=>:BCME,
    )

    for decoherence in [:EDC]
    for bandwidth in [100]
    for M in Ms
    for dt in dts
    for model_id in models
        @info "Model: $model_id"
        for method in methods
            @info "Method: $method"
            for vibrational_number in vibrational_states
                @info "ν = $vibrational_number"
                for gamma in gammas

                    final_state_probability = zeros(length(energies), 21)
                    final_state_count = zeros(Int, length(energies), 21)
                    err = zero(final_state_probability)
                    average = zeros(length(energies))
                    standarddeviation = zeros(length(energies))
                    q0 = zeros(length(energies))
                    q1 = zeros(length(energies))
                    q2 = zeros(length(energies))
                    q3 = zeros(length(energies))
                    q4 = zeros(length(energies))
                    total = zeros(length(energies))
                    for (i, incident_energy) in enumerate(energies)

                        if method == "IESH"
                            params = Dict{String,Any}(
                                "trajectories"=>2000,
                                "dt"=>dt,
                                "M"=>M,
                                "vibrational_number"=>vibrational_number,
                                "incident_energy"=>incident_energy,
                                "temperature"=>300,
                                "method"=>method_symbols[method],
                                "model_id"=>Symbol(model_id),
                                "bandwidth"=>bandwidth,
                                "Γ"=>gamma,
                                "decoherence"=>decoherence,
                            )
                        else
                            params = Dict{String,Any}(
                                "trajectories"=>2000,
                                "dt"=>dt,
                                "M"=>M,
                                "vibrational_number"=>vibrational_number,
                                "incident_energy"=>incident_energy,
                                "temperature"=>300,
                                "method"=>method_symbols[method],
                                "model_id"=>Symbol(model_id),
                                "bandwidth"=>bandwidth,
                                "Γ"=>gamma,
                            )
                        end

                        try
                        vibrational_state = read_data_file(datadir("elbowdynamics", model_id, lowercase(method)), params, "OutputVibrationalState", Int)
                        filter!(x->x!=-10, vibrational_state)
                        filter!(x->x!=-100, vibrational_state)
                        average[i] = mean(vibrational_state)
                        standarddeviation[i] = std(vibrational_state)
                        quantiles = nquantile(vibrational_state, 4)
                        q0[i], q1[i], q2[i], q3[i], q4[i] = quantiles

                        state_count = [count(x->x == state, vibrational_state) for state in 0:20]
                        n = sum(state_count)
                        total[i] = n

                        final_state_count[i,:] .= state_count
                        final_state_probability[i,:] .=  state_count ./ n

                        for s = 0:20
                            p = final_state_probability[i,s+1]
                            err[i,s+1] = sqrt(p * (1-p) / n) # https://stats.stackexchange.com/questions/396952/standard-error-of-individual-groups-in-a-multinomial-distribution
                        end
                        catch
                            final_state_probability[i,:] .= NaN
                            err[i,:] .= NaN
                        end
                    end
                    headers = [:energies]
                    for i=0:20
                        push!(headers, Symbol(:state, i))
                    end
                    for i=0:20
                        push!(headers, Symbol(:err, i))
                    end
                    push!(headers, :mean)
                    push!(headers, :std)
                    push!(headers, :q0)
                    push!(headers, :q1)
                    push!(headers, :q2)
                    push!(headers, :q3)
                    push!(headers, :q4)
                    push!(headers, :total)
                    df = DataFrame([energies final_state_count err average standarddeviation q0 q1 q2 q3 q4 total], headers)

                    if method == "IESH"
                        output = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt-decoherence=$decoherence.csv")
                    else
                        output = projectdir("figure_data", "model=$model_id-method=$method-v=$vibrational_number-gamma=$gamma-M=$M-bandwidth=$bandwidth-dt=$dt.csv")
                    end
                    CSV.write(output, df)
                end
            end
        end
    end
    end
    end
    end
    end
end

process_vibrational_transitions()
