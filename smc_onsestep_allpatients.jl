using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  

include("input/prove_dev.jl")
include("input/population_model_modif.jl")
include("input/ABCmethod.jl")
include("input/modelVE.jl")
include("input/cellSimulation.jl")




function one_step_allpatients(n, file, nsample, patient_data, current_step, n_models, npatients, N)
    df = current_step != 1 ? CSV.read(file, DataFrame) : DataFrame()
    new_file = "all_patients_step$(current_step)_$n.csv"

    # Pre-allocate reusable buffers
    patient_params = Vector{Any}(undef, npatients)  # ← Reused across samples
    results = Vector{Tuple{Float64, Float64, String}}(undef, npatients)  # ← Pre-allocated

    # Precompute patient indices or filters if possible
    if current_step !=1 
        patient_filters = [df.name .== p.name for p in patient_data]
    end

    open(new_file, "w") do io
        # println(io, "p2,p0,Tm,model,dist,CFd,weight,patient")
        # println("running")

        for i in 1:nsample
            while true
                model_choice = current_step == 1 ? mod(i - 1, n_models) + 1 : rand(1:n_models)
                #model_choice = current_step == 1 ? mod(i - 1, n_models) + 1 : rand([1, 3])
                all_pass = true

                # Reuse patient_params and results
                for j in 1:npatients
                    p = patient_data[j]
                    df_patient = current_step == 1 ? DataFrame() : df[patient_filters[j], :]

                    # Assume sample_parameter! is in-place or returns lightweight
                    x = sample_parameter(df_patient, model_choice, current_step, p.L, p.sigma)
                    patient_params[j] = x

                    LTT, CF = run(x, model_choice, current_step, p.dataLTT, p.dataCF, N, p.k, p.L, p.l)

                    diff = piecewise_difference_g(p.dataLTT[1], p.dataLTT[2], LTT[1], LTT[2])
                    dist = area(diff[1], diff[2]) / p.k
                    CFd = abs(CF - p.dataCF)

                    results[j] = (dist, CFd, p.name)

                    if dist ≥ p.LTT_threshold || CFd ≥ p.CF_threshold
                        all_pass = false
                        break
                    end
                end

                if all_pass
                    weight = 1.0
                    for j in 1:npatients
                        p = patient_data[j]
                        df_patient = current_step == 1 ? DataFrame() : df[patient_filters[j], :]
                        x = patient_params[j]
                        weight *= calculate_weight(x, model_choice, current_step, p.L, df_patient, p.sigma)
                    end

                    for j in 1:npatients
                        p = patient_data[j]
                        x = patient_params[j]
                        dist, CFd, name = results[j]
                        println(io, "$(x.p2),$(x.p0),$(x.T_m),$model_choice,$dist,$CFd,$weight,$name")
                    end
                    flush(io)
                    break
                end
            end
        end
    end
end
