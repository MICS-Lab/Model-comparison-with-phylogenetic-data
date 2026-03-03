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

mutable struct point 
    p2::Float64
    p0::Float64
    T_m::Union{Float64,Int}
end

mutable struct perturbation 
    sigma1s::Float64
    sigma1T::Float64
    sigma2s::Float64
    sigma2T::Float64
    sigma3p2::Float64
    sigma3p0::Float64
    sigma3T::Float64
end


    function one_step_threaded(n,file, data, dataCF, nsample, N, k, L, l, LTT_threshold, CF_threshold, current_step, n_models,sigma,name)
        if current_step!=1 
            df = CSV.read(file, DataFrame)
        else
            df=[]
        end
        new_file = "$(name)_step_$(current_step)_$n.csv" 
        println("job n $n started")
        open(new_file, "w") do io
            for i in 1:nsample 
                    while true 
                        if current_step == 1
                            model_choice = mod(i - 1, n_models) + 1
                        else
                            model_choice = rand(1:n_models)
                        end 
                        #println("sampled model: $model_choice")
                        x = sample_parameter(df, model_choice, current_step,L,sigma)
                        # println(io, "$(x.p2),$(x.p0),$(x.T_m),$model_choice,sampled")
                        # flush(io)
                        LTT, CF = run(x, model_choice, current_step, data, dataCF, N, k, L, l)
                        #println("model run")
                        diff = piecewise_difference_g(data[1], data[2], LTT[1], LTT[2])             
                        dist = area(diff[1], diff[2]) / k             
                        CFd = abs(CF - dataCF)   
                        #d_inf=max(CFd ./0.15189758404352943 , dist ./121.0088720943085)   
                        #d_sum=CFd ./0.15189758404352943 + dist ./121.0088720943085 
                        if  dist < LTT_threshold && CFd < CF_threshold #d_sum < LTT_threshold
                            weight = calculate_weight(x, model_choice, current_step, L, df,sigma)
                            p2, p0, T_m = x.p2, x.p0, x.T_m  
                            println(io, "$p2,$p0,$T_m,$model_choice,$dist,$CFd,$weight")
                            flush(io)
                            break 
                        end
                    end
            end
        end
    end


function sample_parameter(df, model_choice, current_step,L,sigma)
    x=point(0.0,0.0,0)
    if current_step ==1 #sample from prior
        if model_choice==1 #VanEgeren
            x.p2=rand()*2
            x.T_m=round(Int,rand(Uniform(0,L-3)))
        elseif model_choice ==2 #Williams
            x.p2=rand()*2
            x.T_m=rand(Uniform(0,L-3))
        elseif model_choice==3
            x.p2=rand()
            x.p0=rand()*min(1-x.p2,x.p2)
            x.T_m=rand(Uniform(0,L-3))
        end 
    else 
        mask = df.model .==model_choice 
        original_indices = findall(mask)
        weighted_sampler=Categorical(Float64.(df.weight[mask]))
        masked_index=rand(weighted_sampler)
        original_index=original_indices[masked_index]
        p2=df.p2[original_index]
        p0=df.p0[original_index]
        T_m=df.T_m[original_index]
        if model_choice ==1 
            sample = -1
            while sample < 0 || sample > L-3
                sample = round(Int,rand(Uniform(T_m - sigma.sigma1T, T_m + sigma.sigma1T)))
            end
            x.T_m=sample
            
            sample = -1
            while sample < 0 || sample > 2 
                sample = rand(Uniform(p2 - sigma.sigma1s, p2 + sigma.sigma1s))
            end 
            x.p2 = sample
        elseif model_choice==2
            sample = -1
            while sample < 0 || sample > L-3
                sample = rand(Normal(T_m, sigma.sigma2T))
            end
            x.T_m = sample
            
            sample = -1
            while sample < 0 || sample > 2 
                sample = rand(Uniform(p2 - sigma.sigma2s, p2 + sigma.sigma2s))
            end 
            x.p2 = sample
        elseif model_choice==3
            sample = -1
            attempts = 0
            while (sample < 0 || sample > L-3) && attempts < max_iterations
                sample = rand(Normal(T_m, sigma.sigma3T))
                attempts += 1
            end
            
            # Fallback: clamp to valid range if sampling fails
            if sample < 0 || sample > L-3
                sample = clamp(T_m + randn() * sigma.sigma3T, 0, L-3)
            end
            x.T_m = sample
            
            # Sample p2 and p0 jointly with bounded attempts
            valid_sample = false
            new_p2, new_p0 = p2, p0
            joint_attempts = 0
            
            while !valid_sample && joint_attempts < max_iterations
                # Sample p2 with bounded attempts
                sample_p2 = -1
                p2_attempts = 0
                p2_lower = max(0.001, p2 - sigma.sigma3p2)
                p2_upper = min(1, p2 + sigma.sigma3p2)
                
                # Check if p2 sampling range is valid
                if p2_lower > p2_upper
                    sample_p2 = clamp(p2, 0.001, 1)
                else
                    while (sample_p2 < 0.001 || sample_p2 > 1) && p2_attempts < max_iterations
                        sample_p2 = rand(Uniform(p2_lower, p2_upper))
                        p2_attempts += 1
                    end
                    # Fallback for p2
                    if sample_p2 < 0.001 || sample_p2 > 1
                        sample_p2 = clamp(p2, 0.001, 1)
                    end
                end
                
                # Sample p0 with bounded attempts
                sample_p0 = -1
                p0_attempts = 0
                p0_lower = max(0.001, p0 - sigma.sigma3p0)
                p0_upper = min(sample_p2, p0 + sigma.sigma3p0)
                
                # Check if p0 sampling range is valid
                if p0_lower > p0_upper
                    sample_p0 = clamp(p0, 0.001, sample_p2)
                else
                    while (sample_p0 < 0.001 || sample_p0 > sample_p2) && p0_attempts < max_iterations
                        sample_p0 = rand(Uniform(p0_lower, p0_upper))
                        p0_attempts += 1
                    end
                    # Fallback for p0
                    if sample_p0 < 0.001 || sample_p0 > sample_p2
                        sample_p0 = clamp(p0, 0.001, sample_p2)
                    end
                end
                
                # Check if joint constraint is satisfied
                if sample_p2 + sample_p0 < 1
                    new_p2, new_p0 = sample_p2, sample_p0
                else
                    # Try to adjust to satisfy constraint
                    total = sample_p2 + sample_p0
                    if total >= 1
                        # Scale down proportionally
                        scale_factor = 0.99 / total
                        sample_p2 *= scale_factor
                        sample_p0 *= scale_factor
                        new_p2, new_p0 = sample_p2, sample_p0
                    end
                end
                
                ε = 1e-6  # your desired threshold — adjust as needed
                if abs(new_p2 - new_p0) <= ε
                    valid_sample = false  # reject this sample
                    continue            # go to next iteration
                else
                    valid_sample = true # accept
                end
                
                joint_attempts += 1
            end
            
            # Final fallback: use clamped original values
            if !valid_sample
                new_p2 = clamp(p2, 0.001, 1)
                new_p0 = clamp(p0, 0.001, min(new_p2, 1 - new_p2))
                @warn "Perturbation kernel failed to find valid sample after $max_iterations attempts, using fallback values"
            end
            
            x.p2 = new_p2
            x.p0 = new_p0
        end
    end
    return x
end

function run(x, model_choice, current_step, data, dataCF, N, k, L, l)
    if model_choice==1
        a = nothing
        while a === nothing  
            a = sim_LTT(round(Int,(L - x.T_m)), N, x.p2, k, L, l)
        end
        return a[1],a[2]
    elseif model_choice ==2
        S_result = nothing 
        while S_result === nothing
            S_result = run_selection_sim(N,x.T_m, L, x.p2)
        end
        tree = construct_mutated_tree(S_result, k, l)
        LTT = LTT_plot(tree, l, L, 0.0)
        CF = S_result.CF
        return LTT,CF 
    elseif model_choice==3
        p = Parameters(1/30, x.p0, x.p2, x.T_m * 365)
        l_result = nothing 
        while l_result === nothing
            l_result = tree_simulation(p, L * 365, 0, 2000, k)
        end
        LTT=[(x.T_m .+ l_result[1][1] ./ 365) .* l, l_result[1][2]]
        return LTT,l_result[2]
    end 
end 

# Add this constant at the top of your file
const max_iterations = 1000

function calculate_weight(x, model_choice, current_step, L, df, sigma) 
    if current_step == 1 
        wt = 1.0
    else 
        if model_choice == 1 
            su = 0.0
            for i in 1:nrow(df)  # Use nrow(df) instead of length(df)
                if df.model[i] == 1
                    try
                        weight_contrib = df.weight[i] * 
                            pdf(Truncated(Uniform(df.p2[i] - sigma.sigma1s, df.p2[i] + sigma.sigma1s), 0, 2), x.p2) *
                            pdf(Truncated(Uniform(df.T_m[i] - sigma.sigma1T, df.T_m[i] + sigma.sigma1T), 0, L), x.T_m)  
                        su += weight_contrib
                    catch e
                        println("Numerical issue in weight calculation: ", e)
                        continue
                    end
                end
            end 
            prior = 1/2 * (L-3) / 2
            wt = prior / su  
            
        elseif model_choice == 2
            su = 0.0
            for i in 1:nrow(df)
                if df.model[i] == 2
                    try
                        weight_contrib = df.weight[i] *
                            pdf(Truncated(Uniform(df.p2[i] - sigma.sigma2s, df.p2[i] + sigma.sigma2s), 0, 2), x.p2) *  # Fixed: sigma1s -> sigma2s
                            pdf(Truncated(Normal(df.T_m[i], sigma.sigma2T), 0, L), x.T_m)
                        su += weight_contrib
                    catch e
                        continue
                    end
                end
            end 
            prior = 1/2 * L / 2
            wt = prior / su 
            
        elseif model_choice == 3
            su = 0.0
            for i in 1:nrow(df)
                if df.model[i] == 3
                    try
                        weight_contrib = df.weight[i] *
                            pdf(Truncated(Uniform(df.p2[i] - sigma.sigma3p2, df.p2[i] + sigma.sigma3p2), 0, 1), x.p2) *
                            pdf(Truncated(Uniform(df.p0[i] - sigma.sigma3p0, df.p0[i] + sigma.sigma3p0), 0, 1), x.p0) *
                            pdf(Truncated(Normal(df.T_m[i], sigma.sigma3T), 0, L), x.T_m)
                        su += weight_contrib
                    catch e
                        println("Numerical issue in weight calculation: ", e)
                        continue
                    end
                end
            end 
            prior = 1/2 * L / 2  
            wt = prior / su  
        end 
    end
    return wt 
end