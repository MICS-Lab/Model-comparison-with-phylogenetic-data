
@everywhere begin
    using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
end

@everywhere include("prove_dev.jl")
@everywhere include("population_model_modif.jl")
@everywhere include("ABCmethod.jl")
@everywhere include("modelVE.jl")
@everywhere include("cellSimulation.jl")

@everywhere abstract type AbstractModel end

@everywhere mutable struct ModelState
    name::String
    results::Union{Dict{String, Any}, Nothing}
    sigmas::Float64
    sigmap0::Float64
    sigmaT::Float64
    id::Int
end

@everywhere mutable struct VEmodel <: AbstractModel
    state::ModelState 
    s::Float64
    g::Int64
end

@everywhere mutable struct Wmodel <: AbstractModel
    state::ModelState
    s::Float64
    T_m::Float64 
end 

@everywhere mutable struct GHmodel <: AbstractModel
    state::ModelState 
    p2::Float64
    p0::Float64
    T_m::Float64 
    alpha::Float64
    Nmax::Int64
end

@everywhere mutable struct prob_distribution
    model::Vector{Int}
    T::Vector{Union{Int,Float64}}
    p2::Vector{Float64}
    p0::Vector{Float64}
    points::Int
    dist::Vector{Float64}
    weight::Vector{Float64}
    acc_rate::Float64
    CFdist::Vector{Float64}
    B::Int
    t::Vector{Float64}
    mem::Vector{Float64}
end

@everywhere function run(model::VEmodel, N, k, L, l)
    a = nothing
    while a === nothing  # Fixed: use === instead of ==
        a = sim_LTT(model.g, N, model.s, k, L, l)
    end
    model.state.results = Dict("LTT" => a[1], "CF" => a[2])
    return nothing
end

@everywhere function run(model::Wmodel, N, k, L, l)
    S_result = nothing 
    while S_result === nothing
        S_result = run_selection_sim(N, model.T_m, L, model.s)
    end
    tree = construct_mutated_tree(S_result, k, l)
    LTT = LTT_plot(tree, l, L, 0.0)
    CF = S_result.CF
    model.state.results = Dict("LTT" => LTT, "CF" => CF)
    return nothing
end

@everywhere function run(model::GHmodel, N, k, L, l)
    p = Parameters(model.alpha, model.p0, model.p2, model.T_m * 365)
    l_result = nothing 
    while l_result === nothing
        l_result = tree_simulation(p, L * 365, 0, model.Nmax, k)
    end
    model.state.results = Dict("LTT" => [(model.T_m .+ l_result[1][1] ./ 365) .* l, l_result[1][2]], "CF" => l_result[2])
    return nothing
end

@everywhere function perturbation_kernel(model::VEmodel, p2, p0, T_m, L)
    sample = 0
    while sample < 3 || sample > L
        sample = round(Int, L - rand(Uniform(T_m - model.state.sigmaT, T_m + model.state.sigmaT)))
    end
    model.g = sample
    
    sample = -1
    while sample < 0 || sample > 2 
        sample = rand(Uniform(p2 - model.state.sigmas, p2 + model.state.sigmas))
    end 
    model.s = sample
    
    return model.s, 0, round(Int, L - model.g)
end

@everywhere function perturbation_kernel(model::Wmodel, p2, p0, T_m, L)
    sample = -1
    while sample < 0 || sample > L
        sample = rand(Normal(T_m, model.state.sigmaT))
    end
    model.T_m = sample
    
    sample = -1
    while sample < 0 || sample > 2 
        sample = rand(Uniform(p2 - model.state.sigmas, p2 + model.state.sigmas))
    end 
    model.s = sample
    
    return model.s, 0, model.T_m
end

@everywhere function perturbation_kernel(model::GHmodel, p2, p0, T_m, L; max_iterations=100)
    # Sample T_m with bounded attempts
    sample = -1
    attempts = 0
    while (sample < 0 || sample > L) && attempts < max_iterations
        sample = rand(Normal(T_m, model.state.sigmaT))
        attempts += 1
    end
    
    # Fallback: clamp to valid range if sampling fails
    if sample < 0 || sample > L
        sample = clamp(T_m + randn() * model.state.sigmaT, 0, L)
    end
    model.T_m = sample
    
    # Sample p2 and p0 jointly with bounded attempts
    valid_sample = false
    new_p2, new_p0 = p2, p0
    joint_attempts = 0
    
    while !valid_sample && joint_attempts < max_iterations
        # Sample p2 with bounded attempts
        sample_p2 = -1
        p2_attempts = 0
        p2_lower = max(0.001, p2 - model.state.sigmas)
        p2_upper = min(1, p2 + model.state.sigmas)
        
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
        p0_lower = max(0.001, p0 - model.state.sigmap0)
        p0_upper = min(sample_p2, p0 + model.state.sigmap0)
        
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
            valid_sample = true
        else
            # Try to adjust to satisfy constraint
            total = sample_p2 + sample_p0
            if total >= 1
                # Scale down proportionally
                scale_factor = 0.99 / total
                sample_p2 *= scale_factor
                sample_p0 *= scale_factor
                new_p2, new_p0 = sample_p2, sample_p0
                valid_sample = true
            end
        end
        
        joint_attempts += 1
    end
    
    # Final fallback: use clamped original values
    if !valid_sample
        new_p2 = clamp(p2, 0.001, 1)
        new_p0 = clamp(p0, 0.001, min(new_p2, 1 - new_p2))
        @warn "Perturbation kernel failed to find valid sample after $max_iterations attempts, using fallback values"
    end
    
    model.p2 = new_p2
    model.p0 = new_p0
    
    return model.p2, model.p0, model.T_m
end

@everywhere function sample_prior(model::VEmodel, L)  # Added L parameter
    model.s = rand() * 2
    model.g = round(Int, rand(Uniform(3, L)))
    return model.s, 0, round(Int, L - model.g)
end

@everywhere function sample_prior(model::Wmodel, L)  # Added L parameter
    model.s = rand() * 2
    model.T_m = rand(Uniform(0, L))
    return model.s, 0, model.T_m
end

@everywhere function sample_prior(model::GHmodel, L)  # Added L parameter
    # d=clamp(rand(Normal(0.017,0.0039)),0,1)
    # q= clamp(1 .-rand(Exponential(0.092)),0,1)
    # model.p2 = clamp(d/(1-q),0,1)
    # model.p0 = clamp(q*d/(1-q),0,min(model.p2, 1 - model.p2)) 
    model.p2 = rand()
    model.p0=rand()* min(model.p2, 1 - model.p2)  
    model.T_m = rand(Uniform(0, L))
    return model.p2, model.p0, model.T_m
end


@everywhere function calculate_weight(model::VEmodel, pd::prob_distribution, L, t, b,new_p2,new_p0,new_T)
    if t == 1
        wt = b
    else 
        su = 0.0
        for i in 1:pd.points
            if pd.model[i] == model.state.id
                try
                    weight_contrib = pd.weight[i] *pdf(Truncated(Uniform(pd.p2[i] - model.state.sigmas, pd.p2[i] + model.state.sigmas), 0, 2),new_p2) *
                        pdf(Truncated(Uniform(pd.T[i] - model.state.sigmaT,pd.T[i] + model.state.sigmaT), 0, L), new_T)
                    su += weight_contrib
                catch e
                    println("Numerical issue in weight calculation: ", e)
                    continue
                end
            end
        end 
        prior = 1/2 * (L-3) / 2
        #wt = su > 0 ? b * prior / su : 0.0
        wt = b * prior / su
    end
    return wt
end

@everywhere function calculate_weight(model::Wmodel, pd::prob_distribution, L, t, b,new_p2,new_p0,new_T)
    if t == 1
        wt = b
    else 
        su = 0.0
        for i in 1:pd.points
            if pd.model[i] == model.state.id
                try
                    weight_contrib = pd.weight[i] *
                    pdf(Truncated(Uniform(pd.p2[i] - model.state.sigmas, pd.p2[i] + model.state.sigmas), 0,2),new_p2)*
                    pdf(Truncated(Normal(pd.T[i], model.state.sigmaT),0,L),new_T)
                    su += weight_contrib
                catch e
                    continue
                end
            end
        end 
        prior = 1/2 * L / 2
        #wt = su > 0 ? b * prior / su : 0.0
        wt = b * prior / su
    end
    return wt
end

@everywhere function calculate_weight(model::GHmodel, pd::prob_distribution, L, t, b,new_p2,new_p0,new_T)
    if t == 1
        wt = b
    else 
        su = 0.0
        for i in 1:pd.points
            if pd.model[i] == model.state.id
                try
                    # weight_contrib = pd.weight[i] *
                    # pdf(Truncated(Uniform(pd.p2[i] - model.state.sigmas, pd.p2[i] + model.state.sigmas), 0, 1),model.p2)*
                    # pdf(Truncated(Uniform(pd.p0[i] - model.state.sigmas, pd.p0[i] + model.state.sigmas), 0,1),model.p0)*
                    #     pdf(Truncated(Normal(pd.T[i], model.state.sigmaT),0,L),model.T_m)
                    weight_contrib = pd.weight[i] *
                    pdf(Truncated(Uniform(pd.p2[i] - model.state.sigmas, pd.p2[i] + model.state.sigmas), 0, 1),new_p2)*
                    pdf(Truncated(Uniform(pd.p0[i] - model.state.sigmap0, pd.p0[i] + model.state.sigmap0), 0,1),new_p0)*
                        pdf(Truncated(Normal(pd.T[i], model.state.sigmaT),0,L),new_T)
                    su += weight_contrib
                catch e
                    println("Numerical issue in weight calculation: ", e)
                    continue
                end
            end
        end 
        
        prior = L / 2
        #wt = su > 0 ? b * prior / su : 0.0
        wt = b * prior / su
    end
    return wt
end


 @everywhere function onesample(i, models::Vector{AbstractModel}, pd::prob_distribution, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    distances = Float64[]     
    CFdistances = Float64[]     
    acc_rate = 0     
    b = 0          
    while true

        n_models = length(models)         
        if reach_convergence == false             
            model_choice = rand(1:n_models)         
        else             
            model_choice = mod(i - 1, n_models) + 1         
        end      
        #n_params = [2, 2, 3]
        # n_params = [1]
        # n_models = length(n_params)

        # # expand into weighted sequence for deterministic cycling
        # expanded_models = reduce(vcat, [fill(m, n_params[m]) for m in 1:n_models])
        # n_expanded = length(expanded_models)

        # if reach_convergence == false
        #     # Random choice, probability proportional to n_params
        #     model_choice = sample(1:n_models, Weights(n_params))
        # else
        #     # Deterministic cycling, but proportional to n_params
        #     model_choice = expanded_models[mod(i - 1, n_expanded) + 1]
        # end           

        # Create a local copy of the model to avoid modifying shared state         
        local_model = deepcopy(models[model_choice])                  
        if t == 1             
            p2, p0, T_m = sample_prior(local_model, L)  # Pass L parameter         
        else              
            mask = pd.model .== model_choice             
            original_indices = findall(mask)             
            if isempty(original_indices)                 
                # If no particles for this model, sample from prior                 
                p2, p0, T_m = sample_prior(local_model, L)             
            else                 
                weighted_sampler = Categorical(Float64.(pd.weight[mask]))                 
                masked_index = rand(weighted_sampler)                 
                original_index = original_indices[masked_index]                 
                p2, p0, T_m = perturbation_kernel(local_model, pd.p2[original_index], pd.p0[original_index], pd.T[original_index], L)             
            end         
        end         
        #println("m: $model_choice, s: $p2, $p0, T: $T_m ")         
        for j in 1:pd.B             
            run(local_model, N, k, L, l)             
            diff = piecewise_difference_g(data[1], data[2], local_model.state.results["LTT"][1], local_model.state.results["LTT"][2])             
            dist = area(diff[1], diff[2]) / k             
            CFd = abs(local_model.state.results["CF"] - dataCF)             
            #println("dist: $dist , CF: $(local_model.state.results["CF"])")             
            if dist < dist_thres && CFd < CFdist                 
                b += 1                 
                push!(distances, dist)                 
                push!(CFdistances, CFd)             
            end             
            acc_rate += 1         
        end                  
        if b > 0              
            models[model_choice]=local_model             
            return [model_choice, p2, p0, T_m, minimum(distances), minimum(CFdistances), acc_rate, b]         
        end     
    end 
end 
@everywhere function update_particle(i::Int, models::Vector{AbstractModel}, pd::prob_distribution, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    time_el = @elapsed f = onesample(i, models, pd, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    new_model = f[1]       
    new_T = f[4]        
    new_p2 = f[2]     
    new_p0 = f[3]       
    new_dist = f[5]     
    new_weight = calculate_weight(models[round(Int64,new_model)], pd, L, t, f[8],new_p2,new_p0,new_T) 
    new_alpha = f[7]     
    new_CF = f[6]          
    return (new_model, new_T, new_p2, new_p0, new_dist, new_CF, new_alpha, new_weight,time_el) 
end 
    
@everywhere function safe_update_particle(i::Int, pd::prob_distribution,
    data, dataCF, N, k, L, l,
    reach_convergence, t,
    dist_thres, CFdist)

    max_tries = 1000
    tries = 0
    result = nothing

    while result === nothing && tries < max_tries
    try
    result = update_particle(i, MODELS, pd, data, dataCF,
    N, k, L, l,
    reach_convergence, t,
    dist_thres, CFdist)
    catch e
    # if something went wrong, retry
    result = nothing
    end
    tries += 1
    end

    if result === nothing
    error("Worker $(myid()): update_particle failed after $max_tries tries")
    end

    return result
end

using Distributed

function onestep(models::Vector{AbstractModel}, pd::prob_distribution,
                 data, dataCF, N, k, L, l,
                 reach_convergence, t, N_points,
                 dist_thres, CFdist)

    n_workers = nprocs() - 1       # exclude master process
    chunk_size = ceil(Int, N_points / n_workers)

    # Split particle indices into chunks
    tasks = [i:min(i+chunk_size-1, N_points) for i in 1:chunk_size:N_points]

    results = pmap(task -> begin
        # Compute particles for this chunk
        map(i -> safe_update_particle(i, pd, data, dataCF,
                                      N, k, L, l,
                                      reach_convergence, t,
                                      dist_thres, CFdist), task)
    end, tasks)

    # Flatten results (list of lists → single list)
    updated = vcat(results...)

    pd.model    = [u[1] for u in updated]
    pd.T        = [u[2] for u in updated]
    pd.p2       = [u[3] for u in updated]
    pd.p0       = [u[4] for u in updated]
    pd.dist     = [u[5] for u in updated]
    pd.CFdist   = [u[6] for u in updated]
    pd.acc_rate = N_points / sum([u[7] for u in updated])
    pd.weight   = [u[8] for u in updated]
    pd.points   = N_points
    pd.t = [u[9] for u in updated]

    normalize_weight(pd)

    return pd
end

using Base.Threads
using StatsBase  # for sample() and Weights()

# Multithreaded version of your onesample function
function onesample_threaded(i, models::Vector{AbstractModel}, pd::prob_distribution, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    distances = Float64[]     
    CFdistances = Float64[]     
    acc_rate = 0     
    b = 0          
    
    while true
        n_models = length(models)         
        if reach_convergence == false             
            model_choice = rand(1:n_models)         
        else             
            model_choice = mod(i - 1, n_models) + 1         
        end   
        # Your weighted model selection logic
        # n_params = [1]
        # n_models = length(n_params)

        # # expand into weighted sequence for deterministic cycling
        # expanded_models = reduce(vcat, [fill(m, n_params[m]) for m in 1:n_models])
        # n_expanded = length(expanded_models)

        # if reach_convergence == false
        #     # Random choice, probability proportional to n_params
        #     model_choice = sample(1:n_models, Weights(n_params))
        # else
        #     # Deterministic cycling, but proportional to n_params
        #     model_choice = expanded_models[mod(i - 1, n_expanded) + 1]
        # end           

        # Create a local copy of the model to avoid modifying shared state         
        local_model = deepcopy(models[model_choice])                  
        
        if t == 1             
            p2, p0, T_m = sample_prior(local_model, L)         
        else              
            mask = pd.model .== model_choice             
            original_indices = findall(mask)             
            if isempty(original_indices)                 
                p2, p0, T_m = sample_prior(local_model, L)             
            else                 
                weighted_sampler = Categorical(Float64.(pd.weight[mask]))                 
                masked_index = rand(weighted_sampler)                 
                original_index = original_indices[masked_index]                 
                p2, p0, T_m = perturbation_kernel(local_model, pd.p2[original_index], pd.p0[original_index], pd.T[original_index], L)             
            end         
        end         
        
        for j in 1:pd.B             
            run(local_model, N, k, L, l)             
            diff = piecewise_difference_g(data[1], data[2], local_model.state.results["LTT"][1], local_model.state.results["LTT"][2])             
            dist = area(diff[1], diff[2]) / k             
            CFd = abs(local_model.state.results["CF"] - dataCF)             
            
            if dist < dist_thres && CFd < CFdist                 
                b += 1                 
                push!(distances, dist)                 
                push!(CFdistances, CFd)             
            end             
            acc_rate += 1         
        end                  
        
        if b > 0              
            # Note: Don't modify shared models array in multithreading
            # Return the local_model information instead
            return [model_choice, p2, p0, T_m, minimum(distances), minimum(CFdistances), acc_rate, b]         
        end     
    end 
end 

# Multithreaded version of update_particle
function update_particle_threaded(i::Int, models::Vector{AbstractModel}, pd::prob_distribution, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    time_el = @elapsed f = onesample_threaded(i, models, pd, data, dataCF, N, k, L, l, reach_convergence, t, dist_thres, CFdist)     
    new_model = f[1]       
    new_T = f[4]        
    new_p2 = f[2]     
    new_p0 = f[3]       
    new_dist = f[5]     
    new_weight = calculate_weight(models[round(Int64,new_model)], pd, L, t, f[8],new_p2,new_p0,new_T)     
    new_alpha = f[7]     
    new_CF = f[6]          
    return (new_model, new_T, new_p2, new_p0, new_dist, new_CF, new_alpha, new_weight,time_el) 
end 

# Safe multithreaded version
function safe_update_particle_threaded(i::Int, models::Vector{AbstractModel}, pd::prob_distribution,
                                      data, dataCF, N, k, L, l,
                                      reach_convergence, t,
                                      dist_thres, CFdist)

    max_tries = 1000
    tries = 0
    result = nothing

    while result === nothing && tries < max_tries
        try
            result = update_particle_threaded(i, models, pd, data, dataCF,
                                            N, k, L, l,
                                            reach_convergence, t,
                                            dist_thres, CFdist)
        catch e
            # if something went wrong, retry
            result = nothing
        end
        tries += 1
    end

    if result === nothing
        error("Thread $(threadid()): update_particle failed after $max_tries tries")
    end

    return result
end

# Main multithreaded onestep function
function onestep_multithreaded(models::Vector{AbstractModel}, pd::prob_distribution,
                              data, dataCF, N, k, L, l,
                              reach_convergence, t, N_points,
                              dist_thres, CFdist)
    
    # Pre-allocate results array
    updated = Vector{Any}(undef, N_points)
    
    # Process all particles in parallel using threads
    @threads for i in 1:N_points
        updated[i] = safe_update_particle_threaded(i, models, pd, data, dataCF,
                                                  N, k, L, l,
                                                  reach_convergence, t,
                                                  dist_thres, CFdist)
    end
    
    # Update particle distribution (same as your original)
    pd.model    = [u[1] for u in updated]
    pd.T        = [u[2] for u in updated]
    pd.p2       = [u[3] for u in updated]
    pd.p0       = [u[4] for u in updated]
    pd.dist     = [u[5] for u in updated]
    pd.CFdist   = [u[6] for u in updated]
    pd.acc_rate = N_points / sum([u[7] for u in updated])
    pd.weight   = [u[8] for u in updated]
    pd.points   = N_points
    pd.t= [u[9] for u in updated]

    normalize_weight(pd)

    return pd
end
@everywhere function normalize_weight(distribution::prob_distribution)
    for m in [1,2,3]
        mask = distribution.model .== m
        summ = sum(distribution.weight[mask])
        if summ > 0
            distribution.weight[mask] ./= summ  # More efficient update
        end
    end
    return distribution
end
function onestep_multithreaded_chunked(models::Vector{AbstractModel}, pd::prob_distribution,
    data, dataCF, N, k, L, l,
    reach_convergence, t, N_points,
    dist_thres, CFdist; chunk_size::Int = 10)

    all_results = Vector{Any}()
    sizehint!(all_results, N_points)

    for chunk_start in 1:chunk_size:N_points
    chunk_end = min(chunk_start + chunk_size - 1, N_points)
    chunk_range = chunk_start:chunk_end
    chunk_results = Vector{Any}(undef, length(chunk_range))

    @threads for (idx, i) in collect(enumerate(chunk_range))
    chunk_results[idx] = safe_update_particle_threaded(i, models, pd, data, dataCF,
                                N, k, L, l,
                                reach_convergence, t,
                                dist_thres, CFdist)
    end

    append!(all_results, chunk_results)

    # Optional: Garbage collection between chunks
    if chunk_end < N_points
    GC.gc(false)
    end
    end

    # Update particle distribution
    pd.model    = [u[1] for u in all_results]
    pd.T        = [u[2] for u in all_results]
    pd.p2       = [u[3] for u in all_results]
    pd.p0       = [u[4] for u in all_results]
    pd.dist     = [u[5] for u in all_results]
    pd.CFdist   = [u[6] for u in all_results]
    pd.acc_rate = N_points / sum([u[7] for u in all_results])
    pd.weight   = [u[8] for u in all_results]
    pd.points   = N_points

    normalize_weight(pd)
    return pd
end

function smc(file_path, models::Vector{AbstractModel}, data, dataCF, N, k, L, l, reach_convergence, N_points, initial_dist, initialCF, maxsteps, tol)
    init_time = time_ns()
    csv_file = file_path
    
    open(csv_file, "w") do file
        println(file, "p2,p0,T_m,weight,dist,CF,model,time(s)")
        flush(file)
        
        pd = prob_distribution(
            Int[],
            Float64[],
            Float64[],
            Float64[],
            N_points,
            Float64[],
            Float64[],
            0.0,
            Float64[],
            1,
            Float64[],
            Float64[]
        )
        
        dist = initial_dist
        CFdist = initialCF
        dist_prec = dist
        CFprec = CFdist
        t = 1
        
        while t <= maxsteps
            #benchmark_approaches(models, pd,
                            # data, dataCF, N, k, L, l, 
                            # reach_convergence, t, N_points, dist, CFdist)
            pd = onestep(models, pd, data, dataCF, N, k, L, l, reach_convergence, t, N_points, dist, CFdist)
            #pd=onestep_multithreaded(models, pd, data, dataCF, N, k, L, l, reach_convergence, t, N_points, dist, CFdist)
            time_elapsed = (time_ns() - init_time) / (1e9 * 60 * 60)
            
            println("Step $t completed in time $time_elapsed hours")
            println(file, "step $t distance $dist CFdistance $CFdist time $time_elapsed alpha $(pd.acc_rate)")
            # if t==maxsteps
            #     println(file, "step $t distance $dist CFdistance $CFdist time $time_elapsed alpha $(pd.acc_rate)")
                for i in 1:N_points
                    println(file, "$(pd.p2[i]),$(pd.p0[i]),$(pd.T[i]),$(pd.weight[i]),$(pd.dist[i]),$(pd.CFdist[i]),$(pd.model[i]), $(pd.t[i])")
                end
                flush(file)
            #end
            
            # Check if any model has too few particles
            for i in 1:length(models)
                if length(pd.model[pd.model .== i]) <= 1
                    println("Model $(models[i].state.name) has no particles, end smc")
                    return pd
                end
            end
            
            # Update model parameters for next iteration
            for i in 1:length(models)
                models[i].state.results = nothing
                model_indices = pd.model .== i
                if sum(model_indices) > 1  # Need at least 2 points for variance
                    models[i].state.sigmas = max(sqrt(var(pd.p2[model_indices])), 0.002)  # Prevent zero variance
                    models[i].state.sigmap0 = max(sqrt(var(pd.p0[model_indices])), 0.002)  
                    models[i].state.sigmaT = max(sqrt(var(pd.T[model_indices])), 1.0)
                else
                    # Use default values if too few particles
                    models[i].state.sigmas = 0.2
                    models[i].state.sigmap0 = 0.2
                    models[i].state.sigmaT = 2.0
                end
                println("Model $(models[i].state.name) sigmas: $(models[i].state.sigmas), sigmap0: $(models[i].state.sigmap0), sigmaT: $(models[i].state.sigmaT)")
            end
            
            # if t >= 5 && (abs(dist - dist_prec) < tol * dist || abs(CFdist - CFprec) < tol * CFdist)
            #     println("Convergence reached at step $t")
            #     break
            # end
            dist = quantile(pd.dist, 0.7)
            CFdist = quantile(pd.CFdist, 0.7)
            t += 1
            if t > 1
                reach_convergence = false
            end
            
            dist_prec = dist
            CFprec = CFdist
        end
        
        return pd
    end
end