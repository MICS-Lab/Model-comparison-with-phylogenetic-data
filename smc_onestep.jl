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
include("cellSimulation.jl")

mutable struct point
    p2::Float64
    p0::Float64
    T_m::Union{Float64,Int}
end

# sigma3p0 is replaced by sigma3u. For current_step > 1 it is recomputed
# automatically from the previous population (see compute_sigma_u), so the value
# passed in by the driver in this slot is ignored.
mutable struct perturbation
    sigma1s::Float64
    sigma1T::Float64
    sigma2s::Float64
    sigma2T::Float64
    sigma3p2::Float64
    sigma3u::Float64
    sigma3T::Float64
end

# ---------------------------------------------------------------
# Model 3 helpers: m(x) = min(x, 1-x),  u = p0 / m(p2)
# ---------------------------------------------------------------
m_half(x) = min(x, 1 - x)

function u_from(p2, p0)
    mp = m_half(p2)
    return mp > 0 ? clamp(p0 / mp, 0.0, 1.0) : 0.0
end

const SIGMA_U_MIN = 0.002

# σ_{u,t} = max{0.002, sd of u over the M3 particles of population t-1}
function compute_sigma_u(df)
    mask = df.model .== 3
    n = count(mask)
    if n < 2
        @warn "Fewer than 2 M3 particles in previous population (n = $n); using σ_u = $SIGMA_U_MIN"
        return SIGMA_U_MIN
    end
    us = u_from.(df.p2[mask], df.p0[mask])
    return max(SIGMA_U_MIN, std(us; corrected = true))   # 1/(n-1) normalization
end

# Local copy of the perturbation widths with σ_u set (avoids mutating a shared struct)
with_sigma_u(s::perturbation, σu) =
    perturbation(s.sigma1s, s.sigma1T, s.sigma2s, s.sigma2T, s.sigma3p2, σu, s.sigma3T)


function one_step_threaded(n, file, data, dataCF, nsample, N, k, L, l, LTT_threshold, CF_threshold,
                           current_step, n_models, sigma, name)
    if current_step != 1
        df = CSV.read(file, DataFrame)
        # σ_u is a deterministic function of the previous population, so every job
        # computes the same value, and it stays fixed for proposal AND weighting.
        sig = with_sigma_u(sigma, compute_sigma_u(df))
        n == 1 && println("step $current_step: σ_u = $(sig.sigma3u)")
    else
        df = []
        sig = sigma
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
                x = sample_parameter(df, model_choice, current_step, L, sig)
                LTT, CF = run(x, model_choice, current_step, data, dataCF, N, k, L, l)
                diff = piecewise_difference_g(data[1], data[2], LTT[1], LTT[2])
                dist = area(diff[1], diff[2]) / k
                CFd = abs(CF - dataCF)
                if dist < LTT_threshold && CFd < CF_threshold
                    weight = calculate_weight(x, model_choice, current_step, L, df, sig)
                    p2, p0, T_m = x.p2, x.p0, x.T_m
                    println(io, "$p2,$p0,$T_m,$model_choice,$dist,$CFd,$weight")
                    flush(io)
                    break
                end
            end
        end
    end
end


function sample_parameter(df, model_choice, current_step, L, sigma)
    x = point(0.0, 0.0, 0)
    if current_step == 1 # sample from prior
        if model_choice == 1 # VanEgeren
            x.p2 = rand() * 2
            x.T_m = round(Int, rand(Uniform(0, L - 3)))
        elseif model_choice == 2 # Williams
            x.p2 = rand() * 2
            x.T_m = rand(Uniform(0, L - 3))
        elseif model_choice == 3
            # p2 ~ U(0,1), u ~ U(0,1), p0 = u * m(p2)
            x.p2 = rand()
            u = rand()
            x.p0 = u * m_half(x.p2)
            x.T_m = rand(Uniform(0, L - 3))
        end
    else
        mask = df.model .== model_choice
        original_indices = findall(mask)
        weighted_sampler = Categorical(Float64.(df.weight[mask]))
        masked_index = rand(weighted_sampler)
        original_index = original_indices[masked_index]
        p2 = df.p2[original_index]
        p0 = df.p0[original_index]
        T_m = df.T_m[original_index]
        if model_choice == 1
            sample = -1
            while sample < 0 || sample > L - 3
                sample = round(Int, rand(Uniform(T_m - sigma.sigma1T, T_m + sigma.sigma1T)))
            end
            x.T_m = sample

            sample = -1
            while sample < 0 || sample > 2
                sample = rand(Uniform(p2 - sigma.sigma1s, p2 + sigma.sigma1s))
            end
            x.p2 = sample
        elseif model_choice == 2
            sample = -1
            while sample < 0 || sample > L - 3
                sample = rand(Normal(T_m, sigma.sigma2T))
            end
            x.T_m = sample

            sample = -1
            while sample < 0 || sample > 2
                sample = rand(Uniform(p2 - sigma.sigma2s, p2 + sigma.sigma2s))
            end
            x.p2 = sample
        elseif model_choice == 3
            #   T'  ~ N(T_j, σT²) truncated to [0, L-3]
            #   p2' ~ U(max(0, p2_j - σ2), min(1, p2_j + σ2))
            #   u'  ~ U(max(0, u_j  - σu), min(1, u_j  + σu))
            #   p0' = u' * m(p2')      (so p0' + p2' ≤ 1 automatically)
            x.T_m = rand(truncated(Normal(T_m, sigma.sigma3T), 0, L - 3))

            p2_new = rand(Uniform(max(0.0, p2 - sigma.sigma3p2), min(1.0, p2 + sigma.sigma3p2)))

            u_j = u_from(p2, p0)
            u_new = rand(Uniform(max(0.0, u_j - sigma.sigma3u), min(1.0, u_j + sigma.sigma3u)))

            x.p2 = p2_new
            x.p0 = u_new * m_half(p2_new)
        end
    end
    return x
end

function run(x, model_choice, current_step, data, dataCF, N, k, L, l)
    if model_choice == 1
        a = nothing
        while a === nothing
            a = sim_LTT(round(Int, (L - x.T_m)), N, x.p2, k, L, l)
        end
        return a[1], a[2]
    elseif model_choice == 2
        S_result = nothing
        while S_result === nothing
            S_result = run_selection_sim(N, x.T_m, L, x.p2)
        end
        tree = construct_mutated_tree(S_result, k, l)
        LTT = LTT_plot(tree, l, L, 0.0)
        CF = S_result.CF
        return LTT, CF
    elseif model_choice == 3
        p = Parameters(1 / 30, x.p0, x.p2, x.T_m * 365)
        l_result = nothing
        while l_result === nothing
            l_result = tree_simulation(p, L * 365, 0, 2000, k)
        end
        LTT = [(x.T_m .+ l_result[1][1] ./ 365) .* l, l_result[1][2]]
        return LTT, l_result[2]
    end
end

# Model 3 kernel density K(θ' | θ_j) in (p2, u, T) space
function kernel_density_m3(p2n, p0n, Tn, p2j, p0j, Tj, sigma, TL)
    a2, b2 = max(0.0, p2j - sigma.sigma3p2), min(1.0, p2j + sigma.sigma3p2)
    (a2 <= p2n <= b2) || return 0.0
    uj, un = u_from(p2j, p0j), u_from(p2n, p0n)
    au, bu = max(0.0, uj - sigma.sigma3u), min(1.0, uj + sigma.sigma3u)
    (au <= un <= bu) || return 0.0
    return (1 / (b2 - a2)) * (1 / (bu - au)) *
           pdf(truncated(Normal(Tj, sigma.sigma3T), 0, TL), Tn)
end

function calculate_weight(x, model_choice, current_step, L, df, sigma)
    if current_step == 1
        wt = 1.0
    else
        if model_choice == 1
            su = 0.0
            for i in 1:nrow(df)
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
            prior = 1 / 2 * (L - 3) / 2
            wt = prior / su

        elseif model_choice == 2
            su = 0.0
            for i in 1:nrow(df)
                if df.model[i] == 2
                    try
                        weight_contrib = df.weight[i] *
                            pdf(Truncated(Uniform(df.p2[i] - sigma.sigma2s, df.p2[i] + sigma.sigma2s), 0, 2), x.p2) *
                            pdf(Truncated(Normal(df.T_m[i], sigma.sigma2T), 0, L), x.T_m)
                        su += weight_contrib
                    catch e
                        continue
                    end
                end
            end
            prior = 1 / 2 * L / 2
            wt = prior / su

        elseif model_choice == 3
            TL = L - 3
            su = 0.0
            for i in 1:nrow(df)
                df.model[i] == 3 || continue
                su += df.weight[i] *
                      kernel_density_m3(x.p2, x.p0, x.T_m, df.p2[i], df.p0[i], df.T_m[i], sigma, TL)
            end
            # Prior density in (p2, u, T): U(0,1) × U(0,1) × U(0, L-3) = 1/(L-3)
            prior = 1 / TL
            if su <= 0
                @warn "Zero kernel mass for new M3 particle; weight set to 0"
                return 0.0
            end
            wt = prior / su
        end
    end
    return wt
end



emp_sd(df, model, col) = std(df[df.model .== model, col])

function compute_sigma(df)
    mask3 = df.model .== 3
    us = u_from.(df.p2[mask3], df.p0[mask3])
    perturbation(
        emp_sd(df, 1, :p2),  emp_sd(df, 1, :T_m),   # model 1
        emp_sd(df, 2, :p2),  emp_sd(df, 2, :T_m),   # model 2
        emp_sd(df, 3, :p2),                         # model 3: σ₂
        max(SIGMA_U_MIN, std(us)),                  # model 3: σ_u (your formula)
        emp_sd(df, 3, :T_m))                        # model 3: σ_T
end
