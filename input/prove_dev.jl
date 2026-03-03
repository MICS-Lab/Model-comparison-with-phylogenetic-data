using SpecialFunctions  # for factorial-related functions
using Random  # for random number generation
using Roots  # for root-finding

# First n terms of the exponential as expressed as a proportion of the full exponential
function expn(n, x)
    sum(exp.(.-x + i * log.(x) .- logfactorial(i)) for i in 0:n)
end

# Cumulative density function with flat prior interval from 0 to T₀
function CDF(x, m=5, λ=20, T₀=50)
    denominator = 1 - expn(m, λ * T₀)
    if isapprox(denominator, 0.0; atol=1e-10)
        error("Denominator in CDF computation is too small, leading to division by zero.")
    end
    [(1 - expn(m, λ * xi)) / denominator for xi in x]
end

# Cumulative density function in λ-space, converting back to time at the end
function CDF2(x, m=5, LAMBDA=50 * 20)
    denominator = 1 .- expn(m, LAMBDA)
    if isapprox(denominator, 0.0; atol=1e-10)
        error("Denominator in CDF2 computation is too small, leading to division by zero.")
    end
    [(1 .- expn(m, xi)) / denominator for xi in x]
end

# Logistic mean function for the standard development model
function logisticMean(L, k, midpoint, a, b)
    ((L / k) * (log(1 + exp(k * (b - midpoint))) - log(1 + exp(k * (a - midpoint))))) / (b - a)
end

# Total λ function based on time intervals
function totalLambda(λ, t, t₀=0)
    (t - t₀) * λ + (t - t₀) * logisticMean(149.2, -50.0, 0.22,t₀, t)
    #(t - t₀) * λ + (t - t₀) * 149.2/(1+exp(50.0*(t - t₀ - 0.225)))
end
# Simulation function to sample from a CDF
function simfromcdf(cdf, start, mutend, size)
    # Granularity of 1/1e4 of the range
    x = range(start, mutend, length=10000)
    cval = cdf(x)
    
    # Remove values from start of CDF where the function is near-zero
    idx = findfirst(x -> x > 1e-8, cval)
    if idx !== nothing && idx > 1
        cval = cval[idx:end]
        x = x[idx:end]
    end

    # Generate random uniform values and find their intervals
    rval = rand(size)
    idxs= map(r -> min(searchsortedfirst(cval, r), length(x) - 1), rval)
    return 0.5 * (x[idxs] .+ x[idxs .+ 1])
end

# Sampling acquisition times
function sample_acq_time(m, λ, T₀, N, t₀=0)
    m = round(Int, m)
    clambda = simfromcdf(x -> CDF2([x], m, λ * (T₀ - t₀))[1], 0, λ * (T₀ - t₀), N)
    
    [find_zero(x -> totalLambda(λ, x, t₀) - clambda_val, (t₀ + 1e-6, T₀)) for clambda_val in clambda]
end

function pdf_from_cdf2(x::Float64, m::Int, λ::Float64, T₀::Float64, t₀::Float64=0.0; epsilon::Float64=1e-6)
    # Define the CDF function
    cdf_function = y -> CDF2([y], m, λ * (T₀ - t₀))[1]
    
    # Compute numerical derivative of the CDF to approximate the PDF
    pdf_value = (cdf_function(x + epsilon) - cdf_function(x - epsilon)) / (2 * epsilon)
    return max(0, pdf_value)  # Ensure PDF is non-negative due to numerical errors
end

function uniform_with_poisson_bounds_pdf(x::Float64, λL::Float64, λU::Float64)
    # Initialize Poisson distributions for L and U
    poisson_L = Poisson(λL)
    poisson_U = Poisson(λU)
    
    # Compute PDF by summing over possible L and U values
    pdf_value = 0.0
    for l in 0:600  # Reasonable limit for Poisson bounds
        for u in (l+1):600  # Ensure u > l
            if l <= x <= u
                # Compute probabilities
                prob_L = pdf(poisson_L, l)
                prob_U = pdf(poisson_U, u)
                uniform_pdf = 1 / (u - l)  # PDF of uniform distribution
                pdf_value += prob_L * prob_U * uniform_pdf
            end
        end
    end
    
    return pdf_value
end
