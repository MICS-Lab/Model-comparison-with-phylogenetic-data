using Distances, Phylo
using CSV, DataFrames
using DelimitedFiles
using MCPhyloTree
using Statistics
using Random
using Distributions 
using Glob

N_points = 600
LTT_threshold = 1000
CF_threshold = 1.0
current_step = 1
nfiles = 20

function read_previous_step(step_number, n_threads, file_pattern)
    col_names = [:p2, :p0, :T_m, :model, :dist, :CFd, :weight]
    
    df = DataFrame(p2=Float64[], p0=Float64[], T_m=Float64[], 
                   model=Int[], dist=Float64[], CFd=Float64[], weight=Float64[])
    
    for i in 1:n_threads
        filename = replace(file_pattern, "INDEX" => string(i))
        
        # Check if file exists and is not empty
        if isfile(filename) && filesize(filename) > 0
            try
                temp_df = CSV.read(filename, DataFrame, header=false)
                
                # Skip if DataFrame is empty
                if nrow(temp_df) > 0
                    rename!(temp_df, col_names)
                    append!(df, temp_df)
                else
                    @warn "File $filename is empty (no rows), skipping"
                end

            catch e
                @warn "Error reading file $filename: $e, skipping"
            end
        elseif isfile(filename)
            @warn "File $filename has size 0, skipping"
        end
    end
    return df
end

# Read first set of files
df1 = read_previous_step(current_step, 260, 
                         "PD5163_step_$(current_step)_INDEX.csv")

# Read second set of files (change the pattern to match your new files)
# df2 = read_previous_step(current_step, 25, 
#                          "/gpfs/workdir/brunom/manual_smc/PD5163_W/PD5163_step2_$(current_step)_INDEX.csv")

# Join the two dataframes
#df = vcat(df1, df2)
df=df1

# Fix infinite weights
for e in eachrow(df)
    if e.weight == Inf 
        e.weight = 1
    end
end

function normalize_weights_by_model!(df)
    """
    Normalize weights within each model so they sum to 1
    """
    for model_id in unique(df.model)
        mask = df.model .== model_id
        total_weight = sum(df.weight[mask])
        
        if total_weight > 0
            df.weight[mask] ./= total_weight
        else
            @warn "Model $model_id has total weight of 0"
        end
    end
    return df
end

df = normalize_weights_by_model!(df)
println("Total rows: ", nrow(df))

# Write to file
CSV.write("PD5163_normalized_step_$(current_step).csv", df, 
          writeheader=true,
          delim=',')

# Calculate statistics
mask = df.model .== 1
println("sigmas1 $(var(df.p2[mask]))")
println("sigmaT1 $((var(df.T_m[mask])))")

mask = df.model .== 2
println("sigmas2 $(var(df.p2[mask]))")
println("sigmaT2 $((var(df.T_m[mask])))")

mask = df.model .== 3
println("sigmas3 $(var(df.p2[mask]))")
println("sigmas3 $(var(df.p0[mask]))")
println("sigmaT3 $((var(df.T_m[mask])))")

println("distLTT $(quantile(df.dist, 0.7))")
println("distCF $(quantile(df.CFd, 0.7))")