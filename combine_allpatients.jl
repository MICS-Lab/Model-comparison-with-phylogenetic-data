using CSV, DataFrames, Statistics, Random, Distributions
using Distances


# ---------------- PARAMETERS ----------------
current_step = 1
nfiles = 305
npatients =4

# ---------------- FILE READER ----------------
function read_previous_step(step_number, n_threads, file_pattern)
    col_names = [:p2, :p0, :T_m, :model, :dist, :CFd, :weight, :name]

    df = DataFrame(p2=Float64[], p0=Float64[], T_m=Float64[],
                   model=Int[], dist=Float64[], CFd=Float64[],
                   weight=Float64[], name=String[])

    for i in 1:n_threads
        filename = replace(file_pattern, "INDEX" => string(i))

        if isfile(filename) && filesize(filename) > 0
            try
                temp_df = CSV.read(filename, DataFrame, header=true)
                if nrow(temp_df) > 0
                    rename!(temp_df, col_names)
                    append!(df, temp_df)
                else
                    @warn "File $filename empty"
                end
            catch e
                @warn "Error reading $filename: $e"
            end
        end
    end
    return df
end


df = read_previous_step(current_step, nfiles,
                         "all_patients_step$(current_step)_INDEX.csv")

println("Total merged rows: ", nrow(df))

# ---------------- FIX INFINITE WEIGHTS ----------------
df.weight .= ifelse.(isinf.(df.weight), 1.0, df.weight)

# ---------------- NORMALIZE WEIGHTS PER PATIENT & MODEL ----------------
function normalize_weights_by_patient_model!(df)
    group_cols = [:name, :model]
    grouped = groupby(df, group_cols)

    for g in grouped
        total_w = sum(g.weight)
        if total_w > 0
            g.weight ./= total_w
        else
            @warn "Zero total weight for patient=$(g.name[1]) model=$(g.model[1])"
        end
    end
    return df
end

normalize_weights_by_patient_model!(df)

# ---------------- SAVE MERGED NORMALIZED FILE ----------------
out_file = "ALL_normalized_step_$(current_step).csv"
CSV.write(out_file, df)
println("Saved → $out_file")

# ---------------- SUMMARY STATISTICS ----------------
println("\n===== PATIENT-SPECIFIC PARAMETER VARIANCES =====")
for pname in unique(df.name)
    dfi = df[df.name .== pname, :]
    println("\nPatient $pname:")
    for m in unique(dfi.model)
        mask = dfi.model .== m
        println("Model $m:")
        println("  var(p2)  = ", var(dfi.p2[mask]))
        println("  var(p0)  = ", var(dfi.p0[mask]))
        println("  var(T_m) = ", var(dfi.T_m[mask]))
    end
    println("LTT 80% quantile = ", quantile(dfi.dist, 0.8))
    println("CF  80% quantile = ", quantile(dfi.CFd, 0.8))
end

println("\n===== MODEL STATISTICS =====")
for m in unique(df.model)
    mask = df.model .== m
    println("Model $m:")
    prob= count(x -> x==m, df.model)./(length(df.model))
    println(prob)
end

# println("\n===== MODEL-SPECIFIC PARAMETER VARIANCES =====")
# for m in unique(df.model)
#     mask = df.model .== m
#     println("Model $m:")
#     println("  var(p2)  = ", var(df.p2[mask]))
#     println("  var(p0)  = ", var(df.p0[mask]))
#     println("  var(T_m) = ", var(df.T_m[mask]))
# end

# println("\n===== DISTANCE QUANTILES (GLOBAL) =====")
# println("LTT 80% quantile = ", quantile(df.dist, 0.8))
# println("CF  80% quantile = ", quantile(df.CFd, 0.8))

# # ---------------- PATIENT-SPECIFIC SUMMARIES ----------------
# println("\n===== PATIENT-SPECIFIC POSTERIOR SUMMARIES =====")

# for pname in unique(df.name)
#     dfi = df[df.name .== pname, :]
#     println("\nPatient $pname:")
#     println("  mean p2  = ", mean(dfi.p2))
#     println("  mean p0  = ", mean(dfi.p0))
#     println("  mean T_m = ", mean(dfi.T_m))
#     println("  ESS ≈ ", 1 / sum(dfi.weight .^ 2))
# end
