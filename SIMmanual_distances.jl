using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("/gpfs/workdir/brunom/manual_smc/distance_comparison/smc_onestep_distances.jl")


df=CSV.read("/gpfs/workdir/brunom/manual_smc/distance_comparison/LTTW.csv", DataFrame)
data= [df.Mutations , df.Linages]

dataCF=0.6867665841781752
N      = 10^5
k      = 15
L      = 53.2
l      = 19.0
N_points = 1000
LTT_threshold=140.77333333333328
CF_threshold=0.18630687909019067
current_step=13
sigma=perturbation(sqrt(0.0014851256044979579),sqrt(14.619284642321155),0.04490213900695595,sqrt(23.028768546369594),0.02043799067406225,0.020832380491794887,sqrt(100.56094656813376))
file = "manual_smc/distance_comparison/SIMWcombined_normalized_step_$(current_step -1).csv"
#n_threads = nthreads()
n_threads =250
nsample= N_points /n_threads 
println("running")
#results = Vector{Any}(undef, n_threads)
# @threads for i in 1:n_threads
#     results[i] = one_step_threaded(i,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma)
# end
name="combined"
task_id = parse(Int, ARGS[1])
results = one_step_threaded(task_id,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,1,sigma,name)