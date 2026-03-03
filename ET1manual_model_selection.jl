using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("/gpfs/workdir/brunom/manual_smc/smc_onestep.jl")


path = "/gpfs/workdir/brunom/snv_patient1.csv"
dat, header = readdlm(path, ';', header=true)
df = DataFrame(dat, vec(header))
M= Matrix(df)
x = Matrix{Float64}(M[1:end,8:48])
cells= Vector{String}(header[8:48])
dm=pairwise(cityblock,x)

tree= neighbor_joining(dm,cells)
NW=newick(tree)
t=parsenewick(NW)
#obtaining LTT plot for the data tree 
nh = nodeheights(t,noleaves=true) 
lh=nodeheights(t,onlyleaves=true) 
h=mean(lh)
index= sortperm(nh)
nh=nh[index]
c=nh.axes[1]
mut=[]
m=1
height=[]
j=0
i=1
for el in nh
    push!(mut,m)
    push!(height,el)
    global m+=length(getchildren(t,c[i]))-1
    push!(mut,m)
    push!(height,el)
    global i+=1
end
push!(mut,m)
push!(height,h)
mask= findfirst(x -> x==20, mut)
height_m=height[mask:end]
mutated_m=mut[mask:end]
pushfirst!(mutated_m,20)
pushfirst!(height_m,0.0)
data=[height_m,mutated_m .- 19]
j=findfirst(x== 2 for x in data[2])
Tmax=height_m[j] 
Tmin=height_m[j-2]
dataCF=0.212195122
N      = 10^5
k      = 22
L      = 34.0
l      = h/L
N_points = 2000
LTT_threshold= 78.04964080909092
CF_threshold=0.581398519282463
current_step=14
sigma=perturbation(sqrt(0.0012534401610438848),sqrt(3.0911815907953994),sqrt(0.003255397152288718),sqrt(9.689584472715826),0.0211972167710249,0.022160898772395247,sqrt(5.931073871760024))
file = "manual_smc/ET1_W/ET1W_normalized_step_$(current_step-1).csv"
n_threads = 250
nsample= N_points /n_threads 
name="ET1"
println("running")
# results = Vector{Any}(undef, n_threads)
# @threads for i in 1:n_threads
#     results[i] = one_step_threaded(i,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma,name)
# end
task_id = parse(Int, ARGS[1])
results = one_step_threaded(task_id,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma,name)