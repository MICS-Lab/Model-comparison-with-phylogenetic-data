using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("smc_onestep.jl")


# df=CSV.read("/gpfs/workdir/brunom/SMC_newmodel/LTTGH3.csv", DataFrame)
# data= [df.Mutations .*19, df.Linages]
path = "data/snv_patient2.csv"
dat, header = readdlm(path, ';', header=true)
df = DataFrame(dat, vec(header))
M= Matrix(df)
x = Matrix{Float64}(M[1:end,8:41])
cells= Vector{String}(header[8:41])
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
mask= findfirst(x -> x==22, mut)
height_m=height[mask:end]
mutated_m=mut[mask:end]
pushfirst!(mutated_m,22)
pushfirst!(height_m,0.0)
data=[height_m,mutated_m .- 21]
j=findfirst(x== 2 for x in data[2])
Tmax=height_m[j] 
Tmin=height_m[j-2]
dataCF=0.2924528302
N      = 10^5
k      = 13
L      = 63.0
l      = h/L
N_points = 2000
LTT_threshold= 1000.0
CF_threshold=  1.0
current_step=1
sigma=perturbation(sqrt( 0.00025637314111458113),sqrt(14.945000250125064),sqrt(0.0005661795780182818),sqrt(18.562703695062464),0.022925129378516713,0.023481790751736727,sqrt( 18.94125864450353))
file = "ET2_normalized_step_$(current_step -1).csv"
#n_threads = nthreads()
n_threads =250
nsample= N_points /n_threads 
println("running")
#results = Vector{Any}(undef, n_threads)
# @threads for i in 1:n_threads
#     results[i] = one_step_threaded(i,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma)
# end
name="ET2"
task_id = parse(Int, ARGS[1])
results = one_step_threaded(task_id,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,1,sigma,name)