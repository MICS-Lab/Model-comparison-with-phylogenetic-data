using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("/smc_onestep.jl")


path = "data/PD5163.csv"
data,header= readdlm(path,',',header=true)
df = DataFrame(data,vec(header))
M=Matrix(df)
dm=pairwise(cityblock,M)
cells= Vector{String}(vec(header))

tree= neighbor_joining(dm,cells)
mutated_cells=Vector()
mut_cells_id=Vector()
for i in 1:length(cells)
    if M[24575,i]==1
        push!(mutated_cells,cells[i])
        push!(mut_cells_id,i)
    end
end
L= leave_incidence_matrix(tree)

node_remove=Vector{AbstractNode}()
for j in 1:137
    rem=0
    for i in mut_cells_id
        if L[i,j]==0 && !(j in mut_cells_id)
            rem+=1
        end
    end
    if rem==length(mut_cells_id)
        push!(node_remove,find_num(tree, j))
    end
end

node_remove=unique(node_remove)
prune_tree!(tree,node_remove)
NW=newick(tree)
t=parsenewick(NW)

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
data=[height,mut]
j=findfirst(x== 2 for x in mut)
L=37.0
l=h/L
dataCF=0.08571428571
N      = 10^5
k      = 6
N_points = 2000
LTT_threshold= 1000.0
CF_threshold=1.0
current_step=1
sigma=perturbation(sqrt(0.009248382416246365),sqrt(5.540760130065033),sqrt(0.0011367290715874044),sqrt(6.650424990910675),0.015513744228198664,0.016026444395609983,sqrt(6.096229408931873))
file = "PD5163_normalized_step_$(current_step -1).csv"
#n_threads = nthreads()
n_threads =250
nsample= N_points /n_threads 
println("running")
#results = Vector{Any}(undef, n_threads)
# @threads for i in 1:n_threads
#     results[i] = one_step_threaded(i,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma)
# end
name="PD5163"
task_id = parse(Int, ARGS[1])
results = one_step_threaded(task_id,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma,name)