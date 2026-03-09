using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("smc_onestep.jl")

file_path = "data/PD7271.csv"
data,header= readdlm(file_path,',',header=true)
df = DataFrame(data,vec(header))
M=Matrix(df)
dm=pairwise(cityblock,M)
cells= Vector{String}(vec(header))

tree= neighbor_joining(dm,cells)
mutated_cells=Vector()
mut_cells_id=Vector()
for i in 1:length(cells)
    if M[20830,i]==1
        push!(mutated_cells,cells[i])
        push!(mut_cells_id,i)
    end
end
L= leave_incidence_matrix(tree)

node_remove=Vector{AbstractNode}()
for j in 1:173
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
dataCF=0.2159090909
N      = 10^5
k      = 19
L      = 23.0
l      = h/L
N_points = 2000
LTT_threshold=1000.0
CF_threshold=1.0
current_step=1
sigma=perturbation(sqrt(0.0015655725838679324),sqrt(0.9670825412706355),sqrt(0.0025681746089748097),sqrt(2.7037102594879507),0.01892285622694091,0.019759485653441468,sqrt(0.6154543444503155))
file = "PD7271_normalized_step_$(current_step -1).csv"
#n_threads = nthreads()
n_threads =250
nsample= N_points /n_threads 
println("running")
#results = Vector{Any}(undef, n_threads)
# @threads for i in 1:n_threads
#     results[i] = one_step_threaded(i,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma)
# end
name="PD7271"
task_id = parse(Int, ARGS[1])
results = one_step_threaded(task_id,file,data,dataCF,round(Int,nsample),N,k,L,l,LTT_threshold,CF_threshold,current_step,3,sigma,name)
