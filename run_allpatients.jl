using Distances, Phylo
    using CSV, DataFrames
    using DelimitedFiles
    using MCPhyloTree
    using Statistics
    using Random
    using Distributions  
using Base.Threads

include("smc_onestep_distances.jl")
include("smc_onsestep_allpatients.jl")

mutable struct patient
    name::String
    dataLTT::Vector{Any}
    dataCF::Float64
    k:: Int
    L::Float64
    l::Float64 
    LTT_threshold::Float64
    CF_threshold::Float64
    sigma::perturbation
end

patient_data= Vector{patient}()

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

push!(patient_data,patient("ET2",data,0.2924528302,13,63.0,h/63.0,67.24000641074481,0.23820230187999994,perturbation(sqrt(0.001643877433970934),sqrt(21.448242564696788),sqrt(0.001988933814458091),sqrt(29.599241148980106),sqrt( 2.9950877430014804e-5),sqrt(2.0739300080941974e-5),sqrt(16.512809057499606))))

path = "data/snv_patient1.csv"
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

push!(patient_data,patient("ET1",data,0.212195122,22,34.0,h/34.0,34.88999696440589,0.15546507320000003,perturbation(sqrt(0.008942605496302358),sqrt(5.586333148605302),sqrt(0.002982488997571501),sqrt(21.895172152435034),0.014198428650927778,0.014720154409003037,sqrt(6.812729318810887))))

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

push!(patient_data,patient("PD7271",data,0.2159090909,19,23.0,h/23.0,33.86080915514873,0.09273425639999999,perturbation(0.0036741332695204446,sqrt(1.3520738240381551),sqrt(0.007327567487932949),sqrt(6.422829614440439),sqrt(5.272540337987596e-6),sqrt(4.125216793550615e-6),sqrt(0.379735107227027))))

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

push!(patient_data,patient("PD5163",data,0.08571428571,6,37.0,h/37.0,47.056476683633626,0.29507571429000007,perturbation(0.00559410837441695,sqrt(6.901595042571415),sqrt(0.004727797449326753),sqrt(8.292500223736122),0.0009794862166547344,0.0010383183163624743,sqrt(5.785017924811733))))
n_threads = 300
nsamples=6 #must be a multiple of n_models
current_step=1
n_models=3
npatients=4
N=10^5
file="ALL_normalized_step_$(current_step -1).csv"
# @time one_step_allpatients(1,file,nsamples,patient_data, current_step, n_models,npatients,N)
# @threads for i in 1:n_threads
#     task_id = i
#     one_step_allpatients(i,file,nsamples,patient_data, current_step, n_models,npatients,N)
# end
task_id = parse(Int, ARGS[1])
results = one_step_allpatients(task_id,file,nsamples,patient_data, current_step, n_models,npatients,N)
