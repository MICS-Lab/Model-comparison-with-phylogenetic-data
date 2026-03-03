using Plots, Distributions,Statistics 
using MCPhyloTree
using DataFrames, DelimitedFiles
using CSV 
using Distances
using Random
using StatsBase


"""
function traj(g, N, s, k) to calculate trajectories (number of mutated cells at each generation) 
inputs: g=number of generations 
        N=maximum number of cells 
        s=fitness of mutation 
        k=number of mutated cells for which tree has to be reconstructed 
output: n= Vector{Int} where n[i]= number of cells at generation i 
"""
function traj(g, N, s, k)
    while true
        n = Vector{Int}(undef, g)
        n[1] = 1
        for i in 1:(g-1)
            p = ((1 + s) * n[i]) / (N + n[i] * s)
            if p < 0
                n[i+1] = 0
            elseif p >1 
                p=1
                n[i+1] = rand(Binomial(N, p))
            else
                n[i+1] = rand(Binomial(N, p))
            end
        end
        #function is recursive until the number of cells in the last generation is at least k (no extinction)
        if n[end] > k
            return n
        end
    end
end
"""
function generations(n,g) to give a name (integer number) to each cell at every generation
inputs: n=Vector{Int} where n[i]= number of cells at generation i
        g=number of generations
output: nodes=Vector{Vector{Int64}}() where nodes[i]=Vector containing the names of cells at generation i 
"""
function generations(n,g)
    nodes=Vector{Vector{Int64}}(undef, g)
    q::Int64=0
    for i in 1:g
        v = Vector{Int64}(undef, n[i])
        for j::Int64 in 1:n[i]
            v[j]=q
            q+=1
        end
        nodes[i] = v 
    end
    return nodes
end

"""
function child_parent(leaves, parents) that builds a dictionary with child-parent relationships at random
inputs: leaves=Vector of names of children cells 
        parents=Vector of names of parent cells 
output: dict = Dict{Int, Int}() where [child] => [parent]
"""
function child_parent(leaves::Vector, parents::Vector)
    dict = Dict{Int, Int}()  
    for leaf in leaves
            x=rand(parents) 
            dict[leaf]=x
    end
    return dict
end

"""
function build_tree(relations) that builds nodes and branches based on child-parent relationships 
inputs: relations=Vector{Dict{Int,Int}} where relations[i]= dictionary [child] => [parent] from generation i-1
output: root=GeneralNode root of the tree built. 
"""

function build_tree(relations::Vector{Dict{Int,Int}})
    root=MCPhyloTree.Node("0")
    child_id=MCPhyloTree.Node()
    dict_tree=Dict{Int,GeneralNode}()
    #dictionary of all nodes of the tree where keys are the names of the cells 
    dict_tree[0]=root
    dict_tree[0].root= true
    dict_tree[0].mother= missing
    MCPhyloTree.initialize_tree!(dict_tree[0])
    for j in 1:length(relations) 
        for (child,parent) in pairs(relations[j])
            dict_tree[0].root= true
            child_id=MCPhyloTree.Node(string(child))
            dict_tree[child]=child_id
            #child node created and added to dictionery of nodes
            add_child!(dict_tree[parent],dict_tree[child])
            #parent node updated with new child 
        end
        MCPhyloTree.update_tree!(dict_tree[0])
    end
    MCPhyloTree.update_tree!(dict_tree[0])
    #print(dict_tree[0].n_desc)
    return root
end

"""
function branch_poisson(r,l) to change the length of branches according to Poisson distribution
inputs: r=root of the tree 
        l= parameter of the Poisson (mutation rate)
output: r= root of the tree with length of branches changed
"""

function branch_poisson(r,l)
    p=Poisson(l)
    for node in level_order(r)
        node.inc_length=rand(p)
        #for each node the length of the incident branch is changed 
    end
    return r
end

"""
function LTT_plot(t,L,g,k) to create LTT plot for a tree 
inputs: t=root of the tree
        L=total number of generations (age of the patient)
        g=number of generations for the mutated lineage
outputs: height= vector containing the height of each node
        mutated= vector containing the number of lineages at each node 
"""
function LTT_plot(t,L,g,k,l)
    x= post_order(t)
    height=Vector(undef,length(x)+1)
    mutated=Vector()
    h=node_age(t)+(L-g)*l
    height[1]=0
    push!(mutated,1)
    m=1
    he=Vector(undef,length(x))
    for i in 1:length(x)
        he[i]=h-node_age(x[i])
    end
    index= sortperm(he)
    x= x[index]
    height[2:end]=he[index]
    j=0
    for i in 1:length(x)
        if x[i].nchild!=0
            push!(mutated,m)
            insert!(height,i+j+1,h-node_age(x[i]))
            m=m+x[i].nchild-1
            j+=1
        end
        #j+=1
        #if h-node_age(x[i])> 1202.5868421852424
            #height[i+j]=1202.5868421852424
        #end
        push!(mutated,m)
    end
    #height[end]=1202.5868421852424
    push!(height,L*l)
    push!(mutated,m)
    return [height,mutated]
end

"""
function sim_LTT(g,N,s,k,L,l) to simulate a tree and obtain a LTT plot 
inputs: g=number of generations for the mutated lineage
        N=maximum number of cells 
        s=fitness of mutation 
        k=number of mutated cells for which tree has to be reconstructed 
        L=total number of generations (age of the patient)
        l= mutation rate of cells 
output: height= vector containing the height of each node
        mutated= vector containing the number of lineages at each node 
"""
function sim_LTT(g,N,s,k,L,l)
    x1=traj(g,N,s,k)
    b=generations(x1,g)
    #print("gen")
    lin=Vector{Dict{Int,Int}}(undef,g-1)
    last_gen=sample(b[g],k, replace=false)
    #random sampling of k cells from last generation
    lin[g-1]=child_parent(last_gen,b[g-1])
    for i in 2:g-1
        lin[g-i]=child_parent(collect(values(lin[g-i+1])),b[g-i])
        #parents chosen at random at this generation are children of the previous generation
    end
    #print("lin")
    tree=build_tree(lin)

    tree1=branch_poisson(tree,l)
    #print(tree)
    #Average(tree1)
    #LTT=LTT_plot(tree1,l,L,0.0)
    LTT=LTT_plot(tree1,L,g,k,l)
    #print("LTT")
    return [LTT,x1[end]/(N)]
end

""" 
function piecewise_difference_g(x1, y1, x2, y2) to calculate the absolute value of the distance between two piecewise constant curves
inputs: x1=vector of x values for curve 1
        y1=vector of y values for curve 1
        x2=vector of x values for curve 2
        y2=vector of y values for curve 2
outputs: difference_x=vector of x values for absolute value of distance curve 
         difference_y=vector of y values for absolute velue of distance curve
"""
function piecewise_difference_g(x1, y1, x2, y2)
    x_union = sort(union(x1, x2))
    difference_x = Vector()
    difference_y = Vector()
    prev_diff = abs(y1[1] - y2[1])
    prev_x = x_union[1]
    
    for x in x_union[2:end]
        y1_val = get_value(x, x1, y1)
        y2_val = get_value(x, x2, y2)
        diff = abs(y1_val - y2_val)
        
        push!(difference_x, prev_x)
        push!(difference_y, prev_diff)
        
        #if diff != prev_diff
            push!(difference_x, x)
            push!(difference_y, prev_diff)
            prev_diff = diff
            prev_x = x
        #end
    end
    
    push!(difference_x, prev_x)
    push!(difference_y, prev_diff)
    
    return [difference_x, difference_y]
end

"""
function get_value(x, x_values, y_values) to return the y value corresponding to the last occurence of a x value in a vector
inputs: x=value to found the corresponding y
        x_values=vector of x values 
        y_values=vector of y values 
outputs: y_values[i]=y value 

"""
function get_value(x, x_values, y_values)
    idx = searchsortedlast(x_values, x)
    if idx == 0
        return y_values[1]
    elseif idx > length(x_values)
        return y_values[end]
    else
        return y_values[idx]
    end
end

"""
function area(x,y) to calculate the area under a LTT plot (just made of squares and triangles)
inputs: x=vector cordinates x for each point of the plot
        y=vector of cordinates y for each point of the plot 
output: A= Float with area under the graph 
"""
function area(x,y)
    A::Float64=0
    for i in 1:length(x)-1
        if y[i+1]!=y[i]
            A=A+(x[i+1]-x[i])*(y[i+1]-y[i])/2
            if y[i]!=0
                A=A+(x[i+1]-x[i])*y[i]
            end
            #area of the trapeziod
        else 
            A=A+(x[i+1]-x[i])*y[i]
            #area of the rectangle
        end
    end
    return A 
end