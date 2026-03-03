using Distances, Phylo
using Plots, CSV, DataFrames
using DelimitedFiles
using MCPhyloTree
using Statistics
using Random

mutable struct Cell 
    id::Int 
    tBirth::Float64 
    niv::Int 
    parent::Union{Cell,Nothing}
    children::Vector{Cell}
end

mutable struct cellSimulation 
    totalpop::Int
    startTime::Float64
    endTime::Float64
    nodes::Dict{Int,Cell}
    currentgeneration::Vector{Cell}
    currentTime::Float64
    status::Int #stochastic or deterministic 
    p2::Float64 #probability of a cell to divide
    p0::Float64 #probability of a cell to die
    p1::Float64 #probability of a cell to stay alive
    alpha::Float64 #rate of cell division
    trajectory::Vector{Tuple}
    idx::Int
end 

mutable struct Parameters 
    alpha::Float64
    p0::Float64
    p2::Float64
    T_m::Float64
end 

function run_simulation(parameters, patientAge,status,maxpop)
    firstCell=Cell(1,parameters.T_m,0,nothing,Cell[])
    hsc= cellSimulation(1,parameters.T_m,patientAge,Dict(1 => firstCell),[firstCell],parameters.T_m,status,parameters.p2,parameters.p0,1- parameters.p2 -parameters.p0 ,parameters.alpha,[],1)
    while hsc.currentTime < hsc.endTime
        prop_sym = hsc.alpha * hsc.p2 * hsc.totalpop
        prop_asym = hsc.alpha * hsc.p1 * hsc.totalpop 
        prop_diff = hsc.alpha * hsc.p0 * hsc.totalpop
        
        #total_propensity = prop_sym  + prop_diff + prop_asym
        total_propensity = hsc.alpha * hsc.totalpop
        τ = -log(rand()) / total_propensity
        hsc.currentTime += τ
        #r = rand() * total_propensity
        r = rand()
        if r <= hsc.p2 #symmetric cell division
            SymDiv!(hsc)
        elseif r <= hsc.p2 +hsc.p1 #asymmetric cell division
            AsymDiv!(hsc)
        else #differentiation
            Differentiation!(hsc)
        end
        if length(hsc.currentgeneration) != hsc.totalpop
            #println("Population mismatch at time: $(hsc.currentTime). Expected: $(hsc.totalpop), Found: $(length(hsc.currentgeneration))")
            break
        end
        push!(hsc.trajectory, (hsc.currentTime, hsc.totalpop))
        if hsc.totalpop == 0 
            #println("Population reached zero at time: $(hsc.currentTime)")
            break 
        end
        if hsc.totalpop > maxpop
            current_hsc = hsc.totalpop
            time= hsc.currentTime
            hsc.status =1 #deterministic 
            net_rate = hsc.alpha *(hsc.p2 - hsc.p0)
            remaining_time = patientAge - time
            if net_rate != 0
                # Add intermediate points for smooth trajectory
                n_points = min(100, Int(ceil(remaining_time / 0.1)))  # Point every 0.1 time units
                for i in 1:n_points
                    t_point = time + (i / n_points) * remaining_time
                    if net_rate > 0
                        # Exponential growth: N(t) = N₀ * exp(r * t)
                        n_point = current_hsc * exp(net_rate * (t_point - time))
                        
                        #n_point = 10^5 / (1 + ((10^5)/current_hsc - 1) * exp(-net_rate * (t_point - time)))
                        
                    else
                        # Exponential decay
                        n_point = max(0, current_hsc * exp(net_rate * (t_point - time)))
                    end
                    push!(hsc.trajectory, (t_point, n_point))
                end
            else
                # Constant population (net_rate = 0)
                #push!(trajectory, (time=t_final, state=Dict("HSC" => current_hsc)))
            end
            #println("Switch to det approximation: $(hsc.totalpop)")
            break
        end
    end 
    #println("Simulation ended at time: $(t_point) with population: $(n_point)")
    return hsc
end

function SymDiv!(hsc::cellSimulation)
    random_index = rand(1:length(hsc.currentgeneration))
    node = hsc.currentgeneration[random_index]
    child1= Cell(hsc.idx+1, hsc.currentTime, node.niv + 1, node, Cell[])
    child2= Cell(hsc.idx+2, hsc.currentTime, node.niv + 1, node, Cell[])
    hsc.nodes[child1.id] = child1
    hsc.nodes[child2.id] = child2
    push!(node.children, child1)
    push!(node.children, child2)
    hsc.currentgeneration[random_index]=child1
    #push!(hsc.nodes, child1)
    push!(hsc.currentgeneration, child2)
    hsc.totalpop += 1    
    hsc.idx += 2
end

function AsymDiv!(hsc::cellSimulation)
    random_index = rand(1:length(hsc.currentgeneration))
    node = hsc.currentgeneration[random_index]
    child= Cell(hsc.idx +1, hsc.currentTime,node.niv+1, node, Cell[])
    hsc.nodes[child.id] = child
    push!(node.children, child)
    hsc.currentgeneration[random_index]=child
    hsc.idx +=1
end

function Differentiation!(hsc::cellSimulation)
    random_index = rand(1:length(hsc.currentgeneration))
    node = hsc.currentgeneration[random_index]
    
    # Always remove from current generation and update population
    hsc.totalpop -= 1
    deleteat!(hsc.currentgeneration, random_index)
    
    # Remove the node and clean up the tree
    remove_node_and_cleanup!(hsc, node)
end

function remove_node_and_cleanup!(hsc::cellSimulation, node::Cell)
    # Remove the node itself
    delete!(hsc.nodes, node.id)
    
    # Remove from parent's children if it has a parent
    if node.parent !== nothing
        node.parent.children = filter(c -> c.id != node.id, node.parent.children)
        
        # If parent now has no children, remove it too
        if isempty(node.parent.children)
            remove_node_and_cleanup!(hsc, node.parent)
        end
    end
end

function to_newick_with_node_age(root::Cell)
    # Recursive helper function to construct the Newick string
    function recurse_tree(node::Cell, parent_height)
        children = node.children
        current_height = node.tBirth
        #if parent_height < current_height 
            branch_length =current_height - parent_height
        #else
            #branch_length=0
        #end
        
        if !isempty(children)
            # Internal node: recurse on all children
            subtree = join([recurse_tree(child, current_height) for child in children], ",")
            return "($subtree)$(node.id):$branch_length"
        else
            # Leaf node: return the node's ID with the correct branch length
            return "$(node.id):$branch_length"
        end
    end
    return recurse_tree(root, 0) * ";"
end

function subsample_mutated_tree(sim,tree,n)
    leaves=getleaves(tree)
    mut_leaves=[el for el in leaves]
    #Random.seed!(11)
    s=sample(mut_leaves, n, replace=false)
    anc=mrca(tree,s)
    r=Phylo.getroot(tree)
    routes=[noderoute(tree, r, s[i]) for i in 1:n] 
    #push!(routes, noderoute(tree,r,anc))
    for el in nodeiter(tree)
        if any(el in route for route in routes)==false && el!=r 
            deletenode!(tree,el)
        end
    end
    return tree
end
function tree_simulation(parameters, patientAge,status,maxpop,k)
    hsc=run_simulation(parameters, patientAge,status,maxpop)
    if length(hsc.currentgeneration) >k 
        Δ= hsc.p2 -hsc.p0 
        final_pop = hsc.trajectory[end][2]
        N=10^5
        #CF=(1 - Δ) * final_pop / ((1 - Δ) * final_pop + N)
        CF=final_pop/(final_pop+N)
        n=collect(values(hsc.nodes))
        root=n[findfirst(x -> x.parent === nothing,n)]
        nw_string=to_newick_with_node_age(root)
        tree=parsenewick(nw_string)
        t=subsample_mutated_tree(hsc,tree,k)
        LTT=LTT_plot(t, 1.0, 60.0, 0.0)
        return LTT,CF
    else 
        return nothing 
    end 
end