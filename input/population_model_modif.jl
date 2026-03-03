using Plots, Distributions,Statistics 
using DataFrames 
using Random 
using StatsBase
using Serialization
using Phylo 
 
mutable struct Node
    tBirth::Float64
    ndiv::Int
    id::Int
    parent::Union{Nothing, Node}
    children::Vector{Node}
    driver::Bool
end
mutable struct CellCompartment
    nodesID::Vector{Int}
    ndiv::Int
    DivisionRate::Float64
    totalpop::Int
    active::Bool
end
mutable struct CellSimulation
    totalpop::Int
    cellCompartments::Vector{CellCompartment}
    ncompart::Int
    compartfitness::Vector{Float64}
    startTime::Float64
    driverAcquisitionRate::Float64
    nodes::Dict{Int,Node}
    currentTime::Float64
    lastSnap::Float64
    populationTrace::Vector{Tuple}
    ndrivers::Int
    T_mut::Float64
    atEquilibrium::Bool
    TotalDeathRate::Float64
    TotalDivisionRate::Float64
    idx::Int
    status::Int
    CF::Float64
end

"""
    addDriver(sim,row_idx,fitness,time)


function that directly modify the simulation structure and the cellCompartment structure by inserting the driver mutation 
    A new compartment (2) become active and the selected cell is pushed from the first compartment to the second, with a new fitness
"""
function addDriver(sim,row_idx,fitness)
    parent_id=sim.cellCompartments[1].nodesID[row_idx]
    if sim.nodes[parent_id].parent == nothing 
        return nothing 
    else
        #parent_id=sim.cellCompartments[1].nodesID[row_idx]
        child1=NewNode(sim,2)
        child2=NewNode(sim,1)
        addChild(sim,parent_id,child1.id)
        addChild(sim,parent_id,child2.id)
        sim.nodes[child1.id].driver=true
        #sim.nodes[child2.id].driver=true
        push!(sim.cellCompartments[2].nodesID,child1.id)
        sim.cellCompartments[1].nodesID[row_idx]=child2.id
        #push!(sim.cellCompartments[2].nodesID,child2.id)
        sim.totalpop+=1
        sim.cellCompartments[2].active=1
        #push!(sim.cellCompartments[2].nodesID,node_id)
        #deleteat!(sim.cellCompartments[1].nodesID,row_idx)
        #sim.cellCompartments[1].totalpop-=1
        sim.cellCompartments[2].totalpop=1
        sim.compartfitness[2]=fitness
        sim.ndrivers+=1
        sim.T_mut=sim.currentTime
        removeChild(sim,parent_id,child2.id)
        #addChild(sim,grandparent_id,sibling_id)
        deleteat!(sim.cellCompartments[1].nodesID,row_idx)
        sim.cellCompartments[1].totalpop-=1
        removeChild(sim,parent_id,child2.id)
        delete!(sim.nodes,child2.id)
        # sim.nodes[child1.id].parent=sim.nodes[parent_id].parent
        # grandparent_id=sim.nodes[parent_id].parent.id
        # removeChild(sim,grandparent_id,parent_id)
        # addChild(sim,grandparent_id,child1.id)
        # #deleteat!(sim.cellCompartments[comp_idx].nodesID,cell_idx)
        #                 #delete!(sim.nodes,node_id)
        # delete!(sim.nodes,parent_id)
        # sim.totalpop-=1
        #sim.cellCompartments[comp_idx].totalpop-=1
        return child1.id
    end
end

function removeChild(sim,parent_id,child_id)
    child=sim.nodes[child_id]
        filter!(x -> x!=child,sim.nodes[parent_id].children)
end
function addChild(sim,parent_id,child_id)
    push!(sim.nodes[parent_id].children,sim.nodes[child_id])
    sim.nodes[child_id].parent=sim.nodes[parent_id]
end
function getSibling(sim::CellSimulation,node_id)
    children=filter(x -> x!=sim.nodes[node_id],sim.nodes[node_id].parent.children)
    return children[1].id
end
"""
    die(sim::CellSimulation,comp_idx::Int,cell_idx::Int)

The selected cell and its parent die so are removed from the list of nodes, the other sibling is assigned as child of the granparent, its tBirth is the same as its parent but the number of symmetric divisions increases
then the compartment population and total population are updated
"""
function die(sim::CellSimulation,comp_idx::Int,cell_idx::Int)
    node_id=sim.cellCompartments[comp_idx].nodesID[cell_idx]
    parent=sim.nodes[node_id].parent
    if parent!==nothing #|| sim.nodes[node_id].tBirth == sim.T_mut #to not kill the root cell and the first driver cell
        parent_id=sim.nodes[node_id].parent.id
        if haskey(sim.nodes, parent_id) && length(sim.nodes[parent_id].children)==2 
            sibling_id=getSibling(sim,node_id)
            sim.nodes[sibling_id].tBirth=sim.nodes[parent_id].tBirth
            sim.nodes[sibling_id].ndiv+=sim.nodes[parent_id].ndiv
            sim.nodes[sibling_id].parent=sim.nodes[parent_id].parent
            #sim.nodes[sibling_id].driver==sim.nodes[parent_id].driver
            if sim.nodes[parent.id].parent !== nothing
                grandparent_id=sim.nodes[parent_id].parent.id
                #if haskey(sim.nodes, grandparent_id)
                    #grandparent_id=sim.nodes[parent_id].parent.id
                    removeChild(sim,grandparent_id,parent_id)
                    addChild(sim,grandparent_id,sibling_id)
                    deleteat!(sim.cellCompartments[comp_idx].nodesID,cell_idx)
                    delete!(sim.nodes,node_id)
                    delete!(sim.nodes,parent_id)
                    sim.totalpop-=1
                    sim.cellCompartments[comp_idx].totalpop-=1
                #end
            end
        end
    end
end
"""
    divide(sim,comp_idx,cell_idx)

function to obtain asymmetric division, two new nodes are born and they are the children of the selected node. Their tBirth is the current time and the number of symmetric divisions is equal to 1. 
    then the compartment population and the total population are updated. 
"""
function divide(sim,comp_idx,cell_idx)
    parent_id=sim.cellCompartments[comp_idx].nodesID[cell_idx]
    child1=NewNode(sim,comp_idx)
    child2=NewNode(sim,comp_idx)
    addChild(sim,parent_id,child1.id)
    addChild(sim,parent_id,child2.id)
    sim.nodes[child1.id].driver=sim.nodes[parent_id].driver
    sim.nodes[child2.id].driver=sim.nodes[parent_id].driver
    sim.cellCompartments[comp_idx].nodesID[cell_idx]=child1.id
    push!(sim.cellCompartments[comp_idx].nodesID,child2.id)
    sim.totalpop+=1
    sim.cellCompartments[comp_idx].totalpop+=1
end
function NewNode(sim::CellSimulation, comp_idx::Int)
    sim.idx+=1
    cell=Node(sim.currentTime,1,sim.idx,nothing,[],false)
    sim.nodes[sim.idx]=cell
    return sim.nodes[sim.idx]
end
function getTotalRate(sim::CellSimulation)
    tot= sim.TotalDivisionRate + sim.TotalDeathRate 
    if sim.atEquilibrium 
        return tot 
    else 
        return sim.TotalDivisionRate
    end
end
function getTotalDivisionRate(sim::CellSimulation,alpha::Float64)
    tot=0.0
    for i in 1:sim.ncompart
        tot+=alpha*(1+sim.compartfitness[i])*sim.cellCompartments[i].totalpop
    end
    return tot 
end
function setRates(sim::CellSimulation,alpha::Float64,maxpop::Int)
    sim.TotalDivisionRate=getTotalDivisionRate(sim,alpha)
    y=(sim.TotalDivisionRate-alpha*(maxpop-sim.totalpop))
    sim.TotalDeathRate=maximum([y,0.0])
    return 
end
"""
    doEvent(sim::CellSimulation)


"""
function doEvent(sim::CellSimulation)
    if sim.atEquilibrium && rand()< sim.TotalDeathRate/(sim.TotalDeathRate+sim.TotalDivisionRate) #death
        prob=[comp.totalpop for comp in sim.cellCompartments]
        comp_idx=sample(1:length(sim.cellCompartments),Weights(prob))
        cell_idx=rand(1:sim.cellCompartments[comp_idx].totalpop)
        die(sim,comp_idx,cell_idx)
    else
        prob=[(1+sim.compartfitness[i])*sim.cellCompartments[i].totalpop for i in 1:sim.ncompart]
        comp_idx=sample(1:length(sim.cellCompartments),Weights(prob))
        cell_idx=rand(1:sim.cellCompartments[comp_idx].totalpop)
        divide(sim,comp_idx,cell_idx)
    end
end    
function snap(sim::CellSimulation)
    push!(sim.populationTrace,(sim.currentTime,sim.cellCompartments[1].totalpop,sim.cellCompartments[2].totalpop))
end
function empty_mutcompartment(sim::CellSimulation)
    sim.cellCompartments[2].totalpop=0
    sim.cellCompartments[2].nodesID=[]
end

"""
    get_default_config()

Gets the initial default configuration for the compartments and the simulation
"""
function get_default_config(fitness::Float64, driver_rate::Float64,ndrivers::Int, rate::Float64)
    C1=CellCompartment([],0,0.0,1,1)
    C2=CellCompartment([],0,0.0,0,0)
    root=Node(0.0,0,1,nothing,[],false)
    push!(C1.nodesID,root.id)
    sim=CellSimulation(1,[C1,C2],2,[0.0,fitness],0.0,driver_rate,Dict(root.id => root),0.0,0.0,[],ndrivers,0.0,0,0.0,rate,1,0,0.0)
    return sim 
end

# function reset_simulation(S::CellSimulation, fitness::Float64, driver_rate::Float64, ndrivers::Int, rate::Float64)
#     # Recreate the default configuration
#     C1 = CellCompartment([], 0, 0.0, 1, 1)
#     C2 = CellCompartment([], 0, 0.0, 0, 0)
#     root = Node(0.0, 0, 1, nothing, [], false)
#     push!(C1.nodesID, root.id)

#     # Reset the simulation to the new default state
#     S.totalpop=1
#     S.cellCompartments = [C1, C2]
#     S.ncompart=2
#     S.compartfitness= [0.0, fitness]
#     S.startTime=0.0
#     S.driverAcquisitionRate = driver_rate
#     S.nodes = Dict(root.id => root)
#     S.currentTime=0.0
#     S.lastSnap=0.0
#     S.populationTrace=[]
#     S.T_mut=0.0
#     S.atEquilibrium=0
#     S.TotalDeathRate=0.0
#     S.TotalDivisionRate=rate
#     S.ndrivers = ndrivers
#     S.status = 0 # Assuming 0 means the simulation is in the default initial state
#     S.idx = 1
#     S.CF=0.0

#     # If there are additional fields or objects that need to be reset, you can handle them here.
#     return S
# end
function reset_simulation(S::CellSimulation, fitness::Float64, driver_rate::Float64, ndrivers::Int, rate::Float64)
    # Instead of creating new compartments, reuse existing ones and clear their contents
    if length(S.cellCompartments) >= 2
        # Clear first compartment
        empty!(S.cellCompartments[1].nodesID)  # Empty nodesID array without reallocation
        S.cellCompartments[1].totalpop = 1   # Reset population to 1
        S.cellCompartments[1].DivisionRate = 1.0
        S.cellCompartments[1].ndiv = 1
        S.cellCompartments[1].active = true
        
        # Clear second compartment
        empty!(S.cellCompartments[2].nodesID)  # Empty nodesID array without reallocation
        S.cellCompartments[2].totalpop = 0
        S.cellCompartments[2].DivisionRate = 0.0
        S.cellCompartments[2].ndiv = 0
        S.cellCompartments[2].active = false
    else
        # If compartments don't exist yet, create them (should be rare)
        S.cellCompartments = [
            CellCompartment(Vector{Int}(), 1, 1.0, 1, true),
            CellCompartment(Vector{Int}(), 0, 0.0, 0, false)
        ]
    end
    
    # Clear and reset nodes dictionary instead of creating a new one
    empty!(S.nodes)
    
    # Create only one new node (the root) instead of recreating all nodes
    root = Node(0.0, 0, 1, nothing, Vector{Node}(), false)
    S.nodes[root.id] = root
    
    # Reset compartment node IDs
    push!(S.cellCompartments[1].nodesID, root.id)
    
    # Reset all simulation parameters (no new allocations)
    S.totalpop = 1
    S.ncompart = 2
    
    # Ensure compartfitness is the right size and reset values
    if length(S.compartfitness) < 2
        S.compartfitness = zeros(Float64, 2)
    end
    S.compartfitness[1] = 0.0
    S.compartfitness[2] = fitness
    
    S.startTime = 0.0
    S.driverAcquisitionRate = driver_rate
    S.currentTime = 0.0
    S.lastSnap = 0.0
    S.T_mut = 0.0
    S.atEquilibrium = false  # Changed from 0 to false to match Bool type
    S.TotalDeathRate = 0.0
    S.TotalDivisionRate = rate
    S.ndrivers = ndrivers
    S.status = 0
    S.idx = 1
    S.CF = 0.0
    
    # Clear the population trace without reallocating
    if isdefined(S, :populationTrace) && !isnothing(S.populationTrace)
        empty!(S.populationTrace)
    else
        S.populationTrace = Vector{Tuple}()  # Use the correct type
    end
    
    return S
end

"""
    run_selection_sim()

Models the evolution of a population of cells with a selective advantage given to a randomly chosen cell, once the population size reaches the specified equilibrium size. 
The selective advantage is introduced as a driver event. The driver is subject to stochastic extinction and so is repeatedly dropped in until it "takes".The prevailing state when the driver is introduced is saved and reinstated with each attempted introduction.

"""
function run_selection_sim(target_pop_size,nyears_driver_acquisistion, nyears, fitness; initial_division_rate=0.1, final_division_rate=1/365, minprop=0.0002,maxtry=100, max_driver_count=1,S=nothing)
    driver_rate=0.0
    #println("simulation started")
    if isnothing(S)
        S = get_default_config(fitness, driver_rate, 0, initial_division_rate)
    else
        reset_simulation(S,fitness, driver_rate, 0, initial_division_rate) # Reset the state for reuse
    end
    simulation(S, fitness,nyears_driver_acquisistion*365,target_pop_size, initial_division_rate,driver_rate, true,minprop,true, max_driver_count,false)
    tries=0
    #println("equilibrium reached")
    if S.status==0 #equilibrium population is reached 
        S.driverAcquisitionRate=0.0
        #phase1= S
        #dc=phase1.status #check if driver mutation population has escaped stochastic extinction 
        while S.status !=1 #cycle ends when the fraction of mutated cells before or at nyears_driver_acquisistion is big enough
            if tries >= maxtry #trying to acquire driver mutation 
                return nothing 
            end
            simulation(S, fitness, nyears_driver_acquisistion*365,target_pop_size, final_division_rate,0.0, false,minprop,true, max_driver_count,false) #simulation until nyears_driver_acquisistion or until the simulation breaks 
            empty_mutcompartment(S) #clean the second compartment
            if S.cellCompartments[1].totalpop <=0
                return nothing 
            end
            row_idx = rand(1:S.cellCompartments[1].totalpop)
            r=addDriver(S, row_idx, fitness) #adding the driver mutation in the empty compartment 
            #dc= temp_phase.status
            if r == nothing
                return nothing
            end
            simulation(S, fitness, nyears*365,target_pop_size, final_division_rate,0.0, false,minprop,true, max_driver_count,false)
            tries+=1
        end
        #println("mut acquired $(S.currentTime)")
        #temp_phase=phase1
        if abs(S.currentTime - nyears*365) > 1
            simulation(S, fitness, nyears*365,target_pop_size, final_division_rate,0.0, false,10.0,false, max_driver_count,false)
        end
        if S.cellCompartments[2].totalpop <=20
            return nothing 
        end
         #simulation until nyears when the driver mutation is acquired
    else 
        println("Driver mutation acquired before equilibrium is reached!")
        return 
    end
    return S
end

# function bulk_update_lineages(sim::CellSimulation, comp_idx::Int, div_events::Int, death_events::Int)
#     comp = sim.cellCompartments[comp_idx]
#     node_ids = comp.nodesID
#     current_time = sim.currentTime

#     # Process divisions
#     new_nodes = Tuple{Int,Int,Int}[]  # (parent_id, child1_id, child2_id)
#     sizehint!(new_nodes, div_events)
    
#     for _ in 1:div_events
#         isempty(node_ids) && break
        
#         # Select random parent
#         parent_idx = rand(1:length(node_ids))
#         parent_id = node_ids[parent_idx]
#         parent_node = sim.nodes[parent_id]
        
#         # Create new children
#         child1_id = sim.idx += 1
#         child2_id = sim.idx += 1
        
#         child1 = Node(current_time, 1, child1_id, parent_node, Node[], parent_node.driver)
#         child2 = Node(current_time, 1, child2_id, parent_node, Node[], parent_node.driver)
        
#         # Update parent's children
#         empty!(parent_node.children)
#         push!(parent_node.children, child1, child2)
        
#         # Store new nodes and update compartment
#         sim.nodes[child1_id] = child1
#         sim.nodes[child2_id] = child2
#         node_ids[parent_idx] = child1_id
#         push!(node_ids, child2_id)
        
#         push!(new_nodes, (parent_id, child1_id, child2_id))
#     end

#     # Process deaths in reverse order
#     death_indices = unique(sort(rand(1:length(node_ids), death_events), rev=true))
#     for idx in death_indices
#         idx > length(node_ids) && continue
        
#         node_id = node_ids[idx]
#         node = sim.nodes[node_id]
#         parent = node.parent
        
#         if parent !== nothing && length(parent.children) == 2
#             # Find surviving sibling
#             surviving_sibling = parent.children[parent.children[1].id == node_id ? 2 : 1]
            
#             # Update grandparent relationship
#             grandparent = parent.parent
#             if grandparent !== nothing
#                 filter!(x -> x.id != parent.id, grandparent.children)
#                 push!(grandparent.children, surviving_sibling)
#                 surviving_sibling.parent = grandparent
#             end
            
#             # Update sibling properties
#             surviving_sibling.tBirth = parent.tBirth
#             surviving_sibling.ndiv += parent.ndiv
#             surviving_sibling.driver = parent.driver
            
#             # Cleanup parent
#             delete!(sim.nodes, parent.id)
#         end
        
#         # Remove node
#         deleteat!(node_ids, idx)
#         delete!(sim.nodes, node_id)
#     end

#     # Update population counts
#     net_change = div_events - death_events
#     comp.totalpop += net_change
#     sim.totalpop += net_change
# end

# function gillespie_step(sim::CellSimulation, comp_idx::Int)
#     comp = sim.cellCompartments[comp_idx]
#     totdeathrate=comp.DivisionRate * comp.totalpop-0.1*(10^5-comp.totalpop)
#     total_rate = comp.DivisionRate * comp.totalpop + maximum([totdeathrate,0.0])
    
#     if total_rate ≈ 0.0
#         sim.currentTime += Inf
#         return
#     end
    
#     # Generate next event time
#     Δt = rand(Exponential(1/total_rate))
#     sim.currentTime += Δt
    
#     # Select event type
#     if rand() < (comp.DivisionRate * comp.totalpop) / total_rate
#         # Division event
#         if !isempty(comp.nodesID)
#             cell_idx = rand(1:length(comp.nodesID))
#             divide(sim, comp_idx, cell_idx)
#         end
#     else
#         # Death event
#         if !isempty(comp.nodesID)
#             cell_idx = rand(1:length(comp.nodesID))
#             die(sim, comp_idx, cell_idx)
#         end
#     end
# end

# function exact_lineage_updates(sim::CellSimulation, comp_idx::Int, div_events::Int, death_events::Int)
#     comp = sim.cellCompartments[comp_idx]
#     node_ids = comp.nodesID
    
#     # Process divisions using existing divide function
#     for _ in 1:div_events
#         isempty(node_ids) && break
#         cell_idx = rand(1:length(node_ids))
#         divide(sim, comp_idx, cell_idx)
#     end
    
#     # Process deaths using existing die function in reverse order
#     death_indices = unique(sort(rand(1:length(node_ids), death_events), rev=true))
#     for idx in death_indices
#         idx > length(node_ids) && continue
#         die(sim, comp_idx, idx)
#     end
# end
function simulation(sim::CellSimulation, fitness::Float64 ,stopTime::Float64,maxpop::Int,alpha::Float64,driverRate::Float64,break_at_equilibrium::Bool, minprop::Float64, break_if_driver::Bool,maxdrivers::Int, tau_approx::Bool)
    threshold=(stopTime - sim.startTime)/100
    snap(sim)
    while sim.currentTime < stopTime
        totrate=getTotalRate(sim)
        if tau_approx==true 
            τ = 1 # Adaptive τ calculation would be better
                div_events = rand(Poisson(sim.TotalDivisionRate * τ))
                death_events = rand(Poisson(sim.TotalDeathRate * τ))
                sim.currentTime+=τ
                # Bulk update population
                net_change = div_events - death_events
                #println(net_change)
                if net_change < 0
                    for _ in 1:net_change
                        prob=[comp.totalpop for comp in sim.cellCompartments]
                        comp_idx=sample(1:length(sim.cellCompartments),Weights(prob))
                        cell_idx=rand(1:sim.cellCompartments[comp_idx].totalpop)
                        die(sim,comp_idx,cell_idx)
                    end
                elseif net_change >0 
                    for _ in 1:net_change 
                        prob=[(1+sim.compartfitness[i])*sim.cellCompartments[i].totalpop for i in 1:sim.ncompart]
                        comp_idx=sample(1:length(sim.cellCompartments),Weights(prob))
                        cell_idx=rand(1:sim.cellCompartments[comp_idx].totalpop)
                        divide(sim,comp_idx,cell_idx)
                    end
                end
        else 
            sim.currentTime+=rand(Exponential(1/totrate))
            doEvent(sim)
        end
        if !sim.atEquilibrium && sim.totalpop >=maxpop
            sim.atEquilibrium=1
            if break_at_equilibrium
                break
            end
        end
        if sim.atEquilibrium && sim.totalpop<maxpop
            sim.atEquilibrium=0
        end
        sim.TotalDivisionRate = max(sim.TotalDivisionRate, 1e-8)
        sim.TotalDeathRate = max(sim.TotalDeathRate, 1e-8)
        setRates(sim,alpha,maxpop)
        if sim.currentTime-sim.lastSnap >threshold
            snap(sim)
            sim.lastSnap=sim.currentTime
        end
        m=(sim.cellCompartments[2].totalpop/(sim.cellCompartments[1].totalpop+sim.cellCompartments[2].totalpop))
        sim.CF=m
        if m>=minprop 
            sim.status=1
            if break_if_driver
                break
            end
        end
    end
    return sim
end


# function simulation(sim::CellSimulation, fitness::Float64 ,stopTime::Float64,maxpop::Int,alpha::Float64,driverRate::Float64,break_at_equilibrium::Bool, minprop::Float64, break_if_driver::Bool,maxdrivers::Int, tau_approx::Bool)
#     threshold=(stopTime - sim.startTime)/100
#     snap(sim)
#     while sim.currentTime < stopTime
#         totrate=getTotalRate(sim)
#         # Process each compartment separately
#         for comp_idx in eachindex(sim.cellCompartments)
#             comp = sim.cellCompartments[comp_idx]
            
#             if comp.totalpop > 10^4  # Tau-leaping threshold
#                 # Calculate tau-leaping parameters
#                 τ = 0.1*365  # Adaptive τ calculation would be better
#                 div_events = rand(Poisson(comp.division_rate * comp.totalpop * τ))
#                 death_events = rand(Poisson(comp.death_rate * comp.totalpop * τ))
                
#                 # Bulk update population
#                 net_change = div_events - death_events
#                 comp.totalpop += net_change
#                 sim.totalpop += net_change
                
#                 # Approximate lineage updates for large populations
#                 if comp_idx == 1  # Wild-type compartment
#                     bulk_update_lineages(sim, comp_idx, div_events, death_events)
#                 else
#                     exact_lineage_updates(sim, comp_idx, div_events, death_events)
#                 end
                
#                 sim.currentTime += τ
#             else
#                 # Use Gillespie for small populations
#                 gillespie_step(sim, comp_idx)
#             end
#         end
        
#         if !sim.atEquilibrium && sim.totalpop >=maxpop
#             sim.atEquilibrium=1
#             if break_at_equilibrium
#                 break
#             end
#         end
#         if sim.atEquilibrium && sim.totalpop<maxpop
#             sim.atEquilibrium=0
#         end
#         sim.TotalDivisionRate = max(sim.TotalDivisionRate, 1e-8)
#         sim.TotalDeathRate = max(sim.TotalDeathRate, 1e-8)
#         setRates(sim,alpha,maxpop)
#         if sim.currentTime-sim.lastSnap >threshold
#             snap(sim)
#             sim.lastSnap=sim.currentTime
#         end
#         m=(sim.cellCompartments[2].totalpop/sim.cellCompartments[1].totalpop)
#         if m>=minprop 
#             sim.status=1
#             if break_if_driver
#                 break
#             end
#         end
#     end
#     return sim
# end

# end
#     while sim.currentTime < stopTime
#         totrate=getTotalRate(sim)
#         #nextStep= max(rand(Exponential(1/(totrate+sim.totalpop*sim.driverAcquisitionRate))),1.0)
#         nextStep=rand(Exponential(1/(totrate+sim.totalpop)))
#         sim.currentTime+=nextStep
#         # if rand()<sim.totalpop*sim.driverAcquisitionRate/(totrate+sim.totalpop*sim.driverAcquisitionRate) && sim.ndrivers <maxdrivers #driver mutation appear 
#         #     row_idx = rand(1:sim.cellCompartments[1].totalpop)
#         #     addDriver(sim, row_idx,fitness)
#         # end 
#         doEvent(sim)
#         if !sim.atEquilibrium && sim.totalpop >=maxpop
#             sim.atEquilibrium=1
#             if break_at_equilibrium
#                 break
#             end
#         end
#         if sim.atEquilibrium && sim.totalpop<maxpop
#             sim.atEquilibrium=0
#         end
#         sim.TotalDivisionRate = max(sim.TotalDivisionRate, 1e-8)
#         sim.TotalDeathRate = max(sim.TotalDeathRate, 1e-8)
#         setRates(sim,alpha,maxpop)
#         if sim.currentTime-sim.lastSnap >threshold
#             snap(sim)
#             sim.lastSnap=sim.currentTime
#         end
#         m=(sim.cellCompartments[2].totalpop/sim.cellCompartments[1].totalpop)
#         if m>=minprop 
#             sim.status=1
#             if break_if_driver
#                 break
#             end
#         end
#     end
#     return sim
# end


# function construct_final_tree(sim::CellSimulation,n_leaves::Int,mutrateperdivision::Float64, odf::Float64, backgroundrate::Float64)
#     n=collect(values(sim.nodes))
#     root=n[findfirst(x -> x.parent === nothing,n)]
#     nw_string=to_newick_with_node_height_dev(root,odf, backgroundrate)
#     tree=parsenewick(nw_string)
#     t=subsample_mutated_tree(sim,tree,n_leaves)
#     return t 
# end


function construct_mutated_tree(sim::CellSimulation,n_leaves::Int, backgroundrate::Float64;  odf=2.0)
    n=collect(values(sim.nodes))
    root=n[findfirst(x -> x.parent === nothing,n)]
    nw_string=to_newick_with_node_height_dev(root,odf,backgroundrate)
    tree=parsenewick(nw_string)
    t=subsample_mutated_tree(sim,tree,n_leaves)
    return t
end
# function construct_time_tree(sim::CellSimulation,n_leaves::Int, backgroundrate::Float64)
#     n=collect(values(sim.nodes))
#     root=n[findfirst(x -> x.parent === nothing,n)]
#     nw_string=to_newick_with_node_height_time(root,backgroundrate)
#     tree=parsenewick(nw_string)
#     t=subsample_mutated_tree(sim,tree,n_leaves)
# end
# function construct_mutation_tree(sim::CellSimulation,n_leaves::Int,odf::Float64, backgroundrate::Float64)
#     n=collect(values(sim.nodes))
#     root=n[findfirst(x -> x.parent === nothing,n)]
#     nw_string=to_newick_with_node_height(root,odf,backgroundrate)
#     tree=parsenewick(nw_string)
#     t=subsample_mutated_tree(sim,tree,n_leaves)
# end

function get_nb(n, meanmuts, od)
    return rand(NegativeBinomial.(meanmuts ./ (od - 1), 1 / od), n)
end
function to_newick_with_node_height_dev(root::Node,odf::Float64, backgroundrate::Float64)
    # Recursive helper function to construct the Newick string
    function recurse_tree(node::Node, parent_height::Int)
        children = node.children
        current_height = node_height(node,odf,backgroundrate)
        #if parent_height < current_height 
            branch_length =current_height #- parent_height
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

    # Function to calculate the height of a node
    function node_height(node::Node,odf::Float64, backgroundrate::Float64)
        if node.parent === nothing
            height=0
        else
            if odf > 1
                if node.parent.tBirth == 0.0 
                    height= get_nb(1, totalLambda(backgroundrate,node.tBirth/365,0), odf)[1]
                else
                    height = get_nb(1, abs(totalLambda(backgroundrate,node.parent.tBirth/365,0) - totalLambda(backgroundrate,node.tBirth/365,0) + 1e-5) , odf)[1]
                end
            else
                if node.parent.tBirth == 0.0 
                    height= rand(Poisson(totalLambda(backgroundrate,node.tBirth/365,0)))
                else
                    height = rand(Poisson(abs(totalLambda(backgroundrate,node.parent.tBirth/365,0) -totalLambda(backgroundrate,node.tBirth/365,0))))
                end
            end
        end
        return height
    end

    # Start recursion from the root node; root has no branch length
    #return recurse_tree(root, node_height(root,odf,backgroundrate)) * ";"
    return recurse_tree(root, 0) * ";"
end


# function to_newick_with_node_height_time(root::Node, backgroundrate::Float64)
#     # Recursive helper function to construct the Newick string
#     function recurse_tree(node::Node, parent_height::Float64)
#         children = node.children
#         current_height = node_height(node,backgroundrate)
#         #if parent_height < current_height 
#             branch_length = current_height #-parent_height
#         #else
#             #branch_length=0
#         #end
        
#         if !isempty(children)
#             # Internal node: recurse on all children
#             subtree = join([recurse_tree(child, current_height) for child in children], ",")
#             return "($subtree)$(node.id):$branch_length"
#         else
#             # Leaf node: return the node's ID with the correct branch length
#             return "$(node.id):$branch_length"
#         end
#     end

#     # Function to calculate the height of a node
#     function node_height(node::Node, backgroundrate::Float64)
#         if node.parent !== nothing
#             height=((node.tBirth -node.parent.tBirth)/365)*backgroundrate
#         else 
#             height=0.0
#         end
#         return height
#     end

#     # Start recursion from the root node; root has no branch length
#     return recurse_tree(root, 0.0) * ";"
# end

# function to_newick_with_node_height(root::Node,odf::Float64, backgroundrate::Float64)
#     # Recursive helper function to construct the Newick string
#     function recurse_tree(node::Node, parent_height::Int)
#         children = node.children
#         current_height = node_height(node,odf,backgroundrate)
#         #if parent_height < current_height 
#             branch_length =current_height #- parent_height
#         #else
#             #branch_length=0
#         #end
        
#         if !isempty(children)
#             # Internal node: recurse on all children
#             subtree = join([recurse_tree(child, current_height) for child in children], ",")
#             return "($subtree)$(node.id):$branch_length"
#         else
#             # Leaf node: return the node's ID with the correct branch length
#             return "$(node.id):$branch_length"
#         end
#     end

#     # Function to calculate the height of a node
#     function node_height(node::Node,odf::Float64, backgroundrate::Float64)
#         if node.tBirth==0.0
#             height=0
#         else
#             if odf > 1
#                 if node.parent.tBirth ==0.0 
#                     height= get_nb(1, backgroundrate*node.tBirth/365, odf)[1]
#                 else
#                     m=abs(totalLambda(backgroundrate,node.parent.tBirth/365,0) - totalLambda(backgroundrate,node.tBirth/365,0))
#                     if m !=NaN 
#                         height = get_nb(1, max , odf)[1]
#                     else 
#                         println("Error in node $(node.id) with birth $(node.tBirth)")
#                     end
#                 end
#             else
#                 if node.parent.tBirth ==0.0 
#                     height= rand(Poisson(backgroundrate*node.tBirth/365))
#                 else
#                     height = rand(Poisson(abs(backgroundrate*node.parent.tBirth/365 -backgroundrate*node.tBirth/365)))
#                 end
#             end
#         end
#         return height
#     end

#     # Start recursion from the root node; root has no branch length
#     return recurse_tree(root, node_height(root,odf,backgroundrate)) * ";"
# end
function subsample_mutated_tree(sim,tree,n)
    leaves=getleaves(tree)
    mut_leaves=[el for el in leaves if sim.nodes[parse(Int,el.name)].driver==true]
    #Random.seed!(11)
    s=rand(mut_leaves, n)
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
