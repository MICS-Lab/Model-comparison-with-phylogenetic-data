using Plots, Distributions,Statistics 
using MCPhyloTree
using DataFrames, DelimitedFiles
using CSV 
using Distances
using Random

"""
function LTT_plot(t,L,g,k) to create LTT plot for a tree 
inputs: t=root of the tree
outputs: height= vector containing the height of each node
        mutated= vector containing the number of lineages at each node 
"""
function LTT_plot(t,l,L,add)
    nh = nodeheights(t,noleaves=true) 
    lh=nodeheights(t,onlyleaves=true) 
    h=mean(lh)
    index= sortperm(nh)
    nh=nh[index]
    c=nh.axes[1]
    mut=[]
    m=1
    height=[]
    i=1
    j=0
    for el in nh
        if el>h
            el=h
        end
        push!(mut,m)
        push!(height,el+add)
        m+=length(getchildren(t,c[i]))-1
        push!(height,el+add)
        push!(mut,m)
        i+=1
    end
    #push!(mut,m)
    #push!(height,(23+T)*h/l)
    pushfirst!(mut,1)
    pushfirst!(height,0.0)
    if length(height) !=length(mut)
        println("ERROR with constructing the LTT plot, nh= $(length(nh)) h= $(length(height)) and m=$(length(mut))")
    end
    return [height,mut]
end

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

function ABC(file_name,data,points,Tmin,Tmax,smin,smax,N,L,k,l, initial_division_rate, final_division_rate,minprop,maxtry,max_driver_count,odf)
    open(file_name, "w") do file
    println(file, "S,T_m,d")
    backgroundrate=l
    i=1
    p=0
    traj_h=[]
    traj_m=[]
    fit=[]
    T=[]
    while i <= points 
        let 
        T_min=sample_acq_time(Tmin, l, L, 1, 0)[1]
        T_max=sample_acq_time(Tmax - Tmin, l, L, 1, T_min)[1]
        s_y=rand(Uniform(smin,smax))
        T_m=rand(Uniform(T_min,T_max))
        s=log(1+s_y)/(365*final_division_rate)
        S=run_selection_sim(initial_division_rate, final_division_rate, N,T_m, L, s, minprop,maxtry, max_driver_count)
        if S !==nothing 
            t=construct_mutated_tree(S,k, odf, backgroundrate)
            a=LTT_plot(t,l, L)
            diff=piecewise_difference_g(a[1],a[2],data[1],data[2])
            d= area(diff[1],diff[2])/k
            push!(traj_h,a[1])
            push!(traj_m,a[2])
            push!(T,T_m)
            push!(fit,s_y)
            i+= 1
            println(i)
            println(file,"$s_y,$T_m,$(abs(d))")
            flush(file)
        end
    end
    end
    return [traj_h,traj_m,T,fit]
end
end
function ABC_prob(file_name,data,points,Tmin,Tmax,smin,smax,N,L,k,l, initial_division_rate, final_division_rate,minprop,maxtry,max_driver_count,odf)
    open(file_name, "w") do file
    println(file, "S,T_m,d")
    name="prob" * file_name
    open(name, "w") do file1
    println(file1, "Extinction probability")
    backgroundrate=l
    i=1
    p=0
    tot=0
    while i <= points 
        let 
        #T_min=sample_acq_time(Tmin, l, L, 1, 0)[1]
        #T_max=sample_acq_time(Tmax - Tmin, l, L, 1, T_min)[1]
        T_min=Tmin/l  
        T_max=Tmax/l
        s_y=rand(Uniform(smin,smax))
        T_m=rand(Uniform(T_min,T_max))
        s=log(1+s_y)/(365*final_division_rate)
        S=run_selection_sim(initial_division_rate, final_division_rate, N,T_m, L, s, minprop,maxtry, max_driver_count)
        if S !==nothing 
            t=construct_mutated_tree(S,k, odf, backgroundrate)
            a=LTT_plot(t,l, L,0.0)
            diff=piecewise_difference_g(a[1],a[2],data[1],data[2])
            d= area(diff[1],diff[2])/k
            i+= 1
            println(i)
            println(file,"$s_y,$T_m,$(abs(d))")
            flush(file)
        else 
            p+=1
            println(file1, "$s_y, $T_m")
            flush(file)
        end
        tot+=1
    end
    end
    println(file1,"probability")
    println(file1, "$(p/tot)")
    end
    return
end
end
function ABC_constantl(file_name,data,points,Tmin,Tmax,smin,smax,N,L,k,l, initial_division_rate, final_division_rate,minprop,maxtry,max_driver_count)
    open(file_name, "w") do file
    println(file, "S,T_m,d")
    backgroundrate=l
    i=1
    p=0
    while i <= points 
        let 
        T_min=sample_acq_time(Tmin, l, L, 1, 0)[1]
        T_max=sample_acq_time(Tmax - Tmin, l, L, 1, T_min)[1]
        s_y=rand(Uniform(smin,smax))
        T_m=rand(Uniform(T_min,T_max))
        s=log(1+s_y)/(365*final_division_rate)
        S=run_selection_sim(initial_division_rate, final_division_rate, N,T_m, L, s, minprop,maxtry, max_driver_count)
        if S !==nothing 
            t=construct_time_tree(S,k, backgroundrate)
            global a=LTT_plot(t,l, L)
            diff=piecewise_difference_g(a[1],a[2],data[1],data[2])
            d= area(diff[1],diff[2])/k
            i+= 1
            println(i)
            println(file,"$s_y,$T_m,$(abs(d))")
            flush(file)
        end
    end
    end
    return 
end
end
function ABC_constant_rate(file_name,data,points,Tmin,Tmax,smin,smax,N,L,k,l, initial_division_rate, final_division_rate,minprop,maxtry,max_driver_count,odf)
    open(file_name, "w") do file
    println(file, "S,T_m,d")
    backgroundrate=l
    i=1
    p=0
    while i <= points 
        let 
        T_min=sample_acq_time(Tmin, l, L, 1, 0)[1]
        T_max=sample_acq_time(Tmax - Tmin, l, L, 1, T_min)[1]
        s_y=rand(Uniform(smin,smax))
        T_m=rand(Uniform(T_min,T_max))
        s=log(1+s_y)/(365*final_division_rate)
        S=run_selection_sim(initial_division_rate, final_division_rate, N,T_m, L, s, minprop,maxtry, max_driver_count)
        if S !==nothing 
            t=construct_mutation_tree(S,k,odf, backgroundrate)
            global a=LTT_plot(t,l, L)
            diff=piecewise_difference_g(a[1],a[2],data[1],data[2])
            d= area(diff[1],diff[2])/k
            i+= 1
            println(i)
            println(file,"$s_y,$T_m,$(abs(d))")
            flush(file)
        end
    end
    end
    return 
end
end