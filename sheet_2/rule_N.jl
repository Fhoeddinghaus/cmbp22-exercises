# Implements rule N as callable functions

# for plotting
using Colors, Plots, LaTeXStrings
plotlyjs()

## convert N into the ruleset/table/map
N₂_map(N) = digits(N, base=2, pad=8) # in order of 1,2,4,8,16,32,64,128

## calculate the new cell value
# get value of b' = f(a,b,c)
function f(N::Int, abc::String)
    abc = parse(Int, abc, base=2)
    return N₂_map(N)[abc + 1] # note: julia counts from 1 not from 0!
end


# get cell value at i from current configuration
# configuration z is vector with indices from 1 to M
function f(N::Int, z::Vector, i::Int)
    if (i < 1) || (i > length(z)) 
        throw(BoundsError)
        return nothing
    end
    
    # get current values of neighbours
    # previous (i-1)
    a = z[end] # periodic boundary condition
    if i > 1
        # get a from configuration
        a = z[i-1]
    end
    # current (i)
    b = z[i]
    # next (i+1)
    c = z[begin] # periodic boundary condition
    if i < length(z)
        c = z[i+1]
    end
    # concat to make abc
    abc = string(a,b,c)
    return f(N, abc)
end

# calculate next configuration from current
function next_z(N::Int, z::Vector)
    z′ = zeros(Int, length(z))
    for i in 1:length(z)
        z′[i] = f(N, z, i)
    end
    return z′
end


# A few global configurations (overwrite in main file)
NUMBER_OF_SITES = 120
N_rule = 184
NUMBER_OF_TIMESTEPS = 50

function calculate_rule_N(
        z_start, # start configuration
        N_rule # wolfram code
    )
    # store all configurations
    zs = zeros(Int, NUMBER_OF_TIMESTEPS+1, NUMBER_OF_SITES)
    zs[1,:] = z_start

    # calculate the next generations
    for t in 2:(NUMBER_OF_TIMESTEPS+1)
        zs[t,:] = next_z(N_rule, zs[t-1,:])
    end
    return zs
end

# plotting
function plot_rule_N(zs)
    heatmap(
        1:1:NUMBER_OF_SITES, # x-axis label
        0:1:NUMBER_OF_TIMESTEPS,  # y-axis label
        zs, # all configurations to plot
        color = reverse(cgrad(:greys)), # reverse the map (0=white, 1=black)
        yflip=true, # t=0 at top
        size=(NUMBER_OF_SITES*10/4*3,NUMBER_OF_TIMESTEPS * 10/4*3), 
        legend=false,
        xlabel="i",
        ylabel="Time t"
    )
end

number_active_cells(zs) = [sum(zs[t+1,:]) for t in 0:50]
function plot_number_active_cells(ns)
    plot(0:NUMBER_OF_TIMESTEPS, ns, xlabel="t", ylabel="n(t)", label="n_alive(t)")
end

;


