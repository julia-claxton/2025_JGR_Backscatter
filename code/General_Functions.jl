using Statistics, LinearAlgebra
# This is a file of general-purpose functions that I find useful, compiled in one
# script so I can use them in all my scripts.

# Written by Julia Claxton. Contact: julia.claxton@colorado.edu
# Released under MIT License.

function print_progress_bar(fraction; bar_length = 30)
# Prints a progress bar to terminal filled to user-specified percentage.
    character_length = 1/bar_length
    number_of_filled_characters = Int(floor(fraction/character_length))
    print("\r")
    print(repeat("█", number_of_filled_characters))
    print(repeat("░", bar_length - number_of_filled_characters))
    print(" [$(round(fraction*100, digits = 1))%]")
    if fraction == 1; println(); end
end

function edges_to_means(edges)
# Calculates means of bins given their edges.
    return [mean([edges[i], edges[i+1]]) for i = 1:length(edges)-1] 
end

function exact_2dhistogram(x, y; weights = ones(length(x)), max_bins = Inf)
# Bin 2D data into a histogram with automatically-selected bin edges.
# This can be very slow, use only for exploratory analysis!
    n_datapoints = length(x)
    xmin = min(x...); xmax = max(x...)
    ymin = min(y...); ymax = max(y...)

    # Freedman–Diaconis rule for bin spacing
    x_bin_space = 2*std(x)/(length(x)^(1/3))
    y_bin_space = 2*std(y)/(length(y)^(1/3))

    # Construct bins
    x_bin_edges = xmin:x_bin_space:xmax
    y_bin_edges = ymin:y_bin_space:ymax

    # Clamp bins if max_bins was provided
    if length(x_bin_edges)-1 > max_bins
        x_bin_edges = LinRange(xmin, xmax, max_bins+1)
    end
    if length(y_bin_edges)-1 > max_bins
        y_bin_edges = LinRange(ymin, ymax, max_bins+1)
    end
    
    # Get frequencies
    frequencies = exact_2dhistogram(x, y, x_bin_edges, y_bin_edges, weights = weights)

    # Return
    return x_bin_edges, y_bin_edges, frequencies
end

function exact_2dhistogram(x, y, x_bin_edges, y_bin_edges; weights = ones(length(x)))
# Bin 2D data into a histogram with bin edges defined by user.
    @assert length(x) == length(y) == length(weights) "x, y, and weight vectors must be same length"

    x_nbins = length(x_bin_edges) - 1
    y_nbins = length(y_bin_edges) - 1
    
    result = zeros(x_nbins, y_nbins)
    for x_idx = 1:x_nbins
        x_slice = x_bin_edges[x_idx] .<= x .< x_bin_edges[x_idx+1]
        result[x_idx,:] += exact_1dhistogram(y[x_slice], y_bin_edges, weights = weights[x_slice])
    end
    return result
end

function exact_1dhistogram(data, bin_edges; weights = ones(length(data)))
# Bin 1D data into a histogram with bin edges defined by user.
    @assert length(data) == length(weights) "Data and weight vectors must be same length"
    nbins = length(bin_edges) - 1
    
    result = zeros(nbins)
    for idx = 1:nbins
        slice = bin_edges[idx] .<= data .< bin_edges[idx+1]
        result[idx] = sum(weights[slice])
    end
    return result
end
    
function meshgrid(x, y)
# Given x and y axis labels, create tuples of coordinates for all points on a grid formed by x and y
    grid = [(x[i], y[j]) for i in eachindex(x), j in eachindex(y)]
    return reshape(grid, :) # Vectorize
end

function box_aspect!(aspect_ratio)
# Set box aspect ratio for a plot using Plots library
    xlims = Plots.xlims(plot!())
    ylims = Plots.ylims(plot!())
    square_aspect = (xlims[2]-xlims[1])/(ylims[2]-ylims[1])

    plot!(aspect_ratio = square_aspect * aspect_ratio)
end