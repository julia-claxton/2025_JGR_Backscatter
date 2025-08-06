const TOP_LEVEL = dirname(@__DIR__)

using Statistics, LinearAlgebra
using BenchmarkTools, Profile, TickTock
using Plots, Plots.PlotMeasures
using DelimitedFiles
using Dates
using NPZ
using NumericalIntegration
using LaTeXStrings

include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Events.jl")
include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Visualization.jl")
include("$(TOP_LEVEL)/code/G4EPP_2.0/Frontend_Functions.jl")
include("$(TOP_LEVEL)/code/General_Functions.jl")

"""
======================================
Main Figures
======================================
"""
function figure_elfin_data()
    event = create_event(DateTime("2021-03-05T03:38:30.084"),DateTime("2021-03-05T03:42:45.159"), "B", maximum_relative_error = .5)
    plot_event(event)
    idxs_to_plot = [5, 52, 67]
    for i in idxs_to_plot
        vline!([event.time[i]],
            label = false,
            linewidth = 2,
            linecolor = RGB(0xff006a)
        )
    end
    p1 = plot!()

    bottom_plots = []
    for i in idxs_to_plot
        pad_heatmap(event, i)
        plot!(
            colorbar_title = "Log10 Electron Flux\n# s⁻¹ cm⁻² str⁻¹ MeV⁻¹",
            colorbar = false,
            clims = (2.5, 6)
        )
        push!(bottom_plots, plot!())
    end

    bottom_plots[begin] = plot(bottom_plots[begin],
        ylabel = "Energy, keV",
        leftmargin = 5mm
    )
    bottom_plots[end] = plot(bottom_plots[end],
        colorbar = true,
        rightmargin = 5mm
    )

    layout = @layout [a; [b{.27w} c{.29w} d{.43w}]]
    plot(p1, bottom_plots...,
        layout = layout,
        size = (1,.55) .* 600,
        dpi = 400,
        thickness_scaling = .7
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/elfin_data.png")
end

function figure_data_coverage()
    data = npzread("$(TOP_LEVEL)/data/figure_data/data_coverage.npz")
    MLT_bin_edges_rad = data["MLT_bin_edges_rad"]
    L_bin_edges = data["L_bin_edges"]
    coverage = data["coverage"]

    to_plot = coverage'
    to_plot = replace(to_plot, -Inf => -100) # Make regions of no coverage black rather than transparent
    heatmap(MLT_bin_edges_rad, L_bin_edges, to_plot,
        projection = :polar,
        axis = false,

        ylims = (0, 10),

        colorbar_title = "# Data Segments",
        clims = (0, 150),
        colormap = :ice,

        size = (1.3, 1) .* 500,
        leftmargin = -5mm,

        background = :white,
        dpi = 400
    )
    plot_ring!.(0:1:10)
    plot_mlt!.(0:3:23)
    plot_earth!()

    annotate!(1.05,  0.0, text("0000", 10, :left))
    annotate!(0.00,  1.1, text("0600", 10, :center))
    annotate!(-1.05, 0.0, text("1200", 10, :right))
    annotate!(0.00, -1.1, text("1800", 10, :center))
    
    ϵ = -0.015
    annotate!(.2 + ϵ, 0, text("2", 10, :left, :white))
    annotate!(.4 + ϵ, 0, text("4", 10, :left, :white))
    annotate!(.6 + ϵ, 0, text("6", 10, :left, :white))
    annotate!(.8 + ϵ, 0, text("8", 10, :left, :white))

    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/data_coverage.png")
end

function figure_filling_vs_backscatter()
    data = npzread("$(TOP_LEVEL)/data/figure_data/backscatter_vs_filling.npz")

    energy_means = data["energy_means"]

    ambient_pa_means = data["ambient_pa_means"]
    ambient_pad = data["ambient_pad"]
    ambient_rn = data["ambient_rn"]
    ambient_lc = data["ambient_lc"]

    filling_pa_means = data["filling_pa_means"]
    filling_pad = data["filling_pad"]
    filling_rn = data["filling_rn"]
    filling_lc = data["filling_lc"]

    between_pa_means = data["between_pa_means"]
    between_pad = data["between_pad"]
    between_rn = data["between_rn"]
    between_lc = data["between_lc"]


    generate_pad_plot(ambient_pa_means, ambient_pad, ambient_lc, 180-ambient_lc, ambient_rn, annotate = true)
    p1 = plot!(
        ylabel = "Energy, keV",
        leftmargin = 10mm,
        rightmargin = 5mm
    )
    generate_pad_plot(filling_pa_means, filling_pad, filling_lc, 180-filling_lc, filling_rn, annotate = true)
    p2 = plot!(
        ylabel = "Energy, keV",
        leftmargin = 5mm,
        rightmargin = 5mm
    )
    generate_pad_plot(between_pa_means, between_pad, between_lc, 180-between_lc, between_rn, annotate = true)
    p3 = plot!(
        ylabel = "Energy, keV",
        leftmargin = 5mm
    )

    plot(p1, p2, p3,
        layout = (1,3),
        size = (2.4, .6) .* 550,
        dpi = 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/filling_vs_backscatter.png")
end

function figure_beam_weighting()
    count_clims = (0, 3.5)
    nflux_clims = (2, 6)

    data = npzread("$(TOP_LEVEL)/data/figure_data/beam_weighting.npz")
    avg_pitch_angles = data["avg_pitch_angles"]
    avg_energies = data["avg_energies"]
    data_nfluence = data["data_nfluence"]
    data_counts = data["data_counts"]
    culled_backscatter_input_distribution = data["culled_backscatter_input_distribution"]
    backscatter_input_distribution = data["backscatter_input_distribution"]
    beam_energies = data["beam_energies"]
    beam_pitch_angles = data["beam_pitch_angles"]
    beam_weights = data["beam_weights"]
    backscatter_counts = data["backscatter_counts"]
    corrected_backscatter_counts = data["corrected_backscatter_counts"]
    cull_idxs = data["cull_idxs"]

    plot_distribution(avg_energies, avg_pitch_angles, data_nfluence, :magma)
    plot!(
        title = "1) In-Situ Measured Fluence",
        colorbar_title = "Log10 Electron Fluence, # cm⁻² str⁻¹ MeV⁻¹",
        clims = nflux_clims,
        leftmargin = 10mm
    )
    p1 = plot!()


    plot_distribution(avg_energies, avg_pitch_angles, data_counts, :ice)
    plot!(
        title = "2) Convert to Counts",
        colorbar_title = "Log10 Electron Count, #",
        clims = count_clims
    )
    p2 = plot!()

    plot_distribution(avg_energies, avg_pitch_angles, backscatter_input_distribution, :ice)
    beams_to_plot = beam_weights .≠ 0
    scatter!(beam_pitch_angles[beams_to_plot], beam_energies[beams_to_plot],
        title = "3) Azimuth Scale, Assign Beam Weights",

        label = false,
        zcolor = log10.(beam_weights[beams_to_plot]),
        markersize = 6,

        colormap = :ice,
        colorbar_title = "Log10 Electron Count, #",
        clims = count_clims .+ 1,
        leftmargin = 10mm
    )
    p3 = plot!()

    plot_distribution(avg_energies, avg_pitch_angles[.!cull_idxs], corrected_backscatter_counts[:, .!cull_idxs], :ice)
    plot!(
        title = "4) Look Up Backscatter, Azimuth De-Scale",
        colorbar_title = "Log10 Electron Fluence, #",
        clims = count_clims
    )
    p4 = plot!()

    plot(p1, p2, p3, p4,
        layout = (2,2),
        size = (1,.8) .* 800,
        thickness_scaling = .65
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/simulation_procedure.png")
end

function figure_elfin_backscatter_vs_nflux()
    data = npzread("$(TOP_LEVEL)/data/figure_data/ELFIN_lifetime_backscatter.npz")
    nflux_bin_edges = data["nflux_bin_edges"]
    backscatter_bin_edges = data["backscatter_bin_edges"]

    n_frequencies = data["n_frequencies"]
    normalized_n_frequencies = data["normalized_n_frequencies"]

    e_frequencies = data["e_frequencies"]
    normalized_e_frequencies = data["normalized_e_frequencies"]

    backscatter_heatmap_plot()
    heatmap!(log10.(nflux_bin_edges), backscatter_bin_edges, log10.(n_frequencies)',
        title = "ELFIN Number Backscatter",

        ylabel = "Number Backscatter Ratio, # ALC / # LC",
        
        colorbar = false,
        colorbar_title = "Log10 Occurrences, # Data Segments",
        colormap = :ice,
        #clims = (0, 2),

        background_color_inside = :black
    )
    p1 = plot!()

    backscatter_heatmap_plot()
    heatmap!(log10.(nflux_bin_edges), backscatter_bin_edges, log10.(e_frequencies)',
        title = "ELFIN Energy Backscatter",

        ylabel = "Energy Backscatter Ratio, keV ALC / keV LC",
        
        colorbar = true,
        colorbar_title = "Log10 Occurrences, # Data Segments",
        colormap = :ice,
        #clims = (0, 2),

        background_color_inside = :black
    )
    p2 = plot!()

    layout = @layout [a{.438w} b]
    plot(p1, p2,
        layout = layout,
        size = (2, .85) .* 400,
        dpi = 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/elfin_backscatter_vs_nflux.png")
end

function figure_elfin_backscatter_vs_precipitation_ratio()
    data = npzread("$(TOP_LEVEL)/data/figure_data/precipitation_ratio_vs_backscatter.npz")

    JoverJ90_bin_edges = data["JoverJ90_bin_edges"]
    backscatter_bin_edges = data["backscatter_bin_edges"]

    energy_frequencies = data["energy_frequencies"]
    number_frequencies = data["number_frequencies"]

    trimmed_energy_frequencies = data["trimmed_energy_frequencies"]
    trimmed_number_frequencies = data["trimmed_number_frequencies"]

    clims = (0, 2.3)
    heatmap(log10.(JoverJ90_bin_edges), backscatter_bin_edges, log10.(number_frequencies'),
        title = "All Eligible Data",
    
        xlabel = "J precip. / J trap",
        xlims = (-2.5, 1),
        xticks = (-2:1, ["10⁻²", "10⁻¹", "10⁰", "10¹"]),
        xminorticks = 4,

        ylabel = "Number Backscatter Ratio",
        ylims = (0, 1),
        yticks = 0:.1:1.1,
        yminorticks = 5,

        colorbar_title = "Log10 Occurrences, # Data Segments",
        colormap = :ice,
        clims = clims,
        colorbar = false,

        leftmargin = 10mm,
        rightmargin = 5mm,

        framestyle = :box,
        tickdirection = :out,
        background_color_inside = :black
    )
    box_aspect!(1)
    p1 = plot!()

    heatmap(log10.(JoverJ90_bin_edges), backscatter_bin_edges, log10.(trimmed_number_frequencies'),
        title = "Low Precipitating Fluxes Excluded",
        
        xlabel = "J precip. / J trap",
        xlims = (-2.5, 1),
        xticks = (-2:1, ["10⁻²", "10⁻¹", "10⁰", "10¹"]),
        xminorticks = 4,

        ylabel = "",
        ylims = (0, 1),
        yticks = 0:.1:1,
        yminorticks = 5,

        colorbar_title = "Log10 Occurrences, # Data Segments",
        colormap = :ice,
        clims = clims,

        leftmargin = 5mm,
        rightmargin = 5mm,

        framestyle = :box,
        tickdirection = :out,
        background_color_inside = :black
    )
    box_aspect!(1)
    p2 = plot!()

    layout = @layout [a{.438w} b]
    plot(p1, p2,
        layout = layout,
        size = (2, .9) .* 400,
        background_color_outside = :transparent,
        thickness_scaling = .6,
        dpi = 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/backscatter_vs_precipitation_ratio.png")
end

function figure_g4epp_backscatter_vs_nflux()
    data = npzread("$(TOP_LEVEL)/data/figure_data/ELFIN_lifetime_backscatter.npz")

    nflux_bin_edges = data["nflux_bin_edges"]
    backscatter_bin_edges = data["backscatter_bin_edges"]

    n_frequencies = data["n_frequencies"]
    e_frequencies = data["e_frequencies"]

    sim_n_frequencies = data["sim_n_frequencies"]
    sim_e_frequencies = data["sim_e_frequencies"]


    backscatter_heatmap_plot()
    heatmap!(log10.(nflux_bin_edges), backscatter_bin_edges, log10.(sim_n_frequencies)',
        title = "Simulation Number Backscatter",

        ylabel = "Backscatter Ratio, # ALC / # LC",
        
        colorbar = false,
        colorbar_title = "Log10 Density",
        colormap = :navia,
        clims = (1, 3)
    )
    p1 = plot!()

    backscatter_heatmap_plot()
    heatmap!(log10.(nflux_bin_edges), backscatter_bin_edges, log10.(sim_e_frequencies)',
        title = "Simulation Energy Backscatter",

        ylabel = "Backscatter Ratio, keV ALC / keV LC",
        
        colorbar = true,
        colorbar_title = "Log10 Density",
        colormap = :navia,
        clims = (1, 3)
    )
    p2 = plot!()

    layout = @layout [a{.4425w} b]
    plot(p1, p2,
        layout = layout,
        size = (2, .85) .* 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/g4epp_backscatter_vs_nflux.png")
end

function figure_1_to_1_event_simulation_residuals()
    # Load data
    data = npzread("$(TOP_LEVEL)/data/figure_data/g4epp_elfin_comparison.npz")
    residual_bin_edges = data["residual_bin_edges"]
    number_residual_frequencies = data["number_residual_frequencies"]
    energy_residual_frequencies = data["energy_residual_frequencies"]
    
    # Plot
    residual_1d_histogram(residual_bin_edges, number_residual_frequencies)
    plot!(xlabel = "Simulation Backscatter/Data Backscatter, #/#")
    p1 = plot!()
    
    residual_1d_histogram(residual_bin_edges, energy_residual_frequencies)
    plot!(xlabel = "Simulation Backscatter/Data Backscatter, keV/keV")
    p2 = plot!()
    
    plot(p1, p2,
        layout = (2,1),
        size = (1, .65) .* 500,
        thickness_scaling = .6,
        background = :transparent,
        leftmargin = 10mm,
        dpi = 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/backscatter_residuals.png")
end

function figure_predicted_backscatter()
    data = npzread("$(TOP_LEVEL)/data/figure_data/simulated_backscatter_rates.npz")
    energies = data["energies"]
    energy_bin_edges = data["energy_bin_edges"]
    coarse_pitch_angles = data["coarse_pitch_angles"]
    fine_pitch_angles = data["fine_pitch_angles"]

    energy_backscatter_coarse_grid = data["energy_backscatter_coarse_grid"]
    number_backscatter_coarse_grid = data["number_backscatter_coarse_grid"]

    energy_backscatter_fine_grid = data["energy_backscatter_fine_grid"]
    number_backscatter_fine_grid = data["number_backscatter_fine_grid"]

    α_to_plot = -5:10:95
    predicted_backscatter_plot(α_to_plot, number_backscatter_coarse_grid)
    p1 = plot!(
        title = "Simulated Number Backscatter",
        xlabel = "",
        xticks = 0:10:90,
    )
    predicted_backscatter_plot(α_to_plot, energy_backscatter_coarse_grid)
    p2 = plot!(
        title = "Simulated Energy Backscatter",
        xlabel = "",
        ylabel = "",
        xticks = 0:10:90,
    )

    α_to_plot = 59.5:1:70.5
    predicted_backscatter_plot(α_to_plot, number_backscatter_fine_grid)
    p3 = plot!(
        xticks = 60:1:70
    )
    predicted_backscatter_plot(α_to_plot, energy_backscatter_fine_grid)
    p4 = plot!(
        ylabel = "",
        xticks = 60:1:70
    )

    plot(p1, p2, p3, p4,
        layout = (2,2),
        size = (1.5, 1.5) .* 800,
        dpi = 300,
        background_color_outside = :transparent
    )
    main_plot = plot!()

    heatmap(rand(2,2) .* NaN,
        colorbar_title = "Backscatter Ratio, %",
        colormap = :navia,
        colorbar = true,
        clims = (0, 100),

        leftmargin = -5mm,  

        axis = false,
        grid = false,
        size = (.2, 1) .* 600,
        dpi = 400,
        background_color = :transparent
    )
    colorbar_plot = plot!()

    layout = @layout [a{.88w} b]
    plot(main_plot, colorbar_plot,
        layout = layout,
        size = (1.1, 1) .* 700,
        thickness_scaling = .6,
        background_color = :transparent,
        dpi = 400
    )

    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/predicted_backscatter.png")
end

function figure_backscattered_pads()
    data = npzread("$(TOP_LEVEL)/data/figure_data/backscattered_pads.npz")
    pad_bin_edges = data["pad_bin_edges"]
    input_pa_edges = data["input_pa_edges"]
    low_energy = data["low_energy"]
    low_energy_pad_weights = data["low_energy_pad_weights"]
    high_energy = data["high_energy"]
    high_energy_pad_weights = data["high_energy_pad_weights"]

    backscattered_pad(pad_bin_edges, input_pa_edges, low_energy_pad_weights, low_energy)
    p1 = plot!(colorbar = false, leftmargin = 5mm)

    backscattered_pad(pad_bin_edges, input_pa_edges, high_energy_pad_weights, high_energy)
    p2 = plot!(ylabel = "")

    layout = @layout [a{.44w} b]
    plot(p1, p2,
        layout = layout,
        size = (2, .9) .* 600,
        dpi = 400
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/backscattered_pads.png")
end

function supplemental_ribbon()
    data = npzread("$(TOP_LEVEL)/data/figure_data/supplemental_curvefit.npz")

    ribbon_plot(data, data["number_25"], data["number_50"], data["number_75"])
    plot!(ylabel = "Number Backscatter Ratio")
    p1 = plot!()

    ribbon_plot(data, data["energy_25"], data["energy_50"], data["energy_75"])
    plot!(ylabel = "Energy Backscatter Ratio")
    p2 = plot!()

    plot(p1, p2,
        layout = (1,2),
        size = (1.8, .9) .* 500,
        dpi = 400,
        background = :transparent
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/supplemental_curves")
end

function supplemental_curvefit()
    data = npzread("$(TOP_LEVEL)/data/figure_data/supplemental_curvefit.npz")
    x = edges_to_means(data["JoverJ90_bin_edges"])

    plot()
    plot_precipitation_curve!(log10.(x), data["number_50"])
    plot!(log10.(x), data["number_a"] .* (x.^data["number_b"]),
        label = "Fit",
        linewidth = 1.5,
        linecolor = :black,
        ylabel = "Number Backscatter Ratio"
    )
    annotate!(-.5, .95, text(latexstring("R^2 = $(round(data["number_r2"], sigdigits = 4))"), :center))
    p1 = plot!()

    plot()
    plot_precipitation_curve!(log10.(x), data["energy_50"])
    plot!(log10.(x), data["energy_a"] .* (x.^data["energy_b"]),
        label = "Fit",
        linewidth = 1.5,
        linecolor = :black,
        ylabel = "Energy Backscatter Ratio"
    )
    annotate!(-.5, .95, text(latexstring("R^2 = $(round(data["energy_r2"], sigdigits = 4))"), :center))
    p2 = plot!()

    plot(p1, p2,
        layout = (1,2),
        size = (1.8, .9) .* 500,
        dpi = 400,
        background = :transparent
    )
    display(plot!())
    png("$(TOP_LEVEL)/paper/figures/supplemental_curvefit")
end

"""
======================================
Recipes
======================================
"""
function backscatter_heatmap_plot()
    plot(
        title = "",

        xlabel = "Loss Cone Number Flux, # cm⁻² s⁻¹",
        xlims = (1, 7),
        xminorticks = 4,
        xticks = (1:7, ["10¹", "10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"]),

        ylabel = "",
        ylims = (0, 1),
        yticks = 0:.1:1,
        yminorticks = 5
    )
    plot!(
        size = (1, .65) .* 400,
        topmargin = 0mm,
        leftmargin = 15mm,
        thickness_scaling = .6,

        framestyle = :box,
        tickdirection = :out,
        minorgrid = true,

        background_color_outside = :transparent,
        background_color_inside = :black,

        dpi = 400
    )
    box_aspect!(1)
    return plot!()
end

function plot_1d_histogram!(bin_edges, data;
    linecolor = :black,
    linestyle = :solid,
    fillcolor = RGBA(0,0,0,.1),
    fillstyle = nothing,
    label = false
    )

    frequencies = exact_1dhistogram(data, bin_edges)
    plot!(bin_edges, [frequencies..., frequencies[end]],
        label = label,
        linestyle = linestyle,
        linetype = :steppost,
        linecolor = linecolor,
        fill = true,
        fillcolor = fillcolor,
        fillstyle = fillstyle
    )
    return plot!()
end

function residual_1d_histogram(bin_edges, frequencies)
    plot(bin_edges, [frequencies..., frequencies[end]],
        label = false,
        linetype = :steppost,
        linecolor = :black,
        fill = true,
        fillcolor = RGBA(0,0,0,.1),

        xlims = (10^-1.5, 10^1.5),
        xticks = 10.0.^(-2:1:2),
        xminorticks = true,
        xminorgrid = true,
        xscale = :log10,

        ylabel = "Occurrences, # Data Segments",
        ylims = (0, 1.05max(frequencies...)),

        leftmargin = 5mm,

        framestyle = :box,
        tickdirection = :out
    )
    box_aspect!(.25)
    return plot!()
end

function plot_event(event::Event)
    xtick_idxs = Int.(round.(LinRange(1, event.n_datapoints, 5)))
    xtick_labels = [Dates.format.(Time.(event.time_datetime[i]), "HH:MM:SS") * "\nL = $(round(event.L[i], digits = 1))" for i in xtick_idxs]

    e_flux, n_flux = integrate_flux(event, pitch_angle = true)
    heatmap(event.time, event.energy_bins_mean, log10.(n_flux'),
        title = "$(Date(event.time_datetime[1])) ELFIN-$(event.satellite)",
        
        xlims = (event.time[begin], event.time[end]),
        xticks = (event.time[xtick_idxs], xtick_labels),
        xminorticks = true,

        ylabel = "Energy, keV",
        ylims = (event.energy_bins_min[begin], event.energy_bins_max[end]),
        yscale = :log10,
        yminorticks = true,

        colorbar_title = "Log10 Electron Flux\n# s⁻¹ cm⁻² MeV⁻¹",
        colormap = :ice,
        clims = (3, 6.5),

        leftmargin = 8mm,
        rightmargin = 5mm,
        topmargin = -5mm,
        
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :box,
        tickdirection = :out,
        thickness_scaling = .8,
        size = (2, .5) .* 350,
        dpi = 400
    )
    box_aspect!(.2)
    return plot!()
end

function pad_heatmap(event::Event, t::Int)
    pad = event.n_flux[t,:,:]
    α_lc = event.loss_cone_angles[t]
    α_alc = event.anti_loss_cone_angles[t]

    generate_pad_plot(event.avg_pitch_angles, pad, α_lc, α_alc, 0, annotate = false)
end

function plot_distribution(e_means, pa_means, values, colormap)
    event = example_event()
    heatmap(pa_means, e_means, log10.(values),
        xlabel = "Pitch Angle, deg",
        xlims = (0, 180),
        xticks = 0:30:180,

        ylabel = "Energy, keV",
        ylims = (event.energy_bins_min[begin], event.energy_bins_max[end]),
        yscale = :log10,
        yminorticks = true,

        colormap = colormap,

        topmargin = -5mm,
        bottommargin = -5mm,
        leftmargin = 5mm,
        rightmargin = 5mm,

        background = :white,
        grid = false,
        framestyle = :box,
        tickdirection = :out,
        thickness_scaling = .8,
        size = (1.3, 1) .* 400,
        dpi = 400
    )
    vline!([67],
        label = false,
        linecolor = RGB(0xff006a),
        style = :solid,
        linewidth = 2.5
    )
    vline!([113],
        label = false,
        linecolor = RGB(0xff006a),
        style = :dash,
        linewidth = 2.5
    )
    box_aspect!(1)
end

function predicted_backscatter_plot(α_to_plot, backscatter)
    data = npzread("$(TOP_LEVEL)/data/figure_data/simulated_backscatter_rates.npz")
    energies = data["energies"]
    energy_bin_edges = data["energy_bin_edges"]
    pitch_angles = edges_to_means(α_to_plot)

    heatmap(α_to_plot, energy_bin_edges, backscatter .* 100,
        xlabel = "Input Pitch Angle, deg",
        xlims = (min(α_to_plot...), max(α_to_plot...)),

        ylabel = "Input Energy, keV",
        yscale = :log10,
        ylims = (energy_bin_edges[begin], energy_bin_edges[end]),
        yminorticks = true,

        colorbar = false,
        clims = (0, 100),
        colormap = :navia, #:navia,

        leftmargin = 8mm,

        background_color_outside = :transparent,
        background_color_inside = :black,
        framestyle = :box,
        tickdirection = :out,
        size = (1,1) .* 600
    )
    for E in eachindex(energies)
        for α in eachindex(pitch_angles)
            color = RGBA(1,1,1,.75)
            if backscatter[E, α] > .7; color = RGBA(.0,.1,.1,.75); end
            annotate!(pitch_angles[α], energies[E], text("$(round(backscatter[E,α] * 100, digits = 1))%", 7, color, :center))
        end
    end
    box_aspect!(1)
    return plot!()
end

function backscattered_pad(pad_bin_edges, input_pa_edges, pad_weights, input_energy_kev)
    heatmap(pad_bin_edges, input_pa_edges, log10.(pad_weights),
        title = "$(Int(round(input_energy_kev))) keV Input",

        xlabel = "Backscattered Pitch Angle, deg",
        xlims = (90, 180),
        xticks = 90:10:180,
        xminorticks = 5,

        ylabel = "Input Pitch Angle, deg",
        ylims = (0, 90),
        yticks = 0:10:90,
        yminorticks = 5,

        colorbar_title = "\nLog10 Pitch Angle Density, electrons/str • deg",
        clims = (-2, -1),
        colormap = :navia,

        background_color_outside = :transparent,
        background_color_inside = cgrad(:navia)[begin],
        framestyle = :box,
        tickdirection = :out,

        rightmargin = 5mm,
        
        legendposition = :topright,
        background_color_legend = :black,
        foreground_color_legend = :white,
        legendfontcolor = :white,

        size = (1.15,1) .* 600,
        dpi = 400
    )
    vline!([180 - 67],
        label = L"$\alpha_{alc}, \alpha_{lc}$",

        linestyle = :dash,
        color = :white,
    )
    hline!([67],
        label = false,

        linestyle = :dash,
        color = :white,
    )
    plot!([90, 180], [90, 0],
        label = L"$\alpha_{in} = \alpha_{out}$",

        linestyle = :solid,
        color = :white,
    )
    annotate!(140, 80, "Untrapped", :white)
    annotate!(100, 30, text("Retrapped", :white, rotation = 90))
    box_aspect!(1)
    return plot!()
end

function plot_ring!(L)
    θ = 0:.01:2π
    plot!(θ, (θ.*0) .+ L,
        label = false,
        linecolor = :white,
        linealpha = 0.2
    )
end

function plot_mlt!(MLT)
    θ = MLT * 2π/24
    plot!([θ, θ], [0, Plots.ylims()[2]],
        label = false,
        linecolor = :white,
        linealpha = 0.2
    )
end

function plot_earth!()
    plot!(Shape(vcat(Plots.partialcircle(0, 2π, 100, 0), reverse(Plots.partialcircle(0, 2π, 100, 1/Plots.ylims()[2])))),
    fillcolor = :white,
    linecolor = RGBA(0,0,0,0),
    label = ""
    )
    plot!(Shape(vcat(Plots.partialcircle(-π/2, π/2, 100, 0), reverse(Plots.partialcircle(-π/2, π/2, 100, 1/(Plots.ylims()[2]+2) )))),
        fillcolor = :black,
        linecolor = RGBA(0,0,0,0),
        label = ""
    )
end

function generate_pad_plot(pitch_angles, pad, α_lc, α_alc, rn; annotate = false)
    upper_clim_Δ = -.1
    lower_clim_Δ = -2.5
    heatmap(pitch_angles, example_event().energy_bins_mean, log10.(pad),
        xlabel = "Pitch Angle, deg",
        xlims = (0, 180),
        xticks = 0:30:180,

        ylims = (example_event().energy_bins_min[begin], example_event().energy_bins_max[end]),
        yscale = :log10,
        yminorticks = true,

        colorbar_title = "Log10 Electron Flux, # s⁻¹ cm⁻² str⁻¹ MeV⁻¹",
        colormap = :ice,
        clims = (log10(max(pad...)) + lower_clim_Δ, log10(max(pad...))+upper_clim_Δ),

        topmargin = -5mm,

        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :box,
        tickdirection = :out,
        thickness_scaling = .8,
        size = (1.3, 1) .* 400,
        dpi = 400
    )
    vline!([α_lc],
        label = false,
        linewidth = 2.5,
        linecolor = RGB(0xff006a)
    )
    vline!([α_alc],
        label = false,
        linestyle = :dash,
        linewidth = 2.5,
        linecolor = RGB(0xff006a)
    )
    if annotate  == true
        annotation = latexstring("r_N = $(round(rn, digits = 2))")
        annotate!(175, 6000, text(annotation, :white, :right, :top, 20))
    end
    box_aspect!(1)
    return plot!()
end

function ribbon_plot(data, lower, middle, upper)
    x = log10.(edges_to_means(data["JoverJ90_bin_edges"]))
    plot(x, middle,
        label = false,
        linecolor = :black,
        ribbon = (middle .- lower, upper .- middle),
        fillcolor = RGB(0,0,0),
        fillalpha = 0.2,

        marker = true,
        markercolor = :black,
        markersize = 3,
    
        xlabel = "J precip. / J trap",
        xlims = (-2, 1),
        xticks = (-2:1, ["10⁻²", "10⁻¹", "10⁰", "10¹"]),
        xminorticks = 4,

        ylims = (0, 1),
        yticks = 0:.1:1,

        framestyle = :box,
        tickdirection = :out,
        margin = 5mm
    )
    box_aspect!(1)
    return plot!()
end

function plot_precipitation_curve!(x, y)
    plot!(x, y,
        label = "Data",
        linecolor = RGBA(0,0,0),
        linealpha = 0.3,
        linewidth = 1.5,

        marker = true,
        markercolor = RGB(0,0,0),
        markerstrokecolor = RGB(0,0,0),
        markeralpha = 0.3,
        markersize = 3,
    
        xlabel = "J precip. / J trap",
        xlims = (-2, 1),
        xticks = (-2:1, ["10⁻²", "10⁻¹", "10⁰", "10¹"]),
        xminorticks = 4,

        ylims = (0, 1),
        yticks = 0:.1:1,

        framestyle = :box,
        tickdirection = :out,
        margin = 5mm,
        legendposition = :topright
    )
    box_aspect!(1)
    return plot!()
end


"""
======================================
Figure Generation
======================================
"""

#=
figure_elfin_data()
figure_data_coverage()

figure_elfin_backscatter_vs_nflux()
figure_elfin_backscatter_vs_precipitation_ratio()
figure_filling_vs_backscatter()

figure_beam_weighting()
figure_1_to_1_event_simulation_residuals()
figure_predicted_backscatter()
figure_backscattered_pads()

supplemental_ribbon()
supplemental_curvefit()

error()
=#



function meow()
    backscatter_ratio_maximum_absolute_error = .025

    data = readdlm("$(TOP_LEVEL)/data/ELFIN_backscatter_and_simulation_no_azimuth_correction_half_fov_standoff.csv", ',', skipstart = 1)
    
    # Event metadata
    start = DateTime.(data[:,1])
    stop = DateTime.(data[:,2])
    sat = data[:,3]
    Δt = [(stop[i] - start[i]).value ./ 1000 for i in eachindex(start)]

    L_start = data[:,4]
    L_stop = data[:,5]
    MLT_start = data[:,6]
    MLT_stop = data[:,7]

    # Data-derived quantities
    lc_n_sectors = data[:,19]
    alc_n_sectors = data[:,20]

    loss_cone_eflux = data[:,8] ./ ((lc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)
    loss_cone_nflux = data[:,9] ./ ((lc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)

    trapped_eflux = data[:,10] ./ (( (16 .-lc_n_sectors .-alc_n_sectors).*Δt./16) .* ELFIN_EPD_AREA) 
    trapped_nflux = data[:,11] ./ (( (16 .-lc_n_sectors .-alc_n_sectors).*Δt./16) .* ELFIN_EPD_AREA) 

    anti_loss_cone_eflux = data[:,12] ./ ((alc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)
    anti_loss_cone_nflux = data[:,13] ./ ((alc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)

    # Simulation quantities
    sim_anti_loss_cone_eflux = data[:,17] ./ ((alc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)
    sim_anti_loss_cone_nflux = data[:,18] ./ ((alc_n_sectors.*Δt./16) .* ELFIN_EPD_AREA)

    # Derived quantities
    energy_backscatter_ratio = anti_loss_cone_eflux ./ loss_cone_eflux
    number_backscatter_ratio = anti_loss_cone_nflux ./ loss_cone_nflux

    sim_energy_backscatter_ratio = sim_anti_loss_cone_eflux ./ loss_cone_eflux
    sim_number_backscatter_ratio = sim_anti_loss_cone_nflux ./ loss_cone_nflux

    JoverJ90 = loss_cone_nflux ./ trapped_nflux

    # Error data
    lc_relative_error = data[:,14]
    alc_relative_error = data[:,16]

    uncertainty_energy_backscatter_ratio = abs.(energy_backscatter_ratio) .* sqrt.(alc_relative_error.^2 .+ lc_relative_error.^2)
    uncertainty_number_backscatter_ratio = abs.(number_backscatter_ratio) .* sqrt.(alc_relative_error.^2 .+ lc_relative_error.^2)

    bad_data_idxs = (-0.5 .< log10.(JoverJ90)) .&& (0.9 .< number_backscatter_ratio)
    good_data_idxs = .!bad_data_idxs

    # Filter the data to only what we want to analyze
    e_idxs_to_analyze = (
        (uncertainty_energy_backscatter_ratio .< backscatter_ratio_maximum_absolute_error) 
        .&& good_data_idxs
        .&& (loss_cone_eflux .> 1)
    )
    
    n_idxs_to_analyze = (
        (uncertainty_number_backscatter_ratio .< backscatter_ratio_maximum_absolute_error) 
        .&& good_data_idxs
    )

    avg_energy = loss_cone_eflux ./ loss_cone_nflux


    data_alc_number = data[:,13]
    sim_alc_number = data[:,18]

    data_alc_flux_edges = 0:.1:10
    sim_alc_flux_edges = 0:.1:10.5

    f = exact_2dhistogram(log10.(data_alc_number[n_idxs_to_analyze]), log10.(sim_alc_number[n_idxs_to_analyze]), data_alc_flux_edges, sim_alc_flux_edges)
    heatmap(data_alc_flux_edges, sim_alc_flux_edges, log10.(f'),
        bg = :black,
        aspect_ratio = 1,
        xlabel = "Data ALC #",
        xlims = (0, 6.5),
        ylabel = "Sim ALC #",
        ylims = (0, 6.5),
    )
    plot!([-100, 100], [-100, 100],
        label = false,
        linecolor = :white,
        linestyle = :dash,
        linewidth = 1.5
    )
    annotate!((3, 1, text("Underestimation", :white, :left)))
    annotate!((1, 5.5, text("Overestimation", :white, :left)))
    display(plot!())

    error()


    data_backscatter_edges = 0:.025:1.5
    sim_backscatter_edges = 0:.025:2
    frequencies = exact_2dhistogram(number_backscatter_ratio[n_idxs_to_analyze], sim_number_backscatter_ratio[n_idxs_to_analyze], data_backscatter_edges, sim_backscatter_edges)

    heatmap(data_backscatter_edges, sim_backscatter_edges, log10.(frequencies'),
        xlabel = "Data r_N",
        xlims = (0, 1),
        xticks = 0:.1:1,

        ylabel = "Sim r_N",
        ylims = (0, 1),
        yticks = 0:.1:1,

        aspect_ratio = 1,
        background_color_inside = :black
    )
    plot!([0, 2], [0, 2],
        label = false,
        linecolor = :white,
        linestyle = :dash,
        linewidth = 1.5
    )
    annotate!((0.5, 0.1, text("Underestimation", :white, :left)))
    annotate!((0.2, 0.9, text("Overestimation", :white, :left)))


    display(plot!())
    error()

    
    red_idxs = (alc_n_sectors .== 6) .&& n_idxs_to_analyze
    blue_idxs = (.!red_idxs) .&& n_idxs_to_analyze

    scatter!([number_backscatter_ratio[red_idxs]], [sim_number_backscatter_ratio[red_idxs]],
        label = false,
        markercolor = :red,
        dpi = 300
    )
    display("image/png", plot!())

    #=
    idxs = n_idxs_to_analyze .&& (sim_number_backscatter_ratio .> 0.3)
    idxs = findall(idxs)
    for i in eachindex(idxs)
        print_progress_bar(i/length(idxs))
        northern_hemisphere = create_event(start[idxs[i]], stop[idxs[i]], sat[idxs[i]], maximum_relative_error = 0.5).avg_loss_cone_angle < 90


        if northern_hemisphere == false
            continue
        end


        scatter!([number_backscatter_ratio[idxs[i]]], [sim_number_backscatter_ratio[idxs[i]]],
            label = false,
            markercolor = :red,
            dpi = 300
        )
        display("image/png", plot!())
    end
    =#

    error()
    to_view = findall(n_idxs_to_analyze .&& (0.2 .< sim_number_backscatter_ratio .< 0.5) .&& (0.4 .< number_backscatter_ratio .< 0.6))
    display(to_view)

    for idx in to_view
        quicklook(create_event(start[idx], stop[idx], sat[idx], maximum_relative_error = 0.5))
    end


end

meow()