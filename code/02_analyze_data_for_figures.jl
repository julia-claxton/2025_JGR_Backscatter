const TOP_LEVEL = dirname(@__DIR__)

# Library includes
using Statistics, LinearAlgebra
using BenchmarkTools, Profile, TickTock
using Plots, Plots.PlotMeasures
using DelimitedFiles
using Dates
using NumericalIntegration
using NPZ
using CurveFit

# Import my tools
include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Events.jl")
include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Visualization.jl")
include("$(TOP_LEVEL)/code/General_Functions.jl")
include("$(TOP_LEVEL)/code/G4EPP_2.0/Frontend_Functions.jl")

const RESULTS_FILENAME = "ELFIN_backscatter_and_simulation.csv"

function save_backscatter_figure_data(;
    backscatter_ratio_maximum_absolute_error = .025
    )
    print("Saving backscatter statistics... ")
    # Load data
    data = readdlm("$(TOP_LEVEL)/data/$(RESULTS_FILENAME)", ',', skipstart = 1)
    
    # Event metadata
    start = DateTime.(data[:,1])
    stop  = DateTime.(data[:,2])
    sat   = data[:,3]
    Δt    = [(stop[i] - start[i]).value ./ 1000 for i in eachindex(start)]

    L_start = data[:,4]
    L_stop  = data[:,5]

    MLT_start = data[:,6]
    MLT_stop  = data[:,7]

    # Data-derived quantities
    lc_n_sectors   = data[:,19]
    alc_n_sectors  = data[:,20]
    trap_n_sectors = 16 .- lc_n_sectors .- alc_n_sectors

    loss_cone_energy      = data[:,8]
    loss_cone_number      = data[:,9]
    trapped_energy        = data[:,10]
    trapped_number        = data[:,11]
    anti_loss_cone_energy = data[:,12]
    anti_loss_cone_number = data[:,13]

    loss_cone_eflux      =      loss_cone_energy ./ ((lc_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA)
    loss_cone_nflux      =      loss_cone_number ./ ((lc_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA)
    trapped_eflux        =        trapped_energy ./ ((trap_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA) 
    trapped_nflux        =        trapped_number ./ ((trap_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA) 
    anti_loss_cone_eflux = anti_loss_cone_energy ./ ((alc_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA)
    anti_loss_cone_nflux = anti_loss_cone_number ./ ((alc_n_sectors .*Δt ./ 16) .* ELFIN_EPD_AREA)

    # Simulation quantities
    sim_anti_loss_cone_energy = data[:,17]
    sim_anti_loss_cone_number = data[:,18]

    sim_anti_loss_cone_eflux = sim_anti_loss_cone_energy ./ ((alc_n_sectors .* Δt ./ 16) .* ELFIN_EPD_AREA)
    sim_anti_loss_cone_nflux = sim_anti_loss_cone_number ./ ((alc_n_sectors .* Δt ./ 16) .* ELFIN_EPD_AREA)

    # Derived quantities
    energy_backscatter_ratio = anti_loss_cone_energy ./ loss_cone_energy
    number_backscatter_ratio = anti_loss_cone_number ./ loss_cone_number

    sim_energy_backscatter_ratio = sim_anti_loss_cone_energy ./ loss_cone_energy
    sim_number_backscatter_ratio = sim_anti_loss_cone_number ./ loss_cone_number

    JoverJ90 = loss_cone_number ./ trapped_number

    # Error data
    lc_relative_error = data[:,14]
    alc_relative_error = data[:,16]

    uncertainty_energy_backscatter_ratio = abs.(energy_backscatter_ratio) .* sqrt.(alc_relative_error.^2 .+ lc_relative_error.^2)
    uncertainty_number_backscatter_ratio = abs.(number_backscatter_ratio) .* sqrt.(alc_relative_error.^2 .+ lc_relative_error.^2)

    # Save figure data
    nflux_bin_edges = 10.0.^(-2:.15:9)
    backscatter_bin_spacing = .025
    backscatter_bin_edges = 0:backscatter_bin_spacing:2

    # Find bad data
    bad_data_idxs = (-0.5 .< log10.(JoverJ90)) .&& (0.9 .< number_backscatter_ratio)
    good_data_idxs = .!bad_data_idxs
    println("Removed $(sum(Δt[bad_data_idxs])/60) minutes of bad data.\n")

    # Filter the data to only what we want to analyze
    e_idxs_to_analyze = (
        (uncertainty_energy_backscatter_ratio .< backscatter_ratio_maximum_absolute_error) 
        .&& good_data_idxs
        .&& (loss_cone_eflux .> 1)
    )
    
    n_idxs_to_analyze = (
        (uncertainty_number_backscatter_ratio .< backscatter_ratio_maximum_absolute_error) 
        .&& good_data_idxs
        .&& (loss_cone_nflux .> 1)
    )

    # Save data coverage
    MLT_bin_edges = 0:1:24
    L_bin_edges = 0:1:10
    coverage = exact_2dhistogram(MLT_start[n_idxs_to_analyze], L_start[n_idxs_to_analyze], MLT_bin_edges, L_bin_edges, weights = Δt[n_idxs_to_analyze])

    npzwrite("$(TOP_LEVEL)/data/figure_data/data_coverage.npz",
        MLT_bin_edges = MLT_bin_edges,
        MLT_bin_edges_rad = MLT_bin_edges .* (2π/24),
        L_bin_edges = L_bin_edges,
        coverage = coverage
    )

    # Create 2D histograms of backscatter vs. loss cone flux for data and simulation
    e_frequencies = exact_2dhistogram(loss_cone_nflux[e_idxs_to_analyze], energy_backscatter_ratio[e_idxs_to_analyze], nflux_bin_edges, backscatter_bin_edges)
    n_frequencies = exact_2dhistogram(loss_cone_nflux[n_idxs_to_analyze], number_backscatter_ratio[n_idxs_to_analyze], nflux_bin_edges, backscatter_bin_edges)

    sim_e_frequencies = exact_2dhistogram(loss_cone_nflux[e_idxs_to_analyze], sim_energy_backscatter_ratio[e_idxs_to_analyze], nflux_bin_edges, backscatter_bin_edges)
    sim_n_frequencies = exact_2dhistogram(loss_cone_nflux[n_idxs_to_analyze], sim_number_backscatter_ratio[n_idxs_to_analyze], nflux_bin_edges, backscatter_bin_edges)

    npzwrite("$(TOP_LEVEL)/data/figure_data/ELFIN_lifetime_backscatter.npz",
        backscatter_ratio_maximum_absolute_error = backscatter_ratio_maximum_absolute_error,
        total_time = sum(Δt[n_idxs_to_analyze]),

        nflux_bin_edges = nflux_bin_edges,
        backscatter_bin_edges = backscatter_bin_edges,

        e_frequencies = e_frequencies,
        n_frequencies = n_frequencies,
        sim_e_frequencies = sim_e_frequencies,
        sim_n_frequencies = sim_n_frequencies,
    )

    # Precipitation strength figure
    JoverJ90_bin_edges = 10.0 .^ (-4:.1:1)
    backscatter_bin_edges = -.025:.025:2

    energy_frequencies = exact_2dhistogram(JoverJ90[e_idxs_to_analyze], energy_backscatter_ratio[e_idxs_to_analyze], JoverJ90_bin_edges, backscatter_bin_edges)
    number_frequencies = exact_2dhistogram(JoverJ90[n_idxs_to_analyze], number_backscatter_ratio[n_idxs_to_analyze], JoverJ90_bin_edges, backscatter_bin_edges)
    
    # Cut off low fluxes to prove a point
    flux_cutoff = 10.0 ^ 2.75
    trimmed_e_idxs_to_analyze = e_idxs_to_analyze .&& (loss_cone_nflux .> flux_cutoff)
    trimmed_n_idxs_to_analyze = n_idxs_to_analyze .&& (loss_cone_nflux .> flux_cutoff)

    trimmed_energy_frequencies = exact_2dhistogram(JoverJ90[trimmed_e_idxs_to_analyze], energy_backscatter_ratio[trimmed_e_idxs_to_analyze], JoverJ90_bin_edges, backscatter_bin_edges)
    trimmed_number_frequencies = exact_2dhistogram(JoverJ90[trimmed_n_idxs_to_analyze], number_backscatter_ratio[trimmed_n_idxs_to_analyze], JoverJ90_bin_edges, backscatter_bin_edges)

    # Write data
    npzwrite("$(TOP_LEVEL)/data/figure_data/precipitation_ratio_vs_backscatter.npz",
        backscatter_ratio_maximum_absolute_error = backscatter_ratio_maximum_absolute_error,
        total_time = sum(Δt[n_idxs_to_analyze]),

        JoverJ90_bin_edges = JoverJ90_bin_edges,
        backscatter_bin_edges = backscatter_bin_edges,
        energy_frequencies = energy_frequencies,
        number_frequencies = number_frequencies,
        trimmed_energy_frequencies = trimmed_energy_frequencies,
        trimmed_number_frequencies = trimmed_number_frequencies
    )

    # Data/model Comparison
    number_residuals = sim_anti_loss_cone_number ./ anti_loss_cone_number
    energy_residuals = sim_anti_loss_cone_energy ./ anti_loss_cone_energy

    residual_bin_edges = 10.0.^(-2.5:.05:2.5)
    energy_residual_frequencies = exact_1dhistogram(energy_residuals[e_idxs_to_analyze], residual_bin_edges)
    number_residual_frequencies = exact_1dhistogram(number_residuals[n_idxs_to_analyze], residual_bin_edges)
        
    backscatter_bin_spacing = .025
    backscatter_bin_edges = 0:backscatter_bin_spacing:2
    
    energy_backscatter_ratio_frequencies = exact_1dhistogram(energy_backscatter_ratio[e_idxs_to_analyze], backscatter_bin_edges)
    number_backscatter_ratio_frequencies = exact_1dhistogram(number_backscatter_ratio[n_idxs_to_analyze], backscatter_bin_edges)
    
    sim_energy_backscatter_ratio_frequencies = exact_1dhistogram(sim_energy_backscatter_ratio[e_idxs_to_analyze], backscatter_bin_edges)
    sim_number_backscatter_ratio_frequencies = exact_1dhistogram(sim_number_backscatter_ratio[n_idxs_to_analyze], backscatter_bin_edges)


    println("\nAll data analyzed")
    @show length(Δt) * 3
    @show sum(Δt)/3600
    @show mean(Δt)
    
    println("\nData for bakscatter assessment")
    @show sum(n_idxs_to_analyze) * 3
    @show sum(Δt[n_idxs_to_analyze])/3600
    @show mean(Δt[n_idxs_to_analyze])
    println()


    # Create histogram of data backscatter ratio vs. sim backscatter ratio
    comparison_sim_backscatter_edges = 0:.025:1
    comparison_data_backscatter_edges = 0:.025:1.1
    comparison_number_frequencies = exact_2dhistogram(number_backscatter_ratio[n_idxs_to_analyze], sim_number_backscatter_ratio[n_idxs_to_analyze], comparison_data_backscatter_edges, comparison_sim_backscatter_edges)
    comparison_energy_frequencies = exact_2dhistogram(energy_backscatter_ratio[e_idxs_to_analyze], sim_energy_backscatter_ratio[e_idxs_to_analyze], comparison_data_backscatter_edges, comparison_sim_backscatter_edges)
    
    # Histogram of data ALC flux vs sim ALC flux
    data_alc_flux_edges = -2:.1:10
    sim_alc_flux_edges = -2:.1:10.5

    data_model_heatmap_number_frequencies = exact_2dhistogram(log10.(anti_loss_cone_nflux[n_idxs_to_analyze]), log10.(sim_anti_loss_cone_nflux[n_idxs_to_analyze]), data_alc_flux_edges, sim_alc_flux_edges)
    data_model_heatmap_energy_frequencies = exact_2dhistogram(log10.(anti_loss_cone_eflux[e_idxs_to_analyze]), log10.(sim_anti_loss_cone_eflux[e_idxs_to_analyze]), data_alc_flux_edges, sim_alc_flux_edges)


    npzwrite("$(TOP_LEVEL)/data/figure_data/g4epp_elfin_comparison.npz",
        backscatter_ratio_maximum_absolute_error = backscatter_ratio_maximum_absolute_error,
        total_time = sum(Δt[n_idxs_to_analyze]),

        residual_bin_edges = residual_bin_edges,

        number_residual_frequencies = number_residual_frequencies,
        energy_residual_frequencies = energy_residual_frequencies,


        backscatter_bin_edges = backscatter_bin_edges,
        
        sim_energy_backscatter_ratio = sim_energy_backscatter_ratio[e_idxs_to_analyze],
        sim_number_backscatter_ratio = sim_number_backscatter_ratio[n_idxs_to_analyze],

        energy_backscatter_ratio = energy_backscatter_ratio[e_idxs_to_analyze],
        number_backscatter_ratio = number_backscatter_ratio[n_idxs_to_analyze],

        comparison_sim_backscatter_edges = comparison_sim_backscatter_edges,
        comparison_data_backscatter_edges = comparison_data_backscatter_edges,
        comparison_number_frequencies = comparison_number_frequencies,
        comparison_energy_frequencies = comparison_energy_frequencies,

        data_alc_flux_edges = data_alc_flux_edges,
        sim_alc_flux_edges = sim_alc_flux_edges,
        data_model_heatmap_number_frequencies = data_model_heatmap_number_frequencies,
        data_model_heatmap_energy_frequencies = data_model_heatmap_energy_frequencies
    )
    println("Done")
end

function save_example_events()
    print("Saving example events... ")
    # Get data
    data = readdlm("$(TOP_LEVEL)/data/$(RESULTS_FILENAME)", ',', skipstart = 1)
    start = DateTime.(data[:,1])

    # Get events
    ambient_start = DateTime("2022-02-06T11:56:30.6")
    ambient_stop = DateTime("2022-02-06T11:56:38")
    ambient_sat = "B"
    ambient = create_event(ambient_start, ambient_stop, ambient_sat, maximum_relative_error = 0.5)
    ambient_pad = get_pad(ambient)
    ambient_lc = ambient.avg_loss_cone_angle
    _, ambient_idx = findmin(abs.(ambient_start .- start))
    ambient_rn = data[ambient_idx, 13] / data[ambient_idx, 9] # 13 = alc nflux, 9 = lc nflux
    ambient_pad = get_pad(ambient)

    filling_start = DateTime("2020-08-14T20:25:42")
    filling_stop = DateTime("2020-08-14T20:25:50")
    filling_sat = "B"
    filling = create_event(filling_start, filling_stop, filling_sat, maximum_relative_error = 0.5)
    filling_pad = get_pad(filling)
    filling_lc = filling.avg_loss_cone_angle
    _, filling_idx = findmin(abs.(filling_start .- start))
    filling_rn = data[filling_idx, 13] / data[filling_idx, 9] # 13 = alc nflux, 9 = lc nflux
    filling_pad = get_pad(filling)

    between_start = DateTime("2021-03-06T16:26:35")
    between_stop = DateTime("2021-03-06T16:26:43")
    between_sat = "B"
    between = create_event(between_start, between_stop, between_sat, maximum_relative_error = 0.5)
    between_pad = get_pad(between)
    between_lc = between.avg_loss_cone_angle
    _, between_idx = findmin(abs.(between_start .- start))
    between_rn = data[between_idx, 13] / data[between_idx, 9] # 13 = alc nflux, 9 = lc nflux
    between_pad = get_pad(between)

    npzwrite("$(TOP_LEVEL)/data/figure_data/backscatter_vs_filling.npz",
        energy_means = ambient.energy_bins_mean,

        ambient_pa_means = ambient.avg_pitch_angles,
        ambient_pad = ambient_pad,
        ambient_rn = ambient_rn,
        ambient_lc = ambient_lc,

        filling_pa_means = filling.avg_pitch_angles,
        filling_pad = filling_pad,
        filling_rn = filling_rn,
        filling_lc = filling_lc,

        between_pa_means = between.avg_pitch_angles,
        between_pad = between_pad,
        between_rn = between_rn,
        between_lc = between_lc
    )
    println("Done")
end

function save_g4epp_predictions()
    print("Saving G4EPP backscatter characteristics... ")
    # Get backscatter statistics
    beam_energies, beam_pitch_angles = _get_beam_locations()
    beam_locations = collect(zip(beam_energies, beam_pitch_angles))

    energies = sort(unique(beam_energies))

    # Get large grid (0º to 90º in 10 degree increments)
    coarse_pitch_angles = float.(0:10:90)
    energy_backscatter_coarse_grid = zeros(length(energies), length(coarse_pitch_angles))
    number_backscatter_coarse_grid = zeros(length(energies), length(coarse_pitch_angles))
    for E in eachindex(energies)
        for α in eachindex(coarse_pitch_angles)
            if (energies[E], coarse_pitch_angles[α]) ∉ beam_locations; continue; end
                # Won't need this once all beams have been simulated

            data, n_input_particles, _ = _find_datafile("processed", "backscatter", "electron", energies[E], coarse_pitch_angles[α])
            energy_backscatter_coarse_grid[E,α] = sum(data["electron_energies"]) / (n_input_particles * energies[E])
            number_backscatter_coarse_grid[E,α] = length(data["electron_energies"]) / n_input_particles
        end
    end

    # Get fine grid (60º to 70º 1 degree increments)
    fine_pitch_angles = float.(60:1:70)
    energy_backscatter_fine_grid = zeros(length(energies), length(fine_pitch_angles))
    number_backscatter_fine_grid = zeros(length(energies), length(fine_pitch_angles))
    for E in eachindex(energies)
        for α in eachindex(fine_pitch_angles)
            if (energies[E], fine_pitch_angles[α]) ∉ beam_locations; continue; end
                # Won't need this once all beams have been simulated

            data, n_input_particles, _ = _find_datafile("processed", "backscatter", "electron", energies[E], fine_pitch_angles[α])
            energy_backscatter_fine_grid[E,α] = sum(data["electron_energies"]) / (n_input_particles * energies[E])
            number_backscatter_fine_grid[E,α] = length(data["electron_energies"]) / n_input_particles
        end
    end

    npzwrite("$(TOP_LEVEL)/data/figure_data/simulated_backscatter_rates.npz",
        energies = energies,
        energy_bin_edges = [example_event().energy_bins_min..., example_event().energy_bins_max[end]],
        coarse_pitch_angles = coarse_pitch_angles,
        fine_pitch_angles = fine_pitch_angles,

        energy_backscatter_coarse_grid = energy_backscatter_coarse_grid,
        number_backscatter_coarse_grid = number_backscatter_coarse_grid,

        energy_backscatter_fine_grid = energy_backscatter_fine_grid,
        number_backscatter_fine_grid = number_backscatter_fine_grid
    )
    println("Done")
end

function save_beam_weighting_procedure()
    print("Saving beam weighting data... ")
    # Get event
    event = create_event(DateTime("2021-03-06T07:05:05"), DateTime("2021-03-06T07:05:15"), "B", maximum_relative_error = 0.5)
    ΔE = (event.energy_bins_max .- event.energy_bins_min) ./ 1000 # MeV

    # Get data fluence
    _, data_nfluence = integrate_flux(event, time = true)

    # Convert to counts
    data_counts = [data_nfluence[E,α] .* ΔE[E] .* ELFIN_GEOMETRIC_FACTOR for E in 1:16, α in 1:16]
    backscatter_input_distribution = copy(data_counts)

    # Cull non-downgoing inputs
    lc_idxs = 0 .< event.avg_pitch_angles .< 90
    culled_backscatter_input_distribution = copy(backscatter_input_distribution)
    culled_backscatter_input_distribution[:, .!lc_idxs] .= 0

    # Get beam locations and weights
    beam_energies, beam_pitch_angles = _get_beam_locations()
    distro = create_distribution(event.energy_bins_min, event.energy_bins_max, event.avg_pitch_angles .- (ELFIN_EPD_FOV/2), event.avg_pitch_angles .+ (ELFIN_EPD_FOV/2), culled_backscatter_input_distribution, "counts")
    beam_weights = _counts_to_beam_weights(distro)

    # Get backscatter
    backscatter_counts = simulate_backscatter(distro)

    # Cull parts we aren't comparing
    cull_idxs = event.avg_pitch_angles .< 90

    npzwrite("$(TOP_LEVEL)/data/figure_data/beam_weighting.npz",
        avg_pitch_angles = event.avg_pitch_angles,
        avg_energies = event.energy_bins_mean,
        data_nfluence = data_nfluence,
        data_counts = data_counts,
        backscatter_input_distribution = backscatter_input_distribution,
        culled_backscatter_input_distribution = culled_backscatter_input_distribution,
        beam_energies = beam_energies,
        beam_pitch_angles = beam_pitch_angles,
        beam_weights = beam_weights,
        backscatter_counts = backscatter_counts,
        cull_idxs = cull_idxs
    )
    println("Done")
end

function save_backscattered_pads()
    print("Saving simulated backscatter PADs... ")
    input_pitch_angles = float([0:5:60..., 61:69..., 70:1:90...])
    input_pa_edges = [-2.5:5:57.5..., 60.5:1:69.5..., 70.5:1:90.5...]

    pad_bin_edges = 0:1:180

    # Collect PADs for each energy
    low_energy = float(63)
    high_energy = float(6500)
    low_energy_pad_weights = get_pads_over_input_pitch_angle(input_pitch_angles, pad_bin_edges, low_energy)
    high_energy_pad_weights = get_pads_over_input_pitch_angle(input_pitch_angles, pad_bin_edges, high_energy)

    npzwrite("$(TOP_LEVEL)/data/figure_data/backscattered_pads.npz",
        pad_bin_edges = pad_bin_edges,
        input_pa_edges = input_pa_edges,

        low_energy = low_energy,
        low_energy_pad_weights = low_energy_pad_weights,

        high_energy = high_energy,
        high_energy_pad_weights = high_energy_pad_weights
    )
    println("Done")
end

function save_data_model_events()
    print("Saving data-model comparison events... ")

    names = ["over", "under", "equal"]
    starts = DateTime.(["2022-09-08T04:24:32.827", "2022-03-01T12:20:52.314", "2022-02-23T11:59:46.174"])
    stops = DateTime.(["2022-09-08T04:24:41.689", "2022-03-01T12:21:18.013", "2022-02-23T11:59:54.630"])
    sats = ["A", "B", "A"]

    to_save = Dict{String, Any}()
    for i in eachindex(names)
        event = create_event(starts[i], stops[i], sats[i], maximum_relative_error = 0.5)
        
        # Get data fluence
        ΔE = (event.energy_bins_max .- event.energy_bins_min) ./ 1000 # MeV
        _, data_nfluence = integrate_flux(event, time = true)

        # Convert to counts
        data_counts = [data_nfluence[E,α] .* ΔE[E] .* ELFIN_GEOMETRIC_FACTOR for E in 1:16, α in 1:16]
        backscatter_input_distribution = copy(data_counts)

        # Cull non-downgoing inputs
        lc_idxs = 0 .< event.avg_pitch_angles .< 90
        culled_backscatter_input_distribution = copy(backscatter_input_distribution)
        culled_backscatter_input_distribution[:, .!lc_idxs] .= 0

        # Get beam locations and weights
        beam_energies, beam_pitch_angles = _get_beam_locations()
        distro = create_distribution(event.energy_bins_min, event.energy_bins_max, event.avg_pitch_angles .- (ELFIN_EPD_FOV/2), event.avg_pitch_angles .+ (ELFIN_EPD_FOV/2), culled_backscatter_input_distribution, "counts")
        beam_weights = _counts_to_beam_weights(distro)

        # Get backscatter
        sim_counts = simulate_backscatter(distro)

        # Cull parts we aren't comparing
        cull_idxs = event.avg_pitch_angles .< 90

        # Get alc amounts
        α_alc = event.avg_anti_loss_cone_angle + (ELFIN_EPD_FOV/2)
        alc_idxs = event.avg_pitch_angles .> α_alc

        to_save["$(names[i])_alc_pitch_angle"] = α_alc
        to_save["$(names[i])_alc_idxs"] = alc_idxs
        to_save["$(names[i])_pitch_angles"] = event.avg_pitch_angles
        to_save["$(names[i])_energies"] = event.energy_bins_mean
        to_save["$(names[i])_data_counts"] = data_counts
        to_save["$(names[i])_sim_counts"] = sim_counts

    end
    npzwrite("$(TOP_LEVEL)/data/figure_data/data_model_events.npz", to_save)
    println("Done")
end

function save_supplemental_curvefit_data()
    print("Saving supplemental curvefit... ")
    data = npzread("$(TOP_LEVEL)/data/figure_data/precipitation_ratio_vs_backscatter.npz")

    JoverJ90_bin_edges = data["JoverJ90_bin_edges"]
    backscatter_bin_edges = data["backscatter_bin_edges"]
    backscatter_bin_means = edges_to_means(backscatter_bin_edges)

    trimmed_energy_frequencies = data["trimmed_energy_frequencies"]
    trimmed_number_frequencies = data["trimmed_number_frequencies"]

    # This code stinks but it doesn't need to be that clean. If you're using this in the future, sorry.
    results = Dict{String, Any}()
    results["JoverJ90_bin_edges"] = JoverJ90_bin_edges
    results["backscatter_bin_edges"] = backscatter_bin_edges
    for quantile in [.25, .50, .75]
        n_values = zeros(length(JoverJ90_bin_edges)-1)
        e_values = zeros(length(JoverJ90_bin_edges)-1)

        for i in eachindex(n_values)
            # Number
            n_frequencies = trimmed_number_frequencies[i,:]

            # Don't do anything if there's not enough data to be useful
            if sum(n_frequencies) < 120
                n_values[i] = NaN
                continue
            end

            # Get the index of the quantile and record the correponding precipitation value
            idx = quantile_index(n_frequencies, quantile)
            n_values[i] = backscatter_bin_means[idx]
        end
        results["number_$(Int(round(quantile*100)))"] = n_values

        # I know this is terrible coding practice. However this is not crucial to the analysis, just supplement info.
        for i in eachindex(e_values)
            # Energy
            e_frequencies = trimmed_energy_frequencies[i,:]

            # Don't do anything if there's not enough data to be useful
            if sum(e_frequencies) < 120
                e_values[i] = NaN
                continue
            end

            # Get the index of the quantile and record the correponding precipitation value
            idx = quantile_index(e_frequencies, quantile)
            e_values[i] = backscatter_bin_means[idx]
        end
        results["energy_$(Int(round(quantile*100)))"] = e_values
    end

    # Fit a powerlaw
    JoverJ90_bin_means = edges_to_means(JoverJ90_bin_edges)
    number_idxs = isfinite.(results["number_50"])
    number_a, number_b = power_fit(JoverJ90_bin_means[number_idxs], results["number_50"][number_idxs])
    
    energy_idxs = isfinite.(results["energy_50"])
    energy_a, energy_b = power_fit(JoverJ90_bin_means[energy_idxs], results["energy_50"][energy_idxs])
        # Fits coefficients a and b for y = a .* (x .^ b)

    results["number_a"] = number_a
    results["number_b"] = number_b
    results["energy_a"] = energy_a
    results["energy_b"] = energy_b

    # Get R^2 for these fits
    # Energy
    energy_y = results["energy_50"][energy_idxs]
    energy_f = energy_a .* (JoverJ90_bin_means[energy_idxs] .^ energy_b)

    energy_ss_residual = sum((energy_y .- energy_f).^2)
    energy_ss_total = sum((energy_y .- mean(energy_y)).^2)
    
    results["energy_r2"] = 1 - (energy_ss_residual / energy_ss_total)

    # Number
    number_y = results["number_50"][number_idxs]
    number_f = number_a .* (JoverJ90_bin_means[number_idxs] .^ number_b)

    number_ss_residual = sum((number_y .- number_f).^2)
    number_ss_total = sum((number_y .- mean(number_y)).^2)
    
    results["number_r2"] = 1 - (number_ss_residual / number_ss_total)

    println(
        """
        \n
        r_N &= $(round(number_a, sigdigits = 8)) \\cdot \\left( \\frac{J_\\text{precip}}{J_\\text{trap}} \\right)^{$(round(number_b, sigdigits = 8))} & R^2 = $(round(results["number_r2"], sigdigits = 8)) \\\\
        r_E &= $(round(energy_a, sigdigits = 8)) \\cdot \\left( \\frac{J_\\text{precip}}{J_\\text{trap}} \\right)^{$(round(energy_b, sigdigits = 8))} & R^2 = $(round(results["energy_r2"], sigdigits = 8))
        """
    )

    npzwrite("$(TOP_LEVEL)/data/figure_data/supplemental_curvefit.npz", results)
    println("Done")
end

function get_pads_over_input_pitch_angle(input_pitch_angles, pad_bin_edges, input_energy_kev)
    # Find existing beam locations (can delete when all beams are done)
    beam_energies, beam_pitch_angles = _get_beam_locations()
    beam_locations = collect(zip(beam_energies, beam_pitch_angles))

    Δα = diff(pad_bin_edges)
    ΔΩ = 2π .* (cosd.(pad_bin_edges[begin:end-1]) .- cosd.(pad_bin_edges[begin+1:end]))
    pad_weights = fill(NaN, length(input_pitch_angles), length(pad_bin_edges)-1)
    for i in eachindex(input_pitch_angles)
        α = input_pitch_angles[i]

        # Ensure this beam exists (can delete when all beams are done)
        if (input_energy_kev, α) ∉ beam_locations; continue; end

        # Get backscattered pitch angle distribution
        data, _, _ = _find_datafile("processed", "backscatter", "electron", input_energy_kev, α)
        pad = exact_1dhistogram(data["electron_pitch_angles"], pad_bin_edges)

        # Get in units of str^-1 and normalize to unity
        pad ./= ΔΩ
        pad ./= sum(Δα .* pad)

        # Add to heatmap
        pad_weights[i, :] = pad
    end
    return pad_weights
end

function get_pad(event::Event)
    _, n_fluence = integrate_flux(event, time = true)
    avg_nflux = n_fluence / (event.n_datapoints * ELFIN_SPIN_PERIOD)
end

function quantile_index(v, q)
    normalized_cdf = cumsum(v) ./ sum(v)
    idx = findfirst(normalized_cdf .> q)
    return idx
end

save_backscatter_figure_data()
error()

save_backscatter_figure_data()
save_example_events()
save_data_model_events()
save_g4epp_predictions()
save_beam_weighting_procedure()
save_backscattered_pads()
save_supplemental_curvefit_data()