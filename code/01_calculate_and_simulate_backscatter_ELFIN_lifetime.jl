include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Events.jl")
include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Visualization.jl")
include("/Users/luna/Documents/Work/Research/Backscatter_Analysis/code/General_Functions.jl")
include("/Users/luna/Documents/Work/Research/G4EPP_2.0/Frontend_Functions.jl")

function write_elfin_lifetime_backscatter_data(;
    slice_length_idx = 3,
    result_filename = "ELFIN_backscatter_and_simulation.csv",
    maximum_relative_error = .5,
    )

    # Create results file
    result_path = "$(dirname(@__DIR__))/data/$(result_filename)"
    file = open(result_path, "w")
    write(file, "start_time,stop_time,sat,start_L,stop_L,start_MLT,stop_MLT,lc_energy,lc_number,trap_energy,trap_number,alc_energy,alc_number,lc_relative_error,trap_relative_error,alc_relative_error,simulation_alc_energy,simulation_alc_number\n")
    close(file)

    # Iterate through ELFIN lifetime
    println("Calculating backscatter statistics and simulating events over ELFIN lifetime (Δidx = $(slice_length_idx))")
    dates, sats = all_elfin_science_dates_and_satellite_ids()
    for i in eachindex(dates)
        print_progress_bar(i/length(dates))
        fullday_event = create_event(dates[i], sats[i])
        if fullday_event == nothing; continue; end
        if fullday_event.data_reliable == false; continue; end

        for obs_idx in 1:fullday_event.n_observations
            start_idx = fullday_event.observation_start_idxs[obs_idx]
            stop_idx = fullday_event.observation_stop_idxs[obs_idx]
            event = create_event(fullday_event.time_datetime[start_idx], fullday_event.time_datetime[stop_idx], fullday_event.satellite,
                maximum_relative_error = maximum_relative_error
            )

            if event == nothing; continue; end
            if event.data_reliable == false; continue; end
            if event.n_datapoints < slice_length_idx; continue; end

            start_idxs = 1:slice_length_idx:(event.n_datapoints - slice_length_idx)
            stop_idxs = slice_length_idx:slice_length_idx:event.n_datapoints
            [log_backscatter_data(event, result_path, start_idxs[i], stop_idxs[i], maximum_relative_error) for i in eachindex(start_idxs)]
        end
        if i % 40 == 0; GC.gc(); end
    end
    println()
end

function log_backscatter_data(event::Event, results_path, start_idx, stop_idx, maximum_relative_error)
    # Only accept data with large pitch angle coverage
    worst_pa_min = max(event.pitch_angles[:,begin]...)
    worst_pa_max = min(event.pitch_angles[:,end]...)

    if worst_pa_min > 15; return; end
    if worst_pa_max < (180-15); return; end

    # Reject data segments with multiple spins per timestep
    if (event.time[stop_idx] + ELFIN_SPIN_PERIOD) - event.time[start_idx] > 10; return; end

    # Get loss cone and anti loss cone values
    α_lc = mean(event.loss_cone_angles[start_idx:stop_idx])
    α_alc = mean(event.anti_loss_cone_angles[start_idx:stop_idx])

    # Get loss cone, trapped, and anti-loss cone region bounds
    cone_standoff_angle = 5 + 11.25 # Minimum distance into the loss/antiloss cone a pitch angle bin center must be in order to be counted 
        # ELFIN EPD FOV = 20º x 20º
    if α_lc < 90 
        # Northern hemisphere
        loss_cone_limits = (0, α_lc - cone_standoff_angle)
        trapped_limits = (α_lc, α_alc)
        anti_loss_cone_limits = (α_alc + cone_standoff_angle, 180)
    else 
        # Southern hemisphere
        loss_cone_limits = (α_lc + cone_standoff_angle, 180)
        trapped_limits = (α_alc, α_lc)
        anti_loss_cone_limits = (0, α_alc - cone_standoff_angle)
    end

    # Get data-derived values
    loss_cone_energy, loss_cone_number = integrate_flux(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = loss_cone_limits,
        energy = true
    )
    trapped_energy, trapped_number = integrate_flux(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = trapped_limits,
        energy = true
    )
    anti_loss_cone_energy, anti_loss_cone_number = integrate_flux(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = anti_loss_cone_limits,
        energy = true
    )

    # Get relative errors for each region
    loss_cone_relative_error = relative_error_of_integration(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = loss_cone_limits,
        energy = true
    )
    trapped_relative_error = relative_error_of_integration(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = trapped_limits,
        energy = true
    )
    anti_loss_cone_relative_error = relative_error_of_integration(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = anti_loss_cone_limits,
        energy = true
    )

    # Simulate event with G4EPP
    _, data_nfluence = integrate_flux(event, time = true, time_idxs = start_idx:stop_idx) # number / cm2 str MeV
    ΔE = (event.energy_bins_max .- event.energy_bins_min) ./ 1000 # MeV
    data_counts = [data_nfluence[E,α] .* ΔE[E] .* ELFIN_GEOMETRIC_FACTOR for E in 1:16, α in 1:16] # number
    backscatter_input_distribution = create_distribution(get_elfin_grid_bin_edges(event, time_idxs = start_idx:stop_idx)..., data_counts, "counts")

    # Flip to northern hemisphere for G4EPP if needed
    lc_idxs = loss_cone_limits[1] .< backscatter_input_distribution.pitch_angle_bins_mean .< loss_cone_limits[2]
    alc_idxs = anti_loss_cone_limits[1] .< backscatter_input_distribution.pitch_angle_bins_mean .< anti_loss_cone_limits[2]
    
    if α_lc > 90
        adiabatically_swap_hemisphere!(backscatter_input_distribution)
        lc_idxs = reverse(lc_idxs)
        alc_idxs = reverse(alc_idxs)
    end

    # G4EPP can see all gyrophases, while ELFIN can only see where it was pointed. So, we need to proportionally
    # scale G4EPP by the amount of solid angle in a ring ELFIN could see at a given pitch angle
    Ω_elfin_epd = 2π * (1 - cosd(ELFIN_EPD_FOV/2))
    Ω_G4EPP = 2π * (cosd.(backscatter_input_distribution.pitch_angle_bins_min) .- cosd.(backscatter_input_distribution.pitch_angle_bins_max))
    G4EPP_to_ELFIN_scaling_factor = Ω_elfin_epd ./ Ω_G4EPP

    [backscatter_input_distribution.values[E,:] ./= G4EPP_to_ELFIN_scaling_factor for E in 1:16]

    # Cull inputs outside of the loss cone due to large errors introduced with ELFIN's FOV near the loss cone edge
    # Without this, we get overestimation on the order of 20x. This is probably because backscatter is so incredibly
    # sensitive near the loss cone edge, thus if a 63º electron is snapped to the 68º beam due to bin size, that's a
    # change from 50% backscatter to 100%.
    backscatter_input_distribution.values[:, .!lc_idxs] .= 0

    # Simulate backscatter
    backscatter_output_distribution = copy(backscatter_input_distribution)
    backscatter_output_distribution.values .= 0
    backscatter_output_distribution.values = simulate_backscatter(backscatter_input_distribution)

    # Convert back to ELFIN FOV
    [backscatter_output_distribution.values[E,:] .*= G4EPP_to_ELFIN_scaling_factor for E in 1:16]
    
    # Flip back to correct hemisphere if needed
    if α_lc > 90
        adiabatically_swap_hemisphere!(backscatter_output_distribution)
        lc_idxs = reverse(lc_idxs)
        alc_idxs = reverse(alc_idxs)
    end

    # Cut out bins that would've been discarded in ELFIN data due to low counts
    # Using δq/q ≈ 1/√counts (shot noise), which is what ELFIN does
    sim_relative_error = 1 ./ sqrt.(backscatter_output_distribution.values)
    backscatter_output_distribution.values[sim_relative_error .> maximum_relative_error] .= 0

    # Get simulated backscatter count and energy
    sim_alc_number = sum(backscatter_output_distribution.values[:, alc_idxs])

    #=
    # <debugging>
    if anti_loss_cone_number > 0
        heatmap(log10.(data_counts[:,alc_idxs]), clims = (3,5), bg = :black)
        box_aspect!(2)
        display(plot!())

        heatmap(log10.(backscatter_output_distribution.values[:,alc_idxs]), clims = (3,5), bg = :black)
        box_aspect!(2)
        display(plot!())

        NH = α_lc < 90
        @show NH
        @show sim_alc_number / anti_loss_cone_number
    end
    # </debugging>
    =#


    convert_distribution!(backscatter_output_distribution, "energy")
    sim_alc_energy = sum(backscatter_output_distribution.values[:, alc_idxs])


    # Write results to file
    event_data = "$(event.time_datetime[start_idx]),$(event.time_datetime[stop_idx]+Microsecond(ELFIN_SPIN_PERIOD*1e6)),$(event.satellite),$(event.L[start_idx]),$(event.L[stop_idx]),$(event.MLT[start_idx]),$(event.MLT[stop_idx])"
    lc_data = "$(loss_cone_energy),$(loss_cone_number)"
    trap_data = "$(trapped_energy),$(trapped_number)"
    alc_data = "$(anti_loss_cone_energy),$(anti_loss_cone_number)"
    error_data = "$(loss_cone_relative_error),$(trapped_relative_error),$(anti_loss_cone_relative_error)"
    sim_data = "$(sim_alc_energy),$(sim_alc_number)"

    file = open(results_path, "a")
    write(file, "$(event_data),$(lc_data),$(trap_data),$(alc_data),$(error_data),$(sim_data)\n")
    close(file)
end

function get_elfin_grid_bin_edges(event::Event; time_idxs = 1:event.n_datapoints)
    # Pitch angle. ELFIN EPD nominal FOV = 22.5º per correspondance with E. Tsai
    avg_pitch_angles = dropdims(mean(event.pitch_angles[time_idxs,:], dims = 1), dims = 1)
    elfin_pitch_angle_bins_min = avg_pitch_angles .- (ELFIN_EPD_FOV/2)
    elfin_pitch_angle_bins_max = avg_pitch_angles .+ (ELFIN_EPD_FOV/2)
    
    elfin_pitch_angle_bins_min = clamp.(elfin_pitch_angle_bins_min, 0, 180)
    elfin_pitch_angle_bins_max = clamp.(elfin_pitch_angle_bins_max, 0, 180)

    return event.energy_bins_min, event.energy_bins_max, elfin_pitch_angle_bins_min, elfin_pitch_angle_bins_max
end

write_elfin_lifetime_backscatter_data()