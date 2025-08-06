include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Events.jl")
include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Visualization.jl")
include("/Users/luna/Documents/Work/Research/Backscatter_Analysis/code/General_Functions.jl")
include("/Users/luna/Documents/Work/Research/Backscatter_Analysis/code/Backscatter_Tools.jl")
using Base.Threads
using Statistics, LinearAlgebra
using BenchmarkTools, Profile, TickTock
using BasicInterpolators

# Run with max threads
# julia --threads 10 /Users/luna/Research/Backscatter_Analysis/code/01_calculate_and_simulate_backscatter_ELFIN_lifetime.jl

results_file_threadlock = ReentrantLock()
days_completed_threadlock = ReentrantLock()
progress_bar_threadlock = ReentrantLock()

function write_elfin_lifetime_backscatter_data(;
    slice_length_idx = 3,
    result_filename = "ELFIN_backscatter_and_simulation.csv",
    maximum_relative_error = .5,
    )

    # Create results file
    result_path = "$(dirname(@__DIR__))/data/$(result_filename)"
    file = open(result_path, "w")
    write(file, "start_time,stop_time,sat,start_L,stop_L,start_MLT,stop_MLT,lc_energy,lc_number,trap_energy,trap_number,alc_energy,alc_number,lc_relative_error,trap_relative_error,alc_relative_error,simulation_alc_energy,simulation_alc_number,lc_n_sectors,alc_n_sectors\n")
    close(file)

    # Iterate through ELFIN lifetime
    println("Calculating backscatter statistics and simulating events over ELFIN lifetime (Δidx = $(slice_length_idx))")
    dates, sats = all_elfin_science_dates_and_satellite_ids()
    days_analyzed = 0
    #Threads.@threads for i in eachindex(dates)
    for i in 1000:1100
        lock(progress_bar_threadlock) do
            print_progress_bar(days_analyzed/length(dates), bar_length = 50)
        end

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

        # Updated completion counter
        lock(days_completed_threadlock) do
           days_analyzed += 1 
        end

        if i % 40 == 0; GC.gc(); end
    end
    println()
end

function log_backscatter_data(event::Event, results_path, start_idx, stop_idx, maximum_relative_error)    
    # Only accept data with large pitch angle coverage
    worst_pa_min = max(event.pitch_angles[start_idx:stop_idx, begin]...)
    worst_pa_max = min(event.pitch_angles[start_idx:stop_idx, end]...)

    if worst_pa_min > 15; return; end
    if worst_pa_max < (180-15); return; end

    # Get loss cone and anti loss cone values
    α_lc = mean(event.loss_cone_angles[start_idx:stop_idx])
    α_alc = mean(event.anti_loss_cone_angles[start_idx:stop_idx])
    α_center = dropdims(mean(event.pitch_angles[start_idx:stop_idx, :], dims = 1), dims = 1)

    # Get northern hemisphere equivalent anti-loss cone angle for simulation
    cone_standoff_angle = ELFIN_EPD_FOV / 2 # Minimum distance into the loss/antiloss cone a pitch angle bin center must be in order to be counted 
    α_backscatter = max(α_alc, 180-α_alc) + cone_standoff_angle
    
    # Get loss cone, trapped, and anti-loss cone region bounds
    if α_lc < 90 
        # Northern hemisphere
        loss_cone_limits = (0, α_lc - cone_standoff_angle)
        trapped_limits = (α_lc, α_alc)
        anti_loss_cone_limits = (α_alc + cone_standoff_angle, 180)
        downgoing_limits = (0, 90)
    else 
        # Southern hemisphere
        loss_cone_limits = (α_lc + cone_standoff_angle, 180)
        trapped_limits = (α_alc, α_lc)
        anti_loss_cone_limits = (0, α_alc - cone_standoff_angle)
        downgoing_limits = (90, 180)
    end

    # Get data-derived values
    loss_cone_energy, loss_cone_number = integrate_flux(event,
        time = true,
        time_idxs = start_idx:stop_idx,
        pitch_angle = true,
        pitch_angle_range = loss_cone_limits,
        energy = true
    )
    
    # Exit early if there is no precipitating flux to avoid wasted work. This is a major speedup!
    if loss_cone_number ≤ 10
        return
    end

    # Continue getting data-derived values
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
    ΔE = (event.energy_bins_max .- event.energy_bins_min) ./ 1000 # Units: MeV
    data_counts = [event.n_flux[t,E,α] * (ELFIN_SPIN_PERIOD/16) * ΔE[E] * ELFIN_GEOMETRIC_FACTOR for t in start_idx:stop_idx, E in 1:16, α in 1:16] # Units: number
    data_counts = dropdims(sum(data_counts, dims = 1), dims = 1)

    # Get masks for loss/anti-loss cone and downgoing pitch angle bins
    lc_mask        =      loss_cone_limits[1] .< α_center .< loss_cone_limits[2]
    alc_mask       = anti_loss_cone_limits[1] .< α_center .< anti_loss_cone_limits[2]
    downgoing_mask =      downgoing_limits[1] .< α_center .< downgoing_limits[2]

    # Flip to northern hemisphere for G4EPP if needed
    if α_lc > 90
        data_counts = reverse(data_counts, dims = 2)
        lc_mask = reverse(lc_mask)
        alc_mask = reverse(alc_mask)
        downgoing_mask = reverse(downgoing_mask)
    end














    
    # Simulate backscatter
    simulation_number_backscatter, simulation_energy_backscatter = simulate_elfin_backscatter(event, α_center, data_counts, α_backscatter)
    r = simulation_number_backscatter / anti_loss_cone_number
    @show r
    return
    
    
















    # Cut out bins that would've been discarded in ELFIN data due to low counts
    # Using δq/q ≈ 1/√counts (shot noise), which is what ELFIN does
    sim_relative_error = 1 ./ sqrt.(simulated_backscatter_counts)
    simulated_backscatter_counts[sim_relative_error .> maximum_relative_error] .= 0

    # Flip back to original hemisphere if needed
    if α_lc > 90
        data_counts = reverse(data_counts, dims = 2)
        lc_mask = reverse(lc_mask)
        alc_mask = reverse(alc_mask)
        downgoing_mask = reverse(downgoing_mask)
    end

    # Get simulated backscatter count and energy
    sim_alc_number = sum(simulated_backscatter_counts[:, alc_mask])
    


    # DEBUGGING
    #=
    if false == (((sim_alc_number/loss_cone_number) > 0.5) .&& ((anti_loss_cone_number/loss_cone_number) < 0.4))
        heatmap(log10.(data_counts), clims = (1,3))
        display(plot!())

        backscatter_output_distribution.values[:, .!alc_mask] .= 0

        heatmap(log10.(backscatter_output_distribution.values), clims = (1,3))
        display(plot!())

        sim_rn = sim_alc_number/loss_cone_number
        rn = anti_loss_cone_number/loss_cone_number

        dt = event.time[stop_idx] - event.time[start_idx]

        println()
        @show sim_rn
        @show loss_cone_number
        @show rn
        @show α_lc

        error()
    end
    =#
    


    
    
    


   










    # Write results to file
    event_data = "$(event.time_datetime[start_idx]),$(event.time_datetime[stop_idx]+Microsecond(ELFIN_SPIN_PERIOD*1e6)),$(event.satellite),$(event.L[start_idx]),$(event.L[stop_idx]),$(event.MLT[start_idx]),$(event.MLT[stop_idx])"
    lc_data = "$(loss_cone_energy),$(loss_cone_number)"
    trap_data = "$(trapped_energy),$(trapped_number)"
    alc_data = "$(anti_loss_cone_energy),$(anti_loss_cone_number)"
    error_data = "$(loss_cone_relative_error),$(trapped_relative_error),$(anti_loss_cone_relative_error)"
    sim_data = "$(sim_alc_energy),$(sim_alc_number)"
    n_idxs = "$(sum(lc_mask)),$(sum(alc_mask))"
    
    lock(results_file_threadlock) do
        file = open(results_path, "a")
        write(file, "$(event_data),$(lc_data),$(trap_data),$(alc_data),$(error_data),$(sim_data),$(n_idxs)\n")
        close(file)
    end

    return
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

function simulate_elfin_backscatter(event::Event, α_center, data_counts, α_backscatter)
    beam_energies, beam_pitch_angles = get_beam_locations()
    beam_energies = sort(unique(beam_energies))
    beam_pitch_angles = sort(unique(beam_pitch_angles))
    


    # Get solid angle span of each beam
    beam_midpoints = beam_pitch_angles[begin:end-1] .+ (diff(beam_pitch_angles) ./ 2)
    beam_edges = [beam_pitch_angles[begin], beam_midpoints..., beam_pitch_angles[end]]
    beam_ΔΩ = 2π .* [cosd(beam_edges[i])  - cosd(beam_edges[i+1]) for i in 1:length(beam_edges)-1]

    # Get backscatter
    number_backscattered = 0
    energy_backscattered = 0
    for E in 1:16
        pad = data_counts[E,:]
        interpolator = LinearInterpolator(α_center, pad, NoBoundaries())
                
        
        for beam_idx in eachindex(beam_pitch_angles)
            #ΔΩ = beam_ΔΩ[beam_idx] # Units: str
            
            number = interpolator.([beam_edges[beam_idx],beam_edges[beam_idx+1]]) # Units: # electrons
            number = clamp.(number, 0, Inf)
            
            Δα = beam_edges[beam_idx+1] - beam_edges[beam_idx]
            beam_weight = 0.5 * Δα * (number[1] + number[2])

            if max(number...) == 0
                continue
            end



            backscatter_data, n_input_particles, _ = find_datafile("processed", "backscatter", "electron", round(event.energy_bins_mean[E]), beam_pitch_angles[beam_idx])
            backscatter_mask = backscatter_data["electron_pitch_angles"] .≥ α_backscatter

            number_backscattered += beam_weight * (sum(backscatter_mask)/n_input_particles)
            energy_backscattered += beam_weight * (sum(backscatter_data["electron_energies"][backscatter_mask])/n_input_particles)
        end
    end
    return number_backscattered, energy_backscattered
end


tick()
write_elfin_lifetime_backscatter_data()
tock()

