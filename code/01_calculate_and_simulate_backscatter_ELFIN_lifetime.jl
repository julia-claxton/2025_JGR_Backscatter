const TOP_LEVEL = dirname(@__DIR__)
const G4EPP_TOP_LEVEL = "$(TOP_LEVEL)/code/G4EPP_2.0/"
include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Events.jl")
include("$(TOP_LEVEL)/code/Julia_ELFIN_Tools/Visualization.jl")
include("$(TOP_LEVEL)/code/General_Functions.jl")
using Base.Threads
using Statistics, LinearAlgebra
using BenchmarkTools, Profile, TickTock


# Run with max threads
# julia --threads 10 /Users/luna/Research/Backscatter_Analysis/code/01_calculate_and_simulate_backscatter_ELFIN_lifetime.jl

write_threadlock  = ReentrantLock()
stdout_threadlock = ReentrantLock()

function main()
    # Parameters
    result_filename = "ELFIN_backscatter_and_simulation.csv"

    # Create results file
    result_path = "$(dirname(@__DIR__))/data/$(result_filename)"
    file = open(result_path, "w")
    write(file, "start_time,stop_time,sat,start_L,stop_L,start_MLT,stop_MLT,lc_energy,lc_number,trap_energy,trap_number,alc_energy,alc_number,lc_relative_error,trap_relative_error,alc_relative_error,simulation_alc_energy,simulation_alc_number,lc_n_sectors,alc_n_sectors\n")
    close(file)

    # Iterate through ELFIN lifetime
    dates, sats = all_elfin_science_dates_and_satellite_ids()
    n_dates = length(dates)
    days_analyzed = zeros(nthreads())
    iteration_assignments = [id:nthreads():n_dates for id in 1:nthreads()]

    Threads.@threads for id in 1:nthreads()
        # Manually designate iterations assigned to this thread to balance workload between threads
        # Use interleaving assignments rather than chunking
        this_thread_idxs =  iteration_assignments[id] # Iterations this thread is responsible for

        # Perform analysis for all days this thread is responsible for
        for i in eachindex(this_thread_idxs)
            # Do analysis
            day_to_analyze = this_thread_idxs[i]
            write_elfin_backscatter_for_day(dates, sats, day_to_analyze, result_path)

            lock(write_threadlock) do
                days_analyzed[id] += 1
                if sum(days_analyzed) % 10 == 0; GC.gc(); end
                unsafe_print_progress_screen(days_analyzed, iteration_assignments)
            end
        end
    end
end

function unsafe_print_progress_screen(days_analyzed, iteration_assignments)
    run(`clear`)
    println(
        """
        ELFIN Lifetime Backscatter Analysis
        --------------------------------------------------------
        Computing backscatter and simulation data...
        """
    )
    for i in 1:nthreads()
        print("Thread $(i)\t")
        print_progress_bar(days_analyzed[i]/length(iteration_assignments[i]), overwrite = false, bar_length = 30)
        println()
    end
end

function write_elfin_backscatter_for_day(dates, sats, i, result_path)
    slice_length_idx = 3
    maximum_relative_error = .5

    fullday_event = create_event(dates[i], sats[i])
    if fullday_event == nothing; return; end
    if fullday_event.data_reliable == false; return; end

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
    sim_backscatter_number, sim_backscatter_energy = simulate_elfin_backscatter(event, α_center, data_counts, alc_mask)

    # Cut out bins that would've been discarded in ELFIN data due to low counts
    # Using δq/q ≈ 1/√counts (shot noise), which is what ELFIN does
    sim_relative_error = 1 ./ sqrt.(sim_backscatter_number)
    sim_backscatter_number[sim_relative_error .> maximum_relative_error] .= 0
    sim_backscatter_energy[sim_relative_error .> maximum_relative_error] .= 0

    # Get totals
    sim_alc_number = sum(sim_backscatter_number)
    sim_alc_energy = sum(sim_backscatter_energy)
    
    # Write results to file
    event_data = "$(event.time_datetime[start_idx]),$(event.time_datetime[stop_idx]+Microsecond(ELFIN_SPIN_PERIOD*1e6)),$(event.satellite),$(event.L[start_idx]),$(event.L[stop_idx]),$(event.MLT[start_idx]),$(event.MLT[stop_idx])"
    lc_data = "$(loss_cone_energy),$(loss_cone_number)"
    trap_data = "$(trapped_energy),$(trapped_number)"
    alc_data = "$(anti_loss_cone_energy),$(anti_loss_cone_number)"
    error_data = "$(loss_cone_relative_error),$(trapped_relative_error),$(anti_loss_cone_relative_error)"
    sim_data = "$(sim_alc_energy),$(sim_alc_number)"
    n_idxs = "$(sum(lc_mask)),$(sum(alc_mask))"
    
    lock(write_threadlock) do
        file = open(results_path, "a")
        write(file, "$(event_data),$(lc_data),$(trap_data),$(alc_data),$(error_data),$(sim_data),$(n_idxs)\n")
        close(file)
    end
    return
end

function simulate_elfin_backscatter(event::Event, α_center, data_counts, alc_mask)
    # Get locations of precomputed beams
    beam_energies, beam_pitch_angles = get_beam_locations()
    beam_energies     = sort(unique(beam_energies))
    beam_pitch_angles = sort(unique(beam_pitch_angles))

    # Iterate over ELFIN's downgoing bins
    sim_alc_counts = zeros(size(data_counts[:,alc_mask]))
    for E_idx in 1:16
        for α_idx in 1:2:16
            # Get energy and pitch angle of this bin
            E = event.energy_bins_mean[E_idx]
            α = α_center[α_idx]

            # Guard block
            if data_counts[E_idx, α_idx] == 0; continue; end
            if !(0 .≤ α .≤ 90); continue; end

            # Get nearest beam
            _, nearest_beam_energy_idx = findmin(abs.(beam_energies .- E))
            nearest_beam_energy = beam_energies[nearest_beam_energy_idx]
            
            _, nearest_beam_pitch_angle_idx = findmin(abs.(beam_pitch_angles .- α))
            nearest_beam_pitch_angle = beam_pitch_angles[nearest_beam_pitch_angle_idx]

            # Get counts
            beam_weight = data_counts[E_idx, α_idx] * azimuth_scaling_factor(α)
            sim_alc_counts .+= get_beam_backscatter(event, beam_weight, nearest_beam_energy, nearest_beam_pitch_angle, α_center[alc_mask])
        end
    end
    sim_alc_energy = [sim_alc_counts[E,α] .* event.energy_bins_mean[E] for E in 1:16, α in 1:size(sim_alc_counts)[2]]
    return sim_alc_counts, sim_alc_energy
end

function azimuth_scaling_factor(α)
    Ω_elfin_epd = 2π * (1 - cosd(ELFIN_EPD_FOV/2))
    α_min = clamp(α - ELFIN_EPD_FOV/2, 0, 180)
    α_max = clamp(α + ELFIN_EPD_FOV/2, 0, 180)
    Ω_ribbon = 2π * (cosd.(α_min) .- cosd.(α_max))
    return Ω_ribbon / Ω_elfin_epd
end

function get_beam_backscatter(event::Event, beam_weight, nearest_beam_energy, nearest_beam_pitch_angle, α_detector)
    backscatter_number = zeros(length(event.energy_bins_mean), length(α_detector))
    for E_idx in eachindex(event.energy_bins_mean)
        for α_idx in eachindex(α_detector)
            backscatter_number[E_idx, α_idx] = beam_weight * detect_backscatter(nearest_beam_energy, nearest_beam_pitch_angle, α_detector[α_idx],
                detector_Emin = event.energy_bins_min[E_idx],
                detector_Emax = event.energy_bins_max[E_idx]
            )
        end
    end
    return backscatter_number
end

function detect_backscatter(beam_energy, beam_pitch_angle, detector_α; detector_Emin = -Inf, detector_Emax = Inf)
    # Get backscatter data
    global backscatter_momentum_data
    magnetic_momenta = backscatter_momentum_data["$(beam_energy)keV_$(beam_pitch_angle)deg_momenta"]
    energies = backscatter_momentum_data["$(beam_energy)keV_$(beam_pitch_angle)deg_energies"]

    # Construct virtual EPD
    EPD_FOV = 22.5 # Detector field of view, deg
    detector_boresight = get_detector_look_direction(detector_α)
    detection_mask = (
        (acosd.(clamp.(dot.([detector_boresight], magnetic_momenta), -1, 1)) .< (EPD_FOV/2))
        .&& (detector_Emin .< energies .< detector_Emax)
    )

    # Return
    number_detected = sum(detection_mask) / 1e5 # Divide by number of input particles for that beam
    return number_detected
end

function find_datafile(data_type, file_prefix, input_particle, beam_energy, beam_pitch_angle; quiet = false)
    # Get the file corresponding to this prefix, energy, and pitch angle. If there are multiple files that have different numbers of input
    # particles, this will select the run with the larger number of input particles. There shouldn't be multiple files like that,
    # but we are being safe just in case.
    if data_type == "raw"
        file_extension = "csv"
    elseif data_type == "processed"
        file_extension = "npz"
    else
        error("Please provide data_type of either \"raw\" or \"processed\". Got $(data_type)")
    end

    available_files = glob("$(file_prefix)_$(input_particle)_input_$(beam_energy)keV_$(beam_pitch_angle)deg_*.$(file_extension)", "$(G4EPP_TOP_LEVEL)/data/$(data_type)")
    
    # If we don't find the file...
    if length(available_files) == 0
        if quiet == false; @warn "$(file_prefix)_$(input_particle)_input_$(beam_energy)keV_$(beam_pitch_angle)deg_*.$(file_extension) not found!"; end
        return nothing, nothing, nothing
    end

    # Find number of input particles each file at this energy and pitch angle has
    matches_n_particles = match.(Regex("deg_(.*?)particles.$(file_extension)"), available_files) # Matches for the pattern containing number of particles
    number_of_particles = [n_particle_match.captures[1] for n_particle_match in matches_n_particles] # Number of particles contained in each match
    number_of_particles = parse.(Int64, number_of_particles) # Cast to integer

    # Select file with largest number of particles
    n_particles = max(number_of_particles...)

    # Construct filename of data file and load it in
    filename = "$(file_prefix)_$(input_particle)_input_$(beam_energy)keV_$(beam_pitch_angle)deg_$(n_particles)particles"
    filepath = "$(G4EPP_TOP_LEVEL)/data/$(data_type)/$(filename).$(file_extension)"

    if data_type == "raw"
        return readdlm(filepath, ',', skipstart = 1), n_particles, filename
    elseif data_type == "processed"
        return npzread(filepath), n_particles, filename
    end
end

function get_beam_locations()
    raw_files = glob("backscatter_electron_input_*", "/Users/luna/Research/G4EPP_2.0/data/processed")
    matches = match.(Regex("input_(.*?)keV_(.*?)deg"), raw_files)
    beam_energies = [pattern_match.captures[1] for pattern_match in matches]
    beam_pitch_angles = [pattern_match.captures[2] for pattern_match in matches]

    # Cast to float and return
    beam_energies = parse.(Float64, beam_energies)
    beam_pitch_angles = parse.(Float64, beam_pitch_angles)

    return beam_energies, beam_pitch_angles
end

function get_backscattered_magnetic_momenta_and_energies(input_energy, input_pitch_angle)
    # Read data
    data = readdlm("$(TOP_LEVEL)/code/G4EPP_2.0/data/raw/backscatter_electron_input_$(float(input_energy))keV_$(float(input_pitch_angle))deg_100000particles.csv", ',', skipstart = 1)
    electron_mask = data[:,1] .== "e-"

    energy = data[electron_mask, 2]
    position = eachrow(data[electron_mask, 7:9])
    momentum = eachrow(data[electron_mask, 4:6])
    B = get_B.(position)
    unit_B = B ./ norm.(B)

    # Define the magnetic frame
    # Define it such that the X axis is always parallel to the ground, the Z axis points along B, and Y completes the set.
    Z_world = [0, 0, 1]
    Z_magnetic = copy(unit_B)

    X_magnetic = Z_magnetic .× [Z_world]
    X_magnetic ./= norm.(X_magnetic)

    Y_magnetic = Z_magnetic .× X_magnetic

    # Convert momentum to the magnetic frame
    momentum_magnetic = [[
        momentum[i] ⋅ X_magnetic[i],
        momentum[i] ⋅ Y_magnetic[i],
        momentum[i] ⋅ Z_magnetic[i]
        ] for i in eachindex(momentum)
    ]
    return momentum_magnetic, energy
end

function get_B(r_origin_to_particle)
    # Constants
    μ0 = 1.257e-6
    MLAT = 65.77 # degrees
    m = [0, -6.4e22 * cosd(MLAT), -6.4e22 * sind(MLAT)] # Earth magnetic moment vector

    # Get position vector to Earth center
    Re = 6371e3 # Earth radius, m
    r = [0, 0, Re] + r_origin_to_particle # Position vector of Earth to particle. Units: m

    # Calculate dipole field
    return (μ0/(4π)) * ( ((3 * dot(m,r) .* r) / (norm(r)^5)) - (m/(norm(r)^3)) )
end

function get_detector_look_direction(detector_α)
    # Get look direction in the magnetic frame for a given pitch angle
    detector_ϕ = 0 # Virtual detector azimuth angle, deg.

    return [
        sind(detector_α)*cosd(detector_ϕ),
        sind(detector_α)*sind(detector_ϕ),
        cosd(detector_α)
    ]
end

function load_all_backscatter_into_memory()
    # Get locations of precomputed beams
    beam_energies, beam_pitch_angles = get_beam_locations()
    beam_energies     = sort(unique(beam_energies))
    beam_pitch_angles = sort(unique(beam_pitch_angles))

    # Iterate and store
    run(`clear`)
    println(
        """
        ELFIN Lifetime Backscatter Analysis
        --------------------------------------------------------
        Loading backscatter momenta into memory...
        """
    )
    beams_finished = 0
    total_beams = length(beam_energies) * length(beam_pitch_angles)
    backscatter_data = Dict()
    Threads.@threads for E in round.(beam_energies)
        for α in round.(beam_pitch_angles)
            magnetic_momenta, energies = get_backscattered_magnetic_momenta_and_energies(E, α)
            
            # Write to dict (safely!)
            lock(write_threadlock) do
                backscatter_data["$(E)keV_$(α)deg_momenta"] = magnetic_momenta
                backscatter_data["$(E)keV_$(α)deg_energies"] = energies
            end

            # Print progress to terminal
            lock(stdout_threadlock) do
                beams_finished += 1
                print_progress_bar(beams_finished/total_beams, bar_length = 30)
            end
        end
    end
    return backscatter_data
end

#################
# Figure stuff
#################
function save_beam_weighting_procedure()
    print("Saving beam weighting data... ")
    # Get event
    maximum_relative_error = 0.5
    event = create_event(DateTime("2021-03-06T07:05:05"), DateTime("2021-03-06T07:05:15"), "B", maximum_relative_error = maximum_relative_error)
    _, n_fluence = integrate_flux(event, time = true)
    avg_nflux = n_fluence / event.duration

    α_lc = mean(event.loss_cone_angles)
    α_alc = mean(event.anti_loss_cone_angles)
    α_center = dropdims(mean(event.pitch_angles, dims = 1), dims = 1)

    cone_standoff_angle = ELFIN_EPD_FOV / 2 # Minimum distance into the loss/antiloss cone a pitch angle bin center must be in order to be counted 
    
    # Get loss cone, trapped, and anti-loss cone region bounds
    # Northern hemisphere
    loss_cone_limits = (0, α_lc - cone_standoff_angle)
    trapped_limits = (α_lc, α_alc)
    anti_loss_cone_limits = (α_alc + cone_standoff_angle, 180)
    downgoing_limits = (0, 90)

    # Simulate event with G4EPP
    ΔE = (event.energy_bins_max .- event.energy_bins_min) ./ 1000 # Units: MeV
    data_counts = [event.n_flux[t,E,α] * (ELFIN_SPIN_PERIOD/16) * ΔE[E] * ELFIN_GEOMETRIC_FACTOR for t in 1:event.n_datapoints, E in 1:16, α in 1:16] # Units: number
    data_counts = dropdims(sum(data_counts, dims = 1), dims = 1)

    # Get masks for loss/anti-loss cone and downgoing pitch angle bins
    lc_mask        = loss_cone_limits[1] .< α_center .< loss_cone_limits[2]
    alc_mask       = anti_loss_cone_limits[1] .< α_center .< anti_loss_cone_limits[2]
    downgoing_mask = downgoing_limits[1] .< α_center .< downgoing_limits[2]

    # Get beam weights
    beam_energies, beam_pitch_angles, beam_weights = get_beam_weights(event::Event, α_center, data_counts, alc_mask)

    # Simulate backscatter
    sim_counts, _ = simulate_elfin_backscatter(event, α_center, data_counts, alc_mask)

    # Cut out bins that would've been discarded in ELFIN data due to low counts
    # Using δq/q ≈ 1/√counts (shot noise), which is what ELFIN does
    sim_relative_error = 1 ./ sqrt.(sim_counts)
    sim_counts[sim_relative_error .> maximum_relative_error] .= 0
    
    npzwrite("$(TOP_LEVEL)/data/figure_data/beam_weighting.npz",
        avg_pitch_angles = α_center,
        avg_energies = event.energy_bins_mean,
        avg_nflux = avg_nflux,
        data_counts = data_counts,
        beam_energies = beam_energies,
        beam_pitch_angles = beam_pitch_angles,
        beam_weights = beam_weights,
        sim_counts = sim_counts
    )
    println("Done")
end

function get_beam_weights(event::Event, α_center, data_counts, alc_mask)
    # Get locations of precomputed beams
    beam_energies, beam_pitch_angles = get_beam_locations()
    beam_energies     = sort(unique(beam_energies))
    beam_pitch_angles = sort(unique(beam_pitch_angles))

    # Iterate over ELFIN's downgoing bins
    weighted_energies = []
    weighted_pitch_angles = []
    beam_weights = []
    for E_idx in 1:16
        for α_idx in 1:2:16
            # Get energy and pitch angle of this bin
            E = event.energy_bins_mean[E_idx]
            α = α_center[α_idx]

            # Guard block
            if data_counts[E_idx, α_idx] == 0; continue; end
            if !(0 .≤ α .≤ 90); continue; end

            # Get nearest beam
            _, nearest_beam_energy_idx = findmin(abs.(beam_energies .- E))
            nearest_beam_energy = beam_energies[nearest_beam_energy_idx]
            
            _, nearest_beam_pitch_angle_idx = findmin(abs.(beam_pitch_angles .- α))
            nearest_beam_pitch_angle = beam_pitch_angles[nearest_beam_pitch_angle_idx]

            # Store
            push!(weighted_energies, nearest_beam_energy)
            push!(weighted_pitch_angles, nearest_beam_pitch_angle)
            push!(beam_weights, data_counts[E_idx, α_idx] * azimuth_scaling_factor(α))
        end
    end
    return float.(weighted_energies), float.(weighted_pitch_angles), float.(beam_weights)
end

#=
tick()
# Preload backscatter data
const global backscatter_momentum_data = load_all_backscatter_into_memory() # Globals are hacky but I want this to be accessible in any scope without passing it as an argument to literally everything
# Do analysis
main()
tock()
=#

# Figure data saving
const global backscatter_momentum_data = load_all_backscatter_into_memory() # Globals are hacky but I want this to be accessible in any scope without passing it as an argument to literally everything
save_beam_weighting_procedure()