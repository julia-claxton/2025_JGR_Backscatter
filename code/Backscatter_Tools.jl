# Library includes
using Statistics, LinearAlgebra         # Basic math functions
using BenchmarkTools, Profile, TickTock # Debugging tools
using Plots, Plots.PlotMeasures         # Visualization tools
using DelimitedFiles                    # Read/write .csv files
using NPZ                               # Read/write .npz files
using Glob                              # File searching and pattern matching
using NumericalIntegration              # Numerical integration functions

# Written by Julia Claxton. Contact: julia.claxton@colorado.edu
# Released under MIT License


function _find_datafile(data_type, file_prefix, input_particle, beam_energy, beam_pitch_angle; quiet = false)
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







function _counts_to_beam_weights(distribution::Distribution)
    """
    Inputs:
        TODO

    Returns:
        TODO

    Description:
    TODO
    """
    @assert distribution.type == "counts"

    # Find nearest beam location for each nonzero input and assign that bin's particles to that beam
    beam_energies, beam_pitch_angles = _get_beam_locations()
    average_ΔE_logspace = mean(diff(sort(unique(log10.(beam_energies))))) # Wow, that's a lot of operations on one line!
    average_Δα = mean(diff(sort(unique(beam_pitch_angles))))

    # Iterate over every energy-pitch angle bin in input and find its nearest beam, then add weight to it
    beam_weights = zeros(length(beam_energies))
    for bin_indices in CartesianIndices(distribution.values)
        # Note: bin_indices = [energy index, pitch angle index]

        # Get energy and pitch angle of this bin
        bin_energy = distribution.energy_bins_mean[bin_indices[1]]
        bin_pitch_angle = distribution.pitch_angle_bins_mean[bin_indices[2]]

        # Skip bins with no flux
        if distribution.values[bin_indices] == 0; continue; end 

        # Don't simulate upgoing particles or erroneous negative pitch angles
        if (0 < bin_pitch_angle < 90) == false; continue; end

        # Find index of beam with closest energy and pitch angle
        nearest_energy, _ = findmin(abs.(beam_energies .- bin_energy))
        nearest_pitch_angle, _ = findmin(abs.(beam_pitch_angles .- bin_pitch_angle))

        ΔE_logspace = abs.(log10.(bin_energy) .- log10.(beam_energies))
        normalized_energy_distance = ΔE_logspace ./ average_ΔE_logspace # Distance between beam and bin normalized to average beam spacing

        Δα = abs.(bin_pitch_angle .- beam_pitch_angles)
        normalized_pitch_angle_distance = Δα ./ average_Δα # Distance between beam and bin normalized to average beam spacing

        normalized_distances = norm.(eachrow(cat(normalized_energy_distance, normalized_pitch_angle_distance, dims = 2)))
        _, nearest_beam_idx = findmin(normalized_distances)

        # TODO assume grid
        # Only allow a beam to grab bins that are within 3/4 the average bin distance
        if (normalized_energy_distance[nearest_beam_idx] > .75) || (normalized_pitch_angle_distance[nearest_beam_idx] > .75)
            #continue
        end

        # Add the number of electrons measured in this bin to the nearest beam
        beam_weights[nearest_beam_idx] += distribution.values[bin_indices]
    end
    @assert sum(beam_weights) - sum(distribution.values) < 1 "Beam weights total = $(sum(beam_weights)), counts total = $(sum(distribution.values))"
    return beam_weights
end

function _get_beam_locations(; data_type = "processed")
    """
    Inputs:
        None

    Returns:
        beam_energies_keV: Every energy that a beam was run at. Units: keV
        beam_pitch_angles_deg: Every pitch angle that a beam was run at. Units: deg

    Description:
    Gets the energies and pitch angles that beams were run at. Beams were run for every combination of
    energy and pitch angle that is returned - i.e., these are the axis labels for the beam grid in
    energy-pitch angle space.
    """

    raw_files = glob("backscatter_electron_input_*", "$(G4EPP_TOP_LEVEL)/data/$(data_type)") # Arbitrarily choose electron backscatter to avoid double counting different types of data
    matches = match.(Regex("input_(.*?)keV_(.*?)deg"), raw_files)
    beam_energies = [pattern_match.captures[1] for pattern_match in matches]
    beam_pitch_angles = [pattern_match.captures[2] for pattern_match in matches]

    # Cast to float and return
    beam_energies = parse.(Float64, beam_energies)
    beam_pitch_angles = parse.(Float64, beam_pitch_angles)

    return beam_energies, beam_pitch_angles
end
