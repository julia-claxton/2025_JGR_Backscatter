const TOP_LEVEL = dirname(@__DIR__)

include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Events.jl")
include("/Users/luna/Documents/Work/Research/Julia_ELFIN_Tools/Visualization.jl")
include("/Users/luna/Documents/Work/Research/Backscatter_Analysis/code/General_Functions.jl")
include("/Users/luna/Documents/Work/Research/G4EPP_2.0/Frontend_Functions.jl")

function main()
    path = "/Users/luna/Research/Backscatter_Analysis/data/ELFIN_backscatter_and_simulation.csv"
    data = readdlm(path, ',', skipstart = 1)

    start = DateTime.(data[:,1])
    stop = DateTime.(data[:,2])
    sat = data[:,3]

    for i in eachindex(start)
        println(i/length(start))
        event = create_event(start[i], stop[i], sat[i])

        if i % 30 == 0; GC.gc(); end
    end
end


main()