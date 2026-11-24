#!/usr/bin/env julia
# A small command-line example. Dimensions are meters; coordinates are mm.
using DelimitedFiles, TOML
include(joinpath(@__DIR__, "1_tiling_coordinates.jl"))
include(joinpath(@__DIR__, "2_add_noise.jl"))
include(joinpath(@__DIR__, "3_output_dots.jl"))

function main()
    rectangle_m = (0.2, 0.12)   # Circumference and height.
    density_per_in2 = 192.72
    noise_percent = 7.0
    tiling_seed, noise_seed = 123, 456
    dpi = 300
    output = joinpath(@__DIR__, "..", "generated", "cylinder")
    mkpath(output)

    points, spacing, info = tiling_coordinates(rectangle_m;
        density_per_in2=density_per_in2, seed=tiling_seed)
    noisy, error_bound = add_noise(points, rectangle_m;
        spacing_mm=spacing, noise_percent=noise_percent, seed=noise_seed)
    # Choose a safely nonoverlapping circle radius relative to the tile edge.
    radius_mm = 0.15spacing
    path = output_dots(noisy, error_bound, rectangle_m, radius_mm;
        spacing_mm=spacing, dpi=dpi, filename=joinpath(output, "pattern.png"))

    open(joinpath(output, "coordinates.csv"), "w") do io
        println(io, "ideal_x_mm,ideal_y_mm,noisy_x_mm,noisy_y_mm")
        writedlm(io, hcat(points, noisy), ',')
    end
    # TOML represents the generator's tuple metadata as arrays.
    report = Dict{String,Any}(string(k)=>(v isa Tuple ? collect(v) : v) for (k,v) in pairs(info))
    merge!(report, Dict("circumference_m"=>rectangle_m[1], "height_m"=>rectangle_m[2],
        "spacing_mm"=>spacing, "noise_percent"=>noise_percent,
        "noise_distance_change_bound_fraction"=>error_bound,
        "tiling_seed"=>tiling_seed, "noise_seed"=>noise_seed,
        "dot_radius_mm"=>radius_mm, "dpi"=>dpi))
    open(joinpath(output, "statistics.toml"), "w") do io
        TOML.print(io, report)
    end
    println("Generated ", size(points,1), " cylinder points: ", path)
    println("Mean-density error: ", info.density_error_percent, "%")
    output
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
