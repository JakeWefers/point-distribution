#!/usr/bin/env julia
# Run with: julia example.jl
# Cylinder dimensions are meters; coordinates, spacing and dot radius are mm.

include(joinpath(@__DIR__, "cylinder", "1_tiling_coordinates.jl"))
include(joinpath(@__DIR__, "cylinder", "2_add_noise.jl"))
include(joinpath(@__DIR__, "cylinder", "3_output_dots.jl"))

function main()
    rectangle_m = (0.20, 0.12) # Circumference and height.
    density_per_in2 = 192.72
    noise_percent = 7.0
    tiling_seed, noise_seed = 123, 456
    dpi = 300

    points, spacing, info = tiling_coordinates(rectangle_m;
        density_per_in2=density_per_in2, seed=tiling_seed)
    noisy_points, error_bound = add_noise(points, rectangle_m;
        spacing_mm=spacing, noise_percent=noise_percent, seed=noise_seed)
    radius_mm = 0.15 * spacing

    image = output_dots(noisy_points, error_bound, rectangle_m, radius_mm;
        spacing_mm=spacing, dpi=dpi, filename=joinpath(@__DIR__, "example.png"))

    println("Generated ", size(points, 1), " points: ", image)
    println("Mean density error: ", info.density_error_percent, "%")
    println("Maximum pair-distance change bound: ", 100 * error_bound, "% of spacing")
    return image
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
