using Random

function add_noise(coordinates_mm::AbstractMatrix, rectangle_size_m; spacing_mm,
                   noise_percent=nothing, max_shift_mm=nothing, seed=20261006)
    n = size(coordinates_mm,1)
    size(coordinates_mm,2)==2 && n>=2 ||
        throw(ArgumentError("coordinates_mm must be an N × 2 matrix with N >= 2."))
    length(rectangle_size_m)==2 || throw(ArgumentError("Use (circumference_m, height_m)."))
    width, height = 1000 .* Float64.(rectangle_size_m)
    isfinite(width) && isfinite(height) && width>0 && height>0 ||
        throw(ArgumentError("Rectangle dimensions must be positive and finite."))
    spacing = Float64(spacing_mm)
    isfinite(spacing) && spacing>0 || throw(ArgumentError("spacing_mm must be positive and finite."))
    (noise_percent===nothing) != (max_shift_mm===nothing) ||
        throw(ArgumentError("Specify exactly one of noise_percent or max_shift_mm."))
    seed isa Integer && seed>=0 || throw(ArgumentError("seed must be a nonnegative integer."))
    shift = max_shift_mm===nothing ? Float64(noise_percent)*spacing/100 : Float64(max_shift_mm)
    isfinite(shift) && shift>=0 || throw(ArgumentError("Noise must be finite and nonnegative."))
    isfinite(width+shift) && isfinite(height+shift) || throw(ArgumentError("Noise displacement is too large."))
    relative_error_bound = 2*(shift/spacing)
    isfinite(relative_error_bound) || throw(ArgumentError("Noise bound is too large."))

    noisy_coordinates_mm = Matrix{Float64}(undef,n,2)
    rng = MersenneTwister(seed)
    @inbounds for i in 1:n
        x, y = coordinates_mm[i,1], coordinates_mm[i,2]
        isfinite(x) && isfinite(y) && 0<=x<=width && 0<=y<=height ||
            throw(ArgumentError("Original coordinates must be finite and inside the rectangle."))
        r = shift*sqrt(rand(rng))
        theta = 2pi*rand(rng)
        noisy_coordinates_mm[i,1] = mod(x+r*cos(theta),width)
        noisy_coordinates_mm[i,2] = y+r*sin(theta)
    end
    return noisy_coordinates_mm, relative_error_bound
end
