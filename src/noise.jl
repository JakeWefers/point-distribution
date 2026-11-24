# Explicit SplitMix64 generator makes seeded output independent of Julia's RNG
# implementation/version. This is reproducible pseudorandom sampling, not a
# cryptographic generator. UInt64 overflow is intentionally modular arithmetic.
function splitmix64(state::UInt64)
    state+=0x9e3779b97f4a7c15
    z=state
    z=(z⊻(z>>30))*0xbf58476d1ce4e5b9
    z=(z⊻(z>>27))*0x94d049bb133111eb
    state,z⊻(z>>31)
end
uniform01(x::UInt64)=Float64(x>>11)*0x1.0p-53

"""One disk-uniform displacement per point. No rejection, optimization or retries.
noise_percent is the maximum radial displacement as a percentage of spacing.
"""
function perturb(points::Vector{Point}; spacing=100.0,noise_percent=10.0,seed=20261006)
    isfinite(noise_percent) && noise_percent>=0 || throw(ArgumentError("noise_percent must be finite and nonnegative"))
    isfinite(spacing) && spacing>0 || throw(ArgumentError("spacing must be finite and positive"))
    seed>=0 || throw(ArgumentError("seed must be nonnegative"))
    state=UInt64(seed);R=spacing*noise_percent/100
    deltas=Point[];noisy=Point[]
    for p in points
        state,u=splitmix64(state);state,v=splitmix64(state)
        radius=R*sqrt(uniform01(u));theta=2pi*uniform01(v)
        d=(radius*cos(theta),radius*sin(theta))
        push!(deltas,d);push!(noisy,add(p,d))
    end
    noisy,deltas
end
