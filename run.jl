#!/usr/bin/env julia
# Also works with VS Code's "Run Julia File": activate this folder explicitly.
pushfirst!(LOAD_PATH, @__DIR__)
using SchlottmannDots

if length(ARGS) > 1 || (length(ARGS) == 1 && ARGS[1] in ("-h", "--help"))
    println("Usage: julia --project=. run.jl [path/to/config.toml]")
else
    config = isempty(ARGS) ? joinpath(@__DIR__, "config.toml") : abspath(ARGS[1])
    SchlottmannDots.run(config)
end
