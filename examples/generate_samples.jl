#!/usr/bin/env julia
pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
using SchlottmannDots, TOML
const S = SchlottmannDots

root = normpath(joinpath(@__DIR__, ".."))
configpath = joinpath(root, "config.toml")
config = TOML.parsefile(configpath)
S.run(configpath)
generated = abspath(joinpath(root, config["output_dir"]))
results = joinpath(@__DIR__, "results")
mkpath(results)

for (source, target) in (("dots-preview.png", "planar-overview.png"),
                         ("tiling-overlay-detail.png", "planar-overlay-detail.png"),
                         ("statistics.toml", "planar-statistics.toml"),
                         ("config-used.toml", "planar-config.toml"))
    cp(joinpath(generated,source), joinpath(results,target); force=true)
end

# Reuse the generated original coordinates for a directly comparable detail.
rows = readlines(joinpath(generated, "coordinates.csv"))[2:end]
original = S.Point[]
for row in rows
    values = parse.(Float64, split(row, ','))
    push!(original, (values[2], values[3]))
end
stats = TOML.parsefile(joinpath(generated,"statistics.toml"))
w, h = stats["image_width_px"], stats["image_height_px"]
side = min(Float64(config["detail_side_px"]), w, h)
ds = Int(config["detail_output_px"])
image = S.render_image(original, original, Tuple{Int,Int}[], ds, ds;
    view=((w-side)/2, (h-side)/2, side, side), dot_radius=config["dot_radius_px"])
S.write_png(joinpath(results,"planar-original-detail.png"), image)

module CylinderExample
include(joinpath(@__DIR__, "..", "cylinder", "run.jl"))
end
cylinder_output = CylinderExample.main()
cp(joinpath(cylinder_output,"pattern.png"), joinpath(results,"cylinder-pattern.png"); force=true)
cp(joinpath(cylinder_output,"statistics.toml"), joinpath(results,"cylinder-statistics.toml"); force=true)
println("Selected samples saved to: ", results)
