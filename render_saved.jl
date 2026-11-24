#!/usr/bin/env julia
# Re-render the exact coordinates from the approved image (no sampling).
pushfirst!(LOAD_PATH,@__DIR__)
using SchlottmannDots
const S=SchlottmannDots
function main(args)
    path=isempty(args) ? joinpath(@__DIR__,"examples","approved-coordinates.csv") : abspath(args[1])
    lines=readlines(path)
    lines[1]=="id,ideal_x_px,ideal_y_px,noisy_x_px,noisy_y_px,dx_px,dy_px" || error("Unexpected CSV header")
    original=S.Point[];noisy=S.Point[]
    for line in lines[2:end]
        a=parse.(Float64,split(line,','))
        push!(original,(a[2],a[3]));push!(noisy,(a[4],a[5]))
    end
    _,edges=S.nearest_and_edges(original;radius=120.,edge_length=100.)
    # The bundled approved image is a 9667 px square.
    width=9667;height=9667
    output=joinpath(@__DIR__,"generated","approved-replay");mkpath(output)
    for overlay in (false,true)
        stem=overlay ? "tiling-overlay" : "dots"
        image=S.render_image(original,noisy,edges,width,height;overlay=overlay)
        S.write_png(joinpath(output,stem*".png"),image)
        S.write_svg(joinpath(output,stem*".svg"),original,noisy,edges,width,height;overlay=overlay)
    end
    println("Replayed $(length(noisy)) saved point coordinates without resampling: $output")
end
main(ARGS)
