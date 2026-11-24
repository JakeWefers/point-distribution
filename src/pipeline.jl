function require_number(c,key,lower;strict=false)
    x=c[key]
    x isa Real && isfinite(x) && (strict ? x>lower : x>=lower) || throw(ArgumentError("Invalid $key"))
    x
end
function run(configpath::AbstractString)
    configpath=abspath(configpath);c=TOML.parsefile(configpath)
    n=c["points"]; n isa Integer && n>=2 || throw(ArgumentError("points must be an integer >= 2"))
    spacing=Float64(require_number(c,"spacing_px",0;strict=true))
    noise=Float64(require_number(c,"noise_percent",0))
    aspect=Float64(require_number(c,"aspect_ratio",0;strict=true))
    dot=Float64(require_number(c,"dot_radius_px",0;strict=true))
    edge=Float64(require_number(c,"edge_width_px",0;strict=true))
    seed=c["seed"];seed isa Integer && seed>=0 || throw(ArgumentError("seed must be a nonnegative integer"))
    gray=c["gray"];gray isa Integer && 0<=gray<=255 || throw(ArgumentError("gray must be an integer in 0:255"))
    center=Tuple(Float64.(c["crop_center"]))::Point
    all(isfinite,center) || throw(ArgumentError("crop_center must be finite"))
    angle=Float64(c["crop_rotation_radians"]);isfinite(angle) || throw(ArgumentError("crop angle must be finite"))
    generations=c["generations"]
    generations isa Integer && 0<=generations<=6 || throw(ArgumentError("generations must be an integer in 0:6"))
    rules=load_rules()
    if generations==0
        # Conservative estimate for a rectangle inside an inflated square seed.
        need=1.6sqrt(n*max(aspect,1/aspect)/(0.5+1/sqrt(3)))+2hypot(center...)
        generations=max(1,ceil(Int,log(need)/log(rules.inflation)))
    end
    println("Generating $n points, $noise% maximum displacement, seed $seed")
    patch=generate_patch(generations;rules=rules)
    allpoints=vertices(patch,rules)
    crop=crop_points(allpoints,n;aspect_ratio=aspect,center=center,angle=angle)
    verification=verify_crop(patch,rules,crop)
    margin=ceil(noise/100*spacing+dot+.5)
    ideal=[((p[1]+crop.halfwidth)*spacing+margin,(p[2]+crop.halfheight)*spacing+margin) for p in crop.points]
    width=ceil(Int,2crop.halfwidth*spacing+2margin);height=ceil(Int,2crop.halfheight*spacing+2margin)
    nn0,edges=nearest_and_edges(ideal;radius=1.2spacing,edge_length=spacing)
    all(abs(d-spacing)<=1e-6 for d in nn0) || error("Crop contains boundary vertices without a unit-distance neighbor; adjust crop_center slightly")
    noisy,deltas=perturb(ideal;spacing=spacing,noise_percent=noise,seed=seed)
    radius=spacing*noise/100
    all(dot<=p[1]<=width-dot && dot<=p[2]<=height-dot for p in noisy) || error("Insufficient image padding")
    nn,_=nearest_and_edges(noisy;radius=spacing+2radius+1e-6)
    shifts=[hypot(d...) for d in deltas]
    maximum(shifts)<=radius+1e-10 || error("Displacement exceeds requested disk radius")
    minimum(nn)>=max(0,spacing-2radius)-1e-6 || error("Separation bound failed")
    report=Dict{String,Any}(
        "construction"=>"Actual five-prototile Schlottmann substitution; no uniform subdivision",
        "dot_count"=>n,"generations"=>generations,"generated_tiles"=>length(patch),"generated_vertices"=>length(allpoints),
        "image_width_px"=>width,"image_height_px"=>height,"ideal_spacing_px"=>spacing,
        "noiseless_nn_min_px"=>minimum(nn0),"noiseless_nn_max_px"=>maximum(nn0),
        "noise_percent"=>noise,"noise_radius_px"=>radius,"noise_distribution"=>"Uniform area in a disk",
        "noise_generator"=>"SplitMix64; two uniforms per point; r=R*sqrt(u), theta=2*pi*v",
        "seed"=>seed,"samples_per_point"=>1,"resampling"=>false,
        "actual_max_shift_px"=>maximum(shifts),"rms_shift_px"=>sqrt(sum(abs2,shifts)/n),
        "noisy_nn_min_px"=>minimum(nn),"noisy_nn_max_px"=>maximum(nn),
        "max_relative_nn_error_percent"=>100maximum(abs.(nn.-spacing))/spacing,
        "verification"=>verification)
    output=abspath(joinpath(dirname(configpath),c["output_dir"]));mkpath(output)
    open(joinpath(output,"statistics.toml"),"w") do io;TOML.print(io,report);end
    open(joinpath(output,"config-used.toml"),"w") do io;TOML.print(io,c);end
    open(joinpath(output,"coordinates.csv"),"w") do io
        println(io,"id,ideal_x_px,ideal_y_px,noisy_x_px,noisy_y_px,dx_px,dy_px")
        for i in eachindex(ideal)
            @printf(io,"%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",i,ideal[i]...,noisy[i]...,deltas[i]...)
        end
    end
    open(joinpath(output,"edges.csv"),"w") do io
        println(io,"point_id_1,point_id_2")
        for (i,j) in edges;println(io,"$i,$j");end
    end
    for overlay in (false,true)
        stem=overlay ? "tiling-overlay" : "dots"
        if c["export_svg"]
            write_svg(joinpath(output,stem*".svg"),ideal,noisy,edges,width,height;overlay=overlay,dot_radius=dot,edge_width=edge,gray=gray)
        end
        if c["export_png"]
            println("Rendering $stem ($width × $height)")
            image=render_image(ideal,noisy,edges,width,height;overlay=overlay,dot_radius=dot,edge_width=edge,gray=gray)
            write_png(joinpath(output,stem*".png"),image)
            image=nothing;GC.gc()
            s=min(1.,c["preview_max_px"]/max(width,height))
            pw=max(1,round(Int,width*s));ph=max(1,round(Int,height*s))
            preview=render_image(ideal,noisy,edges,pw,ph;overlay=overlay,view=(0.,0.,Float64(width),Float64(height)),dot_radius=dot,edge_width=edge,gray=gray)
            write_png(joinpath(output,stem*"-preview.png"),preview)
            side=min(Float64(c["detail_side_px"]),width,height);ds=Int(c["detail_output_px"])
            detail=render_image(ideal,noisy,edges,ds,ds;overlay=overlay,view=((width-side)/2,(height-side)/2,side,side),dot_radius=dot,edge_width=edge,gray=gray)
            write_png(joinpath(output,stem*"-detail.png"),detail)
        end
    end
    @printf("Verified %d points. Noisy nearest-neighbor range: %.5f–%.5f px.\n",n,minimum(nn),maximum(nn))
    println("Saved to: $output")
    report
end
