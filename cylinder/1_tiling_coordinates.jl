"""
    tiling_coordinates((circumference_m, height_m); density_per_in2,
                       refinement_steps=nothing, seed=20261009)
        -> (coordinates_mm, spacing_mm, info)

Generate vertices of a seamless square/equilateral-triangle tiling on a cylinder.
Coordinates are an N × 2 Float64 matrix in MILLIMETERS, with 0 <= x < circumference.
The rectangle dimensions are in METERS. Density input and all reported densities
are in POINTS PER SQUARE INCH. No noise is applied here.

This is a randomized periodic Stampfli approximant, not a strictly aperiodic
Schlottmann tiling. Unrolling a cylinder necessarily repeats after one full turn.
At each refinement, every distinct vertex independently selects one of the two
compatible dodecagon fillings (rotated by 30 degrees), using a reproducible seeded
hash. Periodic copies share the same choice. Repeated local motifs remain.

By default, refinement depth AND coarse triangular-lattice period (p,q) are
chosen from the requested density and circumference. Depths 1 through 8 are
searched; the deepest whose closest compatible MEAN density is within 0.1%
of the target is selected. This is a practical hierarchy-selection rule, not
an optimization or guarantee of minimal visible repetition. If no depth meets
that tolerance, the closest candidate is used with a warning and a flag in info.
An explicit refinement_steps integer overrides the automatic choice.
The tiling is uniformly scaled to the exact circumference. All tile edges have
the returned spacing_mm.
The finite rectangle's point count can differ from its mean-density prediction
because of its top and bottom boundaries. No points are removed to force a count.
info reports both densities and the relative mean-density error in percent.

For circumference 2 m and density 258.064/in², automatic selection uses two steps
and gives mean density 258.115468612/in². Height is freely selectable.
More steps introduce larger hierarchical motifs, not a guarantee of less visible
repetition. Choices vary throughout the full circumference at every step.

Exact integer coordinates in Q(sqrt(3)) prevent duplicate shared vertices.
Only intermediate tiles and the final vertex set are stored. No neighbor search.
A fixed seed gives identical points (including ordering) for identical inputs.

Example:
    points, spacing, info = tiling_coordinates((2.0, 1.0);
        density_per_in2=258.064, seed=123)
"""
function tiling_coordinates(rectangle_size_m; density_per_in2,
                            refinement_steps=nothing, seed=20261009)
    length(rectangle_size_m)==2 || throw(ArgumentError("Use (circumference_m, height_m)."))
    cm, hm = Float64.(rectangle_size_m)
    all(isfinite,(cm,hm)) && cm>0 && hm>0 || throw(ArgumentError("Dimensions must be positive and finite."))
    requested_density = Float64(density_per_in2)
    density = requested_density/6.4516 # Internal geometry uses cm²; 1 in² = 6.4516 cm².
    isfinite(density) && density>0 || throw(ArgumentError("Density must be positive and finite."))
    automatic = refinement_steps===nothing
    automatic || (refinement_steps isa Integer && 1<=refinement_steps<=8) ||
        throw(ArgumentError("refinement_steps must be nothing (automatic), or an integer from 1 to 8."))
    seed isa Integer && 0<=seed<=typemax(UInt64) || throw(ArgumentError("seed must fit a nonnegative UInt64."))
    root3 = sqrt(3.0); inflation = 2+root3
    tolerance = 0.001 # 0.1% relative error in the mean density.
    k = 0; p = 1; q = 0; norm2 = 1; vertices_per_cell = 0
    selected_error = Inf; found_acceptable = false
    triangles, squares = 2, 0
    for depth in 1:(automatic ? 8 : Int(refinement_steps))
        triangles, squares = 7triangles+16squares, 3triangles+7squares
        automatic || depth==refinement_steps || continue
        cell_vertices = triangles÷2+squares
        target_q = density*root3*(100cm)^2/(2cell_vertices)
        isfinite(target_q) && 0<target_q<=1e10 || continue
        best_error = Inf; best_p = 1; best_q = 0; best_norm = 1
        for a in 0:ceil(Int,sqrt(target_q))+1
            bfloat = (-a+sqrt(max(0.0,4target_q-3a^2)))/2
            for b in (clamp(floor(Int,bfloat),0,a),clamp(ceil(Int,bfloat),0,a))
                candidate = a*a+a*b+b*b
                candidate==0 && continue
                discrepancy = abs(candidate-target_q)
                if discrepancy<best_error
                    best_error=discrepancy; best_p=a; best_q=b; best_norm=candidate
                end
            end
        end
        relative_error = best_error/target_q
        acceptable = relative_error<=tolerance
        # Prefer the deepest acceptable level. Until one exists, retain the
        # closest candidate, so an unattainable target is never hidden.
        if k==0 || acceptable || (!found_acceptable && relative_error<selected_error)
            k=depth; p=best_p; q=best_q; norm2=best_norm
            vertices_per_cell=cell_vertices; selected_error=relative_error
        end
        found_acceptable |= acceptable
    end
    k>0 || throw(ArgumentError("Requested circumference/density is outside the supported numerical range."))
    coarse_width = sqrt(Float64(norm2))
    final_width = coarse_width*inflation^k
    spacing_mm = 1000cm/final_width
    final_height = 1000hm/spacing_mm
    mean_density = 2.0*vertices_per_cell*norm2/(root3*(100cm)^2)
    if automatic && !found_acceptable
        @warn "No compatible construction meets the 0.1% mean-density tolerance; using the closest supported result." target_density_per_in2=requested_density mean_density_per_in2=mean_density*6.4516
    end
    estimated_n = mean_density*10000cm*hm
    isfinite(estimated_n) && estimated_n<=1e9 || throw(ArgumentError("Requested rectangle exceeds one billion estimated points."))
    co, si = (p+q/2)/coarse_width, root3*q/(2coarse_width)

    # A vertex (a,b,c,d) means ((a+b*sqrt(3))/2,(c+d*sqrt(3))/2).
    V = NTuple{4,Int64}; F = NTuple{4,V}
    zero_v = (Int64(0),Int64(0),Int64(0),Int64(0))
    directions = V[(2,0,0,0),(0,1,1,0),(1,0,0,1),(0,0,2,0),
                   (-1,0,0,1),(0,-1,1,0),(-2,0,0,0),(0,-1,-1,0),
                   (-1,0,0,-1),(0,0,-2,0),(1,0,0,-1),(0,1,-1,0)]
    templates = (F[],F[])
    template_vertices = (V[],V[])
    for variant in 0:1
        hexagon = [directions[1+variant+2j] for j in 0:5]
        outer_before = V[]; outer_after = V[]
        for j in 1:6
            a,b = hexagon[j],hexagon[mod1(j+1,6)]
            edge = b.-a
            normal = (edge[3],edge[4],-edge[1],-edge[2])
            aa,bb = a.+normal,b.+normal
            push!(templates[variant+1],(zero_v,a,b,zero_v),(a,aa,bb,b))
            push!(outer_after,aa); push!(outer_before,bb)
        end
        for j in 1:6
            a,b,c = hexagon[j],outer_before[mod1(j-1,6)],outer_after[j]
            push!(templates[variant+1],(a,b,c,a))
        end
        append!(template_vertices[variant+1],sort!(collect(Set(v for f in templates[variant+1] for v in f))))
    end

    # Enumerate one representative of each coarse face around the cylinder.
    lower, upper = -3.0, final_height/inflation^k+3
    imin=Inf; imax=-Inf; jmin=Inf; jmax=-Inf
    for u in (-2.0,coarse_width+2), v in (lower-2,upper+2)
        x,y = co*u-si*v,si*u+co*v
        ii,jj = x-y/root3,2y/root3
        imin=min(imin,ii); imax=max(imax,ii); jmin=min(jmin,jj); jmax=max(jmax,jj)
    end
    faces = F[]
    period = (Int64(2p+q),Int64(0),Int64(0),Int64(q))
    for j in floor(Int,jmin)-1:ceil(Int,jmax)+1, i in floor(Int,imin)-1:ceil(Int,imax)+1
        a = (Int64(2i+j),Int64(0),Int64(0),Int64(j))
        b,c,d = a.+directions[1],a.+directions[3],a.+directions[1].+directions[3]
        for face in ((a,b,c,a),(b,d,c,b))
            sx = face[1][1]+face[2][1]+face[3][1]
            sy = face[1][4]+face[2][4]+face[3][4]
            projection = (2p+q)*sx+3q*sy
            0<=projection<12norm2 || continue
            axis = (-si*sx+co*root3*sy)/6
            lower<=axis<=upper && push!(faces,face)
        end
    end
    final_vertices = Set{V}()
    for generation in 1:k
        # Canonicalize shared parent vertices modulo the current exact period.
        parent_vertices = Set{V}()
        period_length = coarse_width*inflation^(generation-1)
        pu = Int128(period[1])^2+3Int128(period[2])^2+Int128(period[3])^2+3Int128(period[4])^2
        pv = 2(Int128(period[1])*period[2]+Int128(period[3])*period[4])
        for face in faces, vertex in face
            x,y = (vertex[1]+root3*vertex[2])/2,(vertex[3]+root3*vertex[4])/2
            t = (co*x+si*y)/period_length
            m = floor(Int,t)
            near = round(Int,t)
            if abs(t-near)<1e-8
                da = Int128(vertex[1])*period[1]+3Int128(vertex[2])*period[2]+Int128(vertex[3])*period[3]+3Int128(vertex[4])*period[4]-near*pu
                db = Int128(vertex[1])*period[2]+Int128(vertex[2])*period[1]+Int128(vertex[3])*period[4]+Int128(vertex[4])*period[3]-near*pv
                sg = da==0 ? sign(db) : db==0 ? sign(da) : sign(da)==sign(db) ? sign(da) : sign(da)*sign(da^2-3db^2)
                m = sg<0 ? near-1 : near
            end
            push!(parent_vertices,vertex.-m.*period)
        end
        period = (2period[1]+3period[2],period[1]+2period[2],2period[3]+3period[4],period[3]+2period[4])
        period_length *= inflation
        pu = Int128(period[1])^2+3Int128(period[2])^2+Int128(period[3])^2+3Int128(period[4])^2
        pv = 2(Int128(period[1])*period[2]+Int128(period[3])*period[4])
        lower,upper = -3.0,final_height/inflation^(k-generation)+3
        new_faces = F[]; candidates = F[]
        # Last step emits points directly; it never materializes all final tiles.
        for source in Iterators.flatten((parent_vertices,faces))
            empty!(candidates)
            if source isa V
                a = source
                center = (2a[1]+3a[2],a[1]+2a[2],2a[3]+3a[4],a[3]+2a[4])
                h = UInt64(seed) ⊻ (UInt64(generation)*0x9e3779b97f4a7c15)
                for coefficient in a
                    h += reinterpret(UInt64,coefficient)+0x9e3779b97f4a7c15
                    h = (h ⊻ (h>>30))*0xbf58476d1ce4e5b9
                    h = (h ⊻ (h>>27))*0x94d049bb133111eb
                    h ⊻= h>>31
                end
                variant = Int(h & 1)+1
                if generation==k
                    for offset in template_vertices[variant]
                        v = center.+offset; push!(candidates,(v,v,v,v))
                    end
                else
                    for f in templates[variant]
                        push!(candidates,(f[1].+center,f[2].+center,f[3].+center,f[4].+center))
                    end
                end
            else
                f = source
                sides = f[1]==f[4] ? 3 : 4
                summed = f[1].+f[2].+f[3]
                sides==4 && (summed=summed.+f[4])
                if sides==3
                    offset = (summed[1]+summed[2],summed[1]÷3+summed[2],summed[3]+summed[4],summed[3]÷3+summed[4])
                    a,b,c = offset.-f[1],offset.-f[2],offset.-f[3]
                    push!(candidates,(a,b,c,a))
                else
                    offset = ((summed[1]+3summed[2])÷4,(summed[1]+summed[2])÷4,(summed[3]+3summed[4])÷4,(summed[3]+summed[4])÷4)
                    square = (f[1].+offset,f[2].+offset,f[3].+offset,f[4].+offset)
                    push!(candidates,square)
                    for j in 1:4
                        a,b = square[j],square[mod1(j+1,4)]
                        e = b.-a; normal = (e[3],e[4],-e[1],-e[2])
                        apex = ((a[1]+b[1]+3normal[2])÷2,(a[2]+b[2]+normal[1])÷2,(a[3]+b[3]+3normal[4])÷2,(a[4]+b[4]+normal[3])÷2)
                        push!(candidates,(a,apex,b,a))
                    end
                end
            end
            for candidate in candidates
                # Process faces at intermediate stages, individual vertices at the last.
                count = generation==k ? (candidate[1]==candidate[2] ? 1 : candidate[1]==candidate[4] ? 3 : 4) : 1
                for j in 1:count
                    f = generation==k ? (candidate[j],candidate[j],candidate[j],candidate[j]) : candidate
                    sides = generation==k ? 1 : f[1]==f[4] ? 3 : 4
                    total = generation==k ? f[1] : f[1].+f[2].+f[3]
                    sides==4 && (total=total.+f[4])
                    x,y = (total[1]+root3*total[2])/(2sides),(total[3]+root3*total[4])/(2sides)
                    axis = -si*x+co*y
                    if generation==k
                        -1e-9<=axis<=final_height+1e-9 || continue
                    else
                        lower<=axis<=upper || continue
                    end
                    t = (co*x+si*y)/period_length
                    m = floor(Int,t); near = round(Int,t)
                    if abs(t-near)<1e-8
                        da = Int128(total[1])*period[1]+3Int128(total[2])*period[2]+Int128(total[3])*period[3]+3Int128(total[4])*period[4]-near*sides*pu
                        db = Int128(total[1])*period[2]+Int128(total[2])*period[1]+Int128(total[3])*period[4]+Int128(total[4])*period[3]-near*sides*pv
                        sg = da==0 ? sign(db) : db==0 ? sign(da) : sign(da)==sign(db) ? sign(da) : sign(da)*sign(da^2-3db^2)
                        m = sg<0 ? near-1 : near
                    end
                    if generation==k
                        push!(final_vertices,f[1].-m.*period)
                    else
                        offset = m.*period
                        push!(new_faces,(f[1].-offset,f[2].-offset,f[3].-offset,f[4].-offset))
                    end
                end
            end
        end
        faces = new_faces
    end
    ordered = sort!(collect(final_vertices))
    coordinates_mm = Matrix{Float64}(undef,length(ordered),2)
    for (i,v) in enumerate(ordered)
        x,y = (v[1]+root3*v[2])/2,(v[3]+root3*v[4])/2
        u = co*x+si*y
        abs(u)<1e-8 && (u=0.0)
        coordinates_mm[i,1] = mod(u*spacing_mm,1000cm)
        coordinates_mm[i,2] = clamp((-si*x+co*y)*spacing_mm,0.0,1000hm)
    end
    info = (target_density_per_in2=requested_density,mean_density_per_in2=mean_density*6.4516,
            density_error_percent=100(mean_density/density-1),
            cropped_density_per_in2=6.4516*length(ordered)/(10000cm*hm),
            number_of_dots=length(ordered),refinement_steps=k,lattice_period=(p,q),
            automatic_refinement=automatic,mean_density_tolerance_percent=100tolerance,
            mean_density_tolerance_met=selected_error<=tolerance,
            seed=seed,rectangle_size_m=(cm,hm))
    return coordinates_mm,spacing_mm,info
end
