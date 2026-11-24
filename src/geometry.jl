const VertexKey = NTuple{2,Int64}
const FaceKey = Tuple{Int,NTuple{4,VertexKey}}
vertexkey(p::Point) = (round(Int64,p[1]*1e7),round(Int64,p[2]*1e7))
function facekey(p::Vector{Point})::FaceKey
    v=sort(vertexkey.(p)); n=length(v)
    (n,ntuple(i->i<=n ? v[i] : (Int64(0),Int64(0)),4))
end
function signature(kind,O,r::Rules)
    # Different marked orientations can be equivalent under the substitution.
    sort([(t.kind,vertexkey(transform(O,t.center))...) for t in r.children[kind]])
end

"""Generate a finite patch using actual overlapping five-prototile substitution.
Duplicate tiles are merged; inconsistent types or replacement orientations fail.
"""
function generate_patch(generations::Integer; rules=load_rules(), seed_type=4, verbose=true)
    0<=generations<=6 || throw(ArgumentError("generations must be in 0:6 (memory grows about 14× each generation)"))
    1<=seed_type<=5 || throw(ArgumentError("seed_type must be in 1:5"))
    patch=[Tile(seed_type,(0.0,0.0),IDENTITY)]
    signatures=Dict{Tuple{Int,NTuple{4,Int64}},Vector{Tuple{Int,Int64,Int64}}}()
    function cached_signature(t)
        key=(t.kind,ntuple(i->round(Int64,t.orientation[i]*1e7),4))
        get!(signatures,key) do
            signature(t.kind,t.orientation,rules)
        end
    end
    for generation in 1:generations
        next=Tile[];seen=Dict{FaceKey,Int}();duplicates=0
        for parent in patch, child in rules.children[parent.kind]
            t=Tile(child.kind,add(scale(parent.center,rules.inflation),transform(parent.orientation,child.center)),compose(parent.orientation,child.orientation))
            key=facekey(corners(t,rules))
            if haskey(seen,key)
                old=next[seen[key]];duplicates+=1
                old.kind==t.kind || error("Conflicting tile types in substitution overlap")
                cached_signature(old)==cached_signature(t) || error("Conflicting substitution orientations")
            else
                push!(next,t);seen[key]=length(next)
            end
        end
        patch=next
        verbose && println("Generation $generation: $(length(patch)) tiles; $duplicates shared copies merged")
    end
    patch
end
function vertices(patch,rules)
    v=Point[]
    for tile in patch, p in corners(tile,rules)
        # Normalize signed zero: Julia treats -0.0 and +0.0 as distinct set keys.
        push!(v,(round(p[1];digits=9)+0.0,round(p[2];digits=9)+0.0))
    end
    sort!(unique!(v))
end

"""Select exactly n vertices in a rectangle of the requested width/height ratio.
The physical extent is chosen from vertex positions, not by deleting extra points.
"""
function crop_points(v::Vector{Point},n::Integer; aspect_ratio=1.0, center=(.13718,.29413),angle=.17320508075688773)
    2<=n<length(v) || throw(ArgumentError("Need at least n+1 vertices; increase generations"))
    isfinite(aspect_ratio) && aspect_ratio>0 || throw(ArgumentError("aspect_ratio must be positive and finite"))
    R=rotation(angle); transformed=[transform(R,sub(p,center)) for p in v]
    a=sqrt(aspect_ratio)
    radii=[max(abs(p[1])/a,abs(p[2])*a) for p in transformed]
    order=sortperm(radii)
    radii[order[n]]<radii[order[n+1]] || error("Crop-boundary tie; change crop_center or crop_rotation_radians slightly")
    half=(radii[order[n]]+radii[order[n+1]])/2
    (points=transformed[order[1:n]],halfwidth=half*a,halfheight=half/a,R=R,center=center)
end
function clip_polygon(poly::Vector{Point},axis,bound,greater)
    isempty(poly) && return poly
    result=Point[]
    for i in eachindex(poly)
        a=poly[mod1(i-1,length(poly))];b=poly[i]
        ai=greater ? a[axis]>=bound : a[axis]<=bound
        bi=greater ? b[axis]>=bound : b[axis]<=bound
        if ai!=bi
            t=(bound-a[axis])/(b[axis]-a[axis])
            push!(result,add(a,scale(sub(b,a),t)))
        end
        bi && push!(result,b)
    end
    result
end
function polygon_area(p)
    length(p)<3 && return 0.0
    abs(sum(p[i][1]*p[mod1(i+1,length(p))][2]-p[i][2]*p[mod1(i+1,length(p))][1] for i in eachindex(p)))/2
end

"""Check unit edges and square coverage, and quantify same-shape tile clusters."""
function verify_crop(patch,rules,crop)
    bounds=(crop.halfwidth,crop.halfheight)
    area=0.0;nt=0;ns=0;interior=Vector{Point}[]
    for tile in patch
        raw=corners(tile,rules)
        all(abs(distance(raw[i],raw[mod1(i+1,length(raw))])-1)<1e-7 for i in eachindex(raw)) || error("Non-unit tile edge")
        p=[transform(crop.R,sub(q,crop.center)) for q in raw]
        any(minimum(q[axis] for q in p)>bounds[axis] || maximum(q[axis] for q in p)<-bounds[axis] for axis in 1:2) && continue
        q=p
        for axis in 1:2
            q=clip_polygon(q,axis,-bounds[axis],true)
            q=clip_polygon(q,axis,bounds[axis],false)
        end
        A=polygon_area(q);area+=A
        if A>1e-9
            length(p)==3 ? (nt+=1) : (ns+=1)
        end
        all(abs(q[1])<bounds[1] && abs(q[2])<bounds[2] for q in p) && push!(interior,p)
    end
    target=4prod(bounds)
    abs(area-target)<=max(1e-7,target*1e-10) || error("Crop not covered by the patch; increase generations")
    faces=Dict{Tuple{VertexKey,VertexKey},Vector{Int}}()
    for (i,p) in enumerate(interior), j in eachindex(p)
        a=vertexkey(p[j]);b=vertexkey(p[mod1(j+1,length(p))]);key=a<b ? (a,b) : (b,a)
        push!(get!(faces,key,Int[]),i)
    end
    parent=collect(eachindex(interior))
    function root(i)
        while parent[i]!=i
            parent[i]=parent[parent[i]];i=parent[i]
        end
        i
    end
    for f in values(faces)
        length(f)<=2 || error("More than two tiles share an edge")
        if length(f)==2 && length(interior[f[1]])==length(interior[f[2]])
            parent[root(f[2])]=root(f[1])
        end
    end
    sizes=Dict{Int,Int}()
    for i in eachindex(interior)
        k=root(i);sizes[k]=get(sizes,k,0)+1
    end
    maxt=0;maxs=0
    for (k,n) in sizes
        length(interior[k])==3 ? (maxt=max(maxt,n)) : (maxs=max(maxs,n))
    end
    Dict("covered_area"=>area,"rectangle_area"=>target,"intersecting_triangles"=>nt,"intersecting_squares"=>ns,
         "largest_triangle_cluster_in_checked_interior"=>maxt,"largest_square_cluster_in_checked_interior"=>maxs)
end

"""Exact local distance search. Fails if no neighbor is within the search radius.
When used on the noiseless tiling, pairs one edge length apart reconstruct its edges.
"""
function nearest_and_edges(points::Vector{Point}; radius=120.0,edge_length=nothing)
    radius>0 || throw(ArgumentError("radius must be positive"))
    cells=Dict{Tuple{Int,Int},Vector{Int}}()
    keys=[(floor(Int,p[1]/radius),floor(Int,p[2]/radius)) for p in points]
    for (i,k) in enumerate(keys);push!(get!(cells,k,Int[]),i);end
    nn=fill(Inf,length(points));edges=Tuple{Int,Int}[]
    for (i,(x,y)) in enumerate(keys), dx in -1:1,dy in -1:1
        for j in get(cells,(x+dx,y+dy),Int[])
            i==j && continue
            d=distance(points[i],points[j]);nn[i]=min(nn[i],d)
            if edge_length!==nothing && j>i && abs(d-edge_length)<=1e-6
                push!(edges,(i,j))
            end
        end
    end
    maximum(nn)<=radius+1e-6 || error("Nearest-neighbor search radius too small or crop contains an isolated vertex")
    nn,edges
end
