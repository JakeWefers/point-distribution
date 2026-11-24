module SchlottmannDots
using TOML, Printf, Zlib_jll
export generate_patch, crop_points, perturb, nearest_and_edges, load_rules, run

const Point = NTuple{2,Float64}
const Matrix2 = NTuple{4,Float64} # row-major
const IDENTITY = (1.0, 0.0, 0.0, 1.0)
add(a::Point,b::Point) = (a[1]+b[1],a[2]+b[2])
sub(a::Point,b::Point) = (a[1]-b[1],a[2]-b[2])
scale(a::Point,s::Real) = (a[1]*s,a[2]*s)
transform(A::Matrix2,p::Point) = (A[1]*p[1]+A[2]*p[2], A[3]*p[1]+A[4]*p[2])
compose(A::Matrix2,B::Matrix2) = (A[1]*B[1]+A[2]*B[3], A[1]*B[2]+A[2]*B[4], A[3]*B[1]+A[4]*B[3], A[3]*B[2]+A[4]*B[4])
distance(a::Point,b::Point) = hypot(a[1]-b[1],a[2]-b[2])
rotation(a::Real) = (cos(a),-sin(a),sin(a),cos(a))

struct Tile
    kind::Int
    center::Point
    orientation::Matrix2
end
struct Rules
    inflation::Float64
    base::Vector{Vector{Point}}
    children::Vector{Vector{Tile}}
end
function load_rules(path=joinpath(@__DIR__,"..","data","substitution.toml"))
    d=TOML.parsefile(path)
    base=[[Tuple(Float64.(p))::Point for p in d["base"][string(k)]] for k in 1:5]
    children=[[Tile(c["type"],Tuple(Float64.(c["center"])),Tuple(Float64.(c["orientation"]))) for c in d["rules"][string(k)]] for k in 1:5]
    Rules(d["inflation"],base,children)
end
corners(t::Tile,r::Rules) = [add(t.center,transform(t.orientation,p)) for p in r.base[t.kind]]
include("geometry.jl")
include("noise.jl")
include("render.jl")
include("pipeline.jl")
end
