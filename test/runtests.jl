pushfirst!(LOAD_PATH,joinpath(@__DIR__,".."))
using Test, SchlottmannDots, TOML, Zlib_jll
const S=SchlottmannDots

@testset "Actual substitution and 10,000-point geometry" begin
    rules=load_rules()
    @test length.(rules.children)==[13,16,14,29,29]
    @test rules.inflation≈2+sqrt(3)
    patch=generate_patch(4;rules=rules,verbose=false)
    @test length(patch)==84198
    v=S.vertices(patch,rules)
    @test length(v)==55321 # also catches signed-zero duplicates
    crop=crop_points(v,10000)
    @test length(crop.points)==10000
    coverage=S.verify_crop(patch,rules,crop)
    @test coverage["covered_area"]≈coverage["rectangle_area"] atol=1e-6
    @test coverage["largest_square_cluster_in_checked_interior"]==2
    @test coverage["largest_triangle_cluster_in_checked_interior"]==6
    original=[S.scale(p,100.) for p in crop.points]
    nn,edges=nearest_and_edges(original;radius=120.,edge_length=100.)
    @test all(abs.(nn.-100).<1e-6)
    @test length(unique(Iterators.flatten(edges)))==10000
    # Independent agreement with the reference geometry from the approved image.
    reference=readlines(joinpath(@__DIR__,"..","examples","approved-coordinates.csv"))
    for i in eachindex(original)
        a=parse.(Float64,split(reference[i+1],','))
        p=((crop.points[i][1]+crop.halfwidth)*100+18,(crop.points[i][2]+crop.halfheight)*100+18)
        @test S.distance(p,(a[2],a[3]))<1e-6
    end
    noisy,d=perturb(original;noise_percent=10.,seed=20261006)
    noisy2,d2=perturb(original;noise_percent=10.,seed=20261006)
    @test noisy==noisy2 && d==d2
    @test maximum(hypot.(first.(d),last.(d)))<=10.
    @test abs(sum(p[1]^2+p[2]^2 for p in d)/length(d)-50.)<2.
    nnp,_=nearest_and_edges(noisy;radius=120.000001)
    @test minimum(nnp)>=80-1e-6 && maximum(nnp)<=120+1e-6
    @test first(perturb(original;noise_percent=0.))==original
    @test_throws ArgumentError perturb(original;noise_percent=-1.)
    @test_throws ArgumentError crop_points(v,0)
    @test_throws ArgumentError crop_points(v,100;aspect_ratio=-1.)
    rect=crop_points(v,1000;aspect_ratio=2.)
    @test rect.halfwidth/rect.halfheight≈2.
    @test length(rect.points)==1000
    @test S.verify_crop(patch,rules,rect)["covered_area"]≈4rect.halfwidth*rect.halfheight
end

include("cylinder.jl")

@testset "Spatial search versus exhaustive distances" begin
    grid=S.Point[(Float64(x*100),Float64(y*100)) for x in 0:8 for y in 0:7]
    noisy,_=perturb(grid;noise_percent=10.,seed=61)
    nn,_=nearest_and_edges(noisy;radius=120.)
    brute=[minimum(S.distance(p,q) for (j,q) in enumerate(noisy) if i!=j) for (i,p) in enumerate(noisy)]
    @test nn≈brute
end

@testset "PNG round trip and SVG coordinates" begin
    mktempdir() do dir
        im=reshape(UInt8.(0:255),16,16)
        path=joinpath(dir,"test.png");S.write_png(path,im)
        io=IOBuffer(read(path))
        @test read(io,8)==UInt8[137,80,78,71,13,10,26,10]
        payload=UInt8[]
        while !eof(io)
            n=ntoh(read(io,UInt32));name=String(read(io,4));data=read(io,n);crc=ntoh(read(io,UInt32))
            tag=Vector{UInt8}(codeunits(name))
            check=ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),0,tag,length(tag))
            check=ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),check,data,length(data))
            @test UInt32(check)==crc
            name=="IDAT" && append!(payload,data)
        end
        raw=Vector{UInt8}(undef,17*16);n=Ref{Culong}(length(raw))
        code=ccall((:uncompress,Zlib_jll.libz),Cint,(Ptr{UInt8},Ref{Culong},Ptr{UInt8},Culong),raw,n,payload,length(payload))
        @test code==0
        for y in 1:16
            @test raw[(y-1)*17+1]==0
            @test raw[(y-1)*17+2:y*17]==im[:,y]
        end
        original=S.Point[(10.,10.),(110.,10.)];noisy=S.Point[(12.,11.),(107.,14.)]
        svg=joinpath(dir,"test.svg")
        S.write_svg(svg,original,noisy,[(1,2)],130,30;overlay=true)
        s=read(svg,String)
        @test length(collect(eachmatch(r"<circle ",s)))==4
        @test occursin("12.000000000",s) && occursin("107.000000000",s)
    end
end
