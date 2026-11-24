module CylinderTests
using Test, Zlib_jll
include(joinpath(@__DIR__, "..", "cylinder", "1_tiling_coordinates.jl"))
include(joinpath(@__DIR__, "..", "cylinder", "2_add_noise.jl"))
include(joinpath(@__DIR__, "..", "cylinder", "3_output_dots.jl"))

@testset "Periodic cylinder geometry and bounded noise" begin
    rectangle = (0.025, 0.025)
    points, spacing, info = tiling_coordinates(rectangle;
        density_per_in2=258.064, refinement_steps=1, seed=123)
    again, d2, i2 = tiling_coordinates(rectangle;
        density_per_in2=258.064, refinement_steps=1, seed=123)
    @test points==again && spacing==d2 && info==i2
    @test size(points,1)>2
    @test all(0 .<= points[:,1] .< 25)
    @test all(0 .<= points[:,2] .<= 25)
    @test length(Set(Tuple(points[i,:]) for i in axes(points,1)))==size(points,1)
    closest = Inf
    for i in axes(points,1), j in i+1:size(points,1)
        dx = abs(points[i,1]-points[j,1])
        closest = min(closest, hypot(min(dx,25-dx),points[i,2]-points[j,2]))
    end
    @test closest≈spacing atol=1e-8
    noisy, bound = add_noise(points, rectangle;
        spacing_mm=spacing, noise_percent=10., seed=77)
    repeated, bound2 = add_noise(points, rectangle;
        spacing_mm=spacing, max_shift_mm=0.1spacing, seed=77)
    @test noisy≈repeated && bound≈bound2
    @test bound≈0.2
    dx = mod.(noisy[:,1]-points[:,1].+12.5,25).-12.5
    dy = noisy[:,2]-points[:,2]
    @test maximum(hypot.(dx,dy))<=0.1spacing+1e-10
    @test all(0 .<= noisy[:,1] .< 25)
    zero, b0 = add_noise(points, rectangle; spacing_mm=spacing, noise_percent=0.)
    @test zero==points && b0==0
    @test_throws ArgumentError add_noise(points, rectangle;
        spacing_mm=spacing, noise_percent=-1.)
    @test_throws ArgumentError add_noise(points, rectangle;
        spacing_mm=spacing, noise_percent=10., max_shift_mm=1.)
end

function decoded_png(path)
    io=IOBuffer(read(path)); @test read(io,8)==UInt8[137,80,78,71,13,10,26,10]
    width=height=0; depth=color=-1; payload=UInt8[]
    while !eof(io)
        n=ntoh(read(io,UInt32)); name=String(read(io,4)); data=read(io,n)
        read(io,UInt32)
        if name=="IHDR"
            header=IOBuffer(data)
            width=Int(ntoh(read(header,UInt32))); height=Int(ntoh(read(header,UInt32)))
            depth=Int(read(header,UInt8)); color=Int(read(header,UInt8))
        elseif name=="IDAT"
            append!(payload,data)
        end
    end
    raw=Vector{UInt8}(undef,height*(cld(width,8)+1)); size_out=Ref{Culong}(length(raw))
    code=ccall((:uncompress,Zlib_jll.libz),Cint,
        (Ptr{UInt8},Ref{Culong},Ptr{UInt8},Culong),raw,size_out,payload,length(payload))
    @test code==0 && size_out[]==length(raw)
    (;width,height,depth,color,raw)
end

@testset "Packed 1-bit PNG and cylinder join" begin
    mktempdir() do dir
        # A circle centered just to the right of the join must also appear at
        # the right edge. Another point provides a separate ordinary dot.
        points = [0.05 5.0; 5.0 2.0]
        a=output_dots(points,0.2,(0.01,0.01),0.3;spacing_mm=1.,dpi=254,
            filename=joinpath(dir,"a.png"),stripe_rows=1)
        b=output_dots(points,0.2,(0.01,0.01),0.3;spacing_mm=1.,dpi=254,
            filename=joinpath(dir,"b.png"),stripe_rows=17)
        first=decoded_png(a); second=decoded_png(b)
        @test first==second
        @test (first.width,first.height,first.depth,first.color)==(100,100,1,0)
        stride=cld(first.width,8)+1
        isblack(x,y) = (first.raw[y*stride+2+div(x,8)] & (0x80>>mod(x,8)))==0
        @test isblack(0,49) && isblack(99,49)
        @test !isblack(20,49)
        saved=read(a)
        @test_throws ArgumentError output_dots(points,0.2,(0.01,0.01),0.5;
            spacing_mm=1.,dpi=254,filename=a)
        @test read(a)==saved
    end
end
end
