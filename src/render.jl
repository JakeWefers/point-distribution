# Small grayscale rasterizer and PNG encoder. Zlib_jll ships with official Julia
# distributions, so PNG output needs no plotting framework or Python installation.
function png_chunk(io,name,data::Vector{UInt8})
    tag=Vector{UInt8}(codeunits(name))
    write(io,hton(UInt32(length(data))));write(io,tag);write(io,data)
    crc=ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),0,tag,length(tag))
    crc=ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),crc,data,length(data))
    write(io,hton(UInt32(crc)))
end
function write_png(path,image::Matrix{UInt8})
    width,height=size(image) # first index is x, for contiguous PNG scanlines
    raw=Vector{UInt8}(undef,(width+1)*height)
    for y in 1:height
        start=(y-1)*(width+1)+1;raw[start]=0 # PNG filter: None
        copyto!(raw,start+1,image,(y-1)*width+1,width)
    end
    capacity=ccall((:compressBound,Zlib_jll.libz),Culong,(Culong,),length(raw))
    compressed=Vector{UInt8}(undef,capacity);n=Ref{Culong}(capacity)
    code=ccall((:compress2,Zlib_jll.libz),Cint,(Ptr{UInt8},Ref{Culong},Ptr{UInt8},Culong,Cint),compressed,n,raw,length(raw),6)
    code==0 || error("zlib compression failed: $code")
    resize!(compressed,n[])
    open(path,"w") do io
        write(io,UInt8[137,80,78,71,13,10,26,10])
        header=IOBuffer();write(header,hton(UInt32(width)),hton(UInt32(height)),UInt8[8,0,0,0,0])
        png_chunk(io,"IHDR",take!(header))
        png_chunk(io,"IDAT",compressed)
        png_chunk(io,"IEND",UInt8[])
    end
end
@inline function blend!(im,x,y,color,coverage)
    alpha=clamp(coverage,0.,1.)
    @inbounds im[x,y]=round(UInt8,(1-alpha)*im[x,y]+alpha*color)
end
function disk!(im,c::Point,radius,color)
    width,height=size(im)
    for y in max(1,floor(Int,c[2]-radius)):min(height,ceil(Int,c[2]+radius+1))
        for x in max(1,floor(Int,c[1]-radius)):min(width,ceil(Int,c[1]+radius+1))
            blend!(im,x,y,color,radius+.5-hypot(x-.5-c[1],y-.5-c[2]))
        end
    end
end
function segment!(im,a::Point,b::Point,line_width,color)
    w,h=size(im);r=line_width/2;dx=b[1]-a[1];dy=b[2]-a[2];den=dx*dx+dy*dy
    den>0 || return
    for y in max(1,floor(Int,min(a[2],b[2])-r)):min(h,ceil(Int,max(a[2],b[2])+r+1))
        for x in max(1,floor(Int,min(a[1],b[1])-r)):min(w,ceil(Int,max(a[1],b[1])+r+1))
            t=clamp(((x-.5-a[1])*dx+(y-.5-a[2])*dy)/den,0.,1.)
            d=hypot(x-.5-a[1]-t*dx,y-.5-a[2]-t*dy)
            d<=r+.5 && blend!(im,x,y,color,r+.5-d)
        end
    end
end
"""Render a viewport from geometric coordinates, rather than enlarging a bitmap."""
function render_image(ideal,noisy,edges,width,height;overlay=false,view=(0.,0.,Float64(width),Float64(height)),dot_radius=7.5,edge_width=2.,gray=205)
    width>0 && height>0 || throw(ArgumentError("Image dimensions must be positive"))
    width*height<=400_000_000 || error("Image exceeds 400 megapixels; reduce spacing_px or disable PNG")
    image=fill(UInt8(255),width,height)
    sx=width/view[3];sy=height/view[4];s=min(sx,sy)
    project(p)=((p[1]-view[1])*sx,(p[2]-view[2])*sy)
    if overlay
        for (i,j) in edges;segment!(image,project(ideal[i]),project(ideal[j]),edge_width*s,gray);end
        for p in ideal;disk!(image,project(p),dot_radius*s,gray);end
    end
    for p in noisy;disk!(image,project(p),dot_radius*s,0);end
    image
end
function write_svg(path,ideal,noisy,edges,width,height;overlay=false,dot_radius=7.5,edge_width=2.,gray=205)
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io,"<title>Schlottmann square-triangle tiling with sampled point noise</title>")
        println(io,"<rect width=\"$width\" height=\"$height\" fill=\"white\"/>")
        if overlay
            println(io,"<g stroke=\"rgb($gray,$gray,$gray)\" stroke-width=\"$edge_width\" fill=\"none\">")
            for (i,j) in edges
                a=ideal[i];b=ideal[j]
                @printf(io,"<path d=\"M %.9f,%.9f L %.9f,%.9f\"/>\n",a[1],a[2],b[1],b[2])
            end
            println(io,"</g>")
        end
        for (points,color) in (overlay ? [(ideal,"rgb($gray,$gray,$gray)"),(noisy,"black")] : [(noisy,"black")])
            println(io,"<g fill=\"$color\">")
            for p in points;@printf(io,"<circle cx=\"%.9f\" cy=\"%.9f\" r=\"%.4f\"/>\n",p[1],p[2],dot_radius);end
            println(io,"</g>")
        end
        println(io,"</svg>")
    end
end
