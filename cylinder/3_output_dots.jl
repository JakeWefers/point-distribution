using Zlib_jll 

function output_dots(noisy_coordinates_mm::AbstractMatrix, relative_error_bound,
                     rectangle_size_m, radius_mm; spacing_mm, dpi=2000,
                     filename="dots.tif", stripe_rows=256, crop_mm=nothing)
    n = size(noisy_coordinates_mm,1)
    size(noisy_coordinates_mm,2)==2 && n>=1 ||
        throw(ArgumentError("noisy_coordinates_mm must be an N × 2 matrix."))
    all(isfinite,noisy_coordinates_mm) || throw(ArgumentError("Coordinates must be finite."))
    length(rectangle_size_m)==2 || throw(ArgumentError("Use (circumference_m, height_m)."))
    wm, hm = Float64.(rectangle_size_m)
    isfinite(wm) && isfinite(hm) && wm>0 && hm>0 ||
        throw(ArgumentError("Rectangle dimensions must be positive and finite."))
    resolution = Float64(dpi)
    isfinite(resolution) && resolution>0 || throw(ArgumentError("dpi must be positive and finite."))
    pixels_per_mm = resolution/25.4
    wfloat, hfloat = wm*1000pixels_per_mm, hm*1000pixels_per_mm
    0.5<=wfloat<=typemax(Int32) && 0.5<=hfloat<=typemax(Int32) ||
        throw(ArgumentError("Pixel dimensions must lie between 1 and 2^31-1."))
    width, height = round(Int,wfloat), round(Int,hfloat)
    width>=1 && height>=1 || throw(ArgumentError("Image dimensions round to zero pixels."))
    full_width, full_height = width, height
    xscale, yscale = full_width/(1000wm), full_height/(1000hm)
    x_offset, y_offset = 0, 0
    if crop_mm!==nothing
        length(crop_mm)==4 || throw(ArgumentError("Use crop_mm=(x_mm, y_mm, width_mm, height_mm)."))
        crop_x, crop_y, crop_width, crop_height = Float64.(crop_mm)
        all(isfinite,(crop_x,crop_y,crop_width,crop_height)) &&
            crop_x>=0 && crop_y>=0 && crop_width>0 && crop_height>0 &&
            crop_x+crop_width<=1000wm && crop_y+crop_height<=1000hm ||
            throw(ArgumentError("The cutout must have positive dimensions and fit inside the full rectangle."))
        x_offset, y_offset = round(Int,crop_x*xscale), round(Int,crop_y*yscale)
        width = min(round(Int,crop_width*xscale),full_width-x_offset)
        height = min(round(Int,crop_height*yscale),full_height-y_offset)
        width>=1 && height>=1 || throw(ArgumentError("Cutout dimensions round to zero pixels."))
    end
    radius, spacing, err = Float64(radius_mm), Float64(spacing_mm), Float64(relative_error_bound)
    isfinite(radius) && radius>0 || throw(ArgumentError("radius_mm must be positive and finite."))
    isfinite(spacing) && spacing>0 || throw(ArgumentError("spacing_mm must be positive and finite."))
    isfinite(err) && err>=0 || throw(ArgumentError("relative_error_bound must be a nonnegative finite fraction."))
    safe_radius = max(0.0,spacing*(1-err-1e-7)/2)
    radius<=safe_radius || throw(ArgumentError(
        "Nonoverlap cannot be guaranteed from this noise bound. Use radius_mm <= $safe_radius."))
    stripe_rows isa Integer && stripe_rows>0 || throw(ArgumentError("stripe_rows must be a positive integer."))
    extension = lowercase(splitext(filename)[2])
    extension in (".png",".tif",".tiff") || throw(ArgumentError("filename must end in .png, .tif or .tiff."))
    tiff_output = extension!=".png"
    vips_executable = Sys.which("vips")
    if tiff_output && vips_executable===nothing
        for candidate in ("/opt/homebrew/bin/vips","/usr/local/bin/vips")
            if isfile(candidate)
                vips_executable=candidate
                break
            end
        end
        vips_executable===nothing && throw(ArgumentError("TIFF export requires libvips. Install with: brew install vips"))
    end
    destination = abspath(filename)
    isdir(dirname(destination)) || throw(ArgumentError("The output directory must already exist."))
    stride = cld(width,8)+1 # One filter byte followed by packed monochrome pixels.
    rows_per_stripe = min(Int(stripe_rows),height,max(1,div(64*1024*1024,stride)))
    stride*rows_per_stripe <= 512*1024*1024 || throw(ArgumentError("One scanline is too large."))
    order = sortperm(view(noisy_coordinates_mm,:,2))
    pixel_radius = radius*yscale
    radius<500wm || throw(ArgumentError("radius_mm must be smaller than half the circumference."))
    isfinite(pixel_radius^2) || throw(ArgumentError("Pixel radius is too large."))

    gzip_path, temp_io = mktemp(dirname(destination))
    close(temp_io)
    png_path = ""
    tiff_path = ""
    gzip_handle = Ptr{Cvoid}(C_NULL)
    gzinput = nothing
    output = nothing
    try
        gzip_handle = ccall((:gzopen,Zlib_jll.libz),Ptr{Cvoid},(Cstring,Cstring),gzip_path,"wb6")
        gzip_handle!=C_NULL || error("Cannot open temporary compression stream.")
        raw = Vector{UInt8}(undef,stride*rows_per_stripe)
        adler = Culong(1)
        lo, hi = 1, 0
        for first_row in 1:rows_per_stripe:height
            last_row = min(height,first_row+rows_per_stripe-1)
            resize!(raw,stride*(last_row-first_row+1))
            fill!(raw,0xff)
            @inbounds for row in first_row:last_row
                raw[(row-first_row)*stride+1] = 0 # PNG filter type 0.
            end
            @inbounds while lo<=n && noisy_coordinates_mm[order[lo],2]*yscale-y_offset+pixel_radius<first_row-0.5
                lo += 1
            end
            hi = max(hi,lo-1)
            @inbounds while hi<n && noisy_coordinates_mm[order[hi+1],2]*yscale-y_offset-pixel_radius<=last_row-0.5
                hi += 1
            end
            @inbounds for index in lo:hi
                point = order[index]
                cx = mod(noisy_coordinates_mm[point,1],1000wm)*xscale-x_offset
                cy = noisy_coordinates_mm[point,2]*yscale-y_offset
                ylo = max(first_row,ceil(Int,clamp(cy-pixel_radius+0.5,0.0,height+1.0)))
                yhi = min(last_row,floor(Int,clamp(cy+pixel_radius+0.5,0.0,height+1.0)))
                for y in ylo:yhi
                    half_span = sqrt(max(0.0,pixel_radius^2-(y-0.5-cy)^2))*xscale/yscale
                    for wrapped_cx in (cx-full_width,cx,cx+full_width)
                        xlo = max(1,ceil(Int,clamp(wrapped_cx-half_span+0.5,0.0,width+1.0)))
                        xhi = min(width,floor(Int,clamp(wrapped_cx+half_span+0.5,0.0,width+1.0)))
                        xlo<=xhi || continue
                        offset = (y-first_row)*stride+1
                        bfirst, blast = div(xlo-1,8)+1, div(xhi-1,8)+1
                        left_mask = UInt8(0xff) >> mod(xlo-1,8)
                        right_mask = UInt8(0xff) << (7-mod(xhi-1,8))
                        if bfirst==blast
                            raw[offset+bfirst] &= ~(left_mask & right_mask)
                        else
                            raw[offset+bfirst] &= ~left_mask
                            for b in bfirst+1:blast-1
                                raw[offset+b] = 0
                            end
                            raw[offset+blast] &= ~right_mask
                        end
                    end
                end
            end
            adler = ccall((:adler32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),adler,raw,length(raw))
            written = ccall((:gzwrite,Zlib_jll.libz),Cint,(Ptr{Cvoid},Ptr{UInt8},Cuint),gzip_handle,raw,length(raw))
            written==length(raw) || error("PNG compression write failed; check available disk space.")
        end
        status = ccall((:gzclose,Zlib_jll.libz),Cint,(Ptr{Cvoid},),gzip_handle)
        gzip_handle = C_NULL
        status==0 || error("PNG compression close failed.")

        # gzip and PNG's zlib stream both contain raw DEFLATE. Replace gzip's
        # fixed header/trailer with a zlib header and the scanlines' Adler-32.
        gzinput = open(gzip_path,"r")
        gzip_header = read(gzinput,10)
        gzip_header[1:4]==UInt8[0x1f,0x8b,0x08,0x00] || error("Unexpected gzip header.")
        remaining = filesize(gzip_path)-18
        remaining>=0 || error("Incomplete compression stream.")
        header = IOBuffer()
        write(header,hton(UInt32(width)),hton(UInt32(height)),UInt8[1,0,0,0,0])
        physical = IOBuffer()
        ppm = round(Int,resolution/0.0254)
        0<ppm<=typemax(UInt32) || throw(ArgumentError("DPI exceeds PNG metadata range."))
        write(physical,hton(UInt32(ppm)),hton(UInt32(ppm)),UInt8(1))
        prefix = [("IHDR",take!(header)),("pHYs",take!(physical)),
                  ("tEXt",Vector{UInt8}(codeunits("relative_error_bound\0$(repr(err))"))),
                  ("IDAT",UInt8[0x78,0x9c])]
        trailer = IOBuffer()
        write(trailer,hton(UInt32(adler)))
        suffix = [("IDAT",take!(trailer)),("IEND",UInt8[])]
        png_path, output = mktemp(dirname(destination))
        write(output,UInt8[137,80,78,71,13,10,26,10])
        prefix_index, suffix_index = 1, 1
        while true
            if prefix_index<=length(prefix)
                name,data = prefix[prefix_index]
                prefix_index += 1
            elseif remaining>0
                name = "IDAT"
                count = min(remaining,1024*1024)
                data = read(gzinput,count)
                length(data)==count || error("Incomplete compressed data.")
                remaining -= count
            elseif suffix_index<=length(suffix)
                name,data = suffix[suffix_index]
                suffix_index += 1
            else
                break
            end
            tag = Vector{UInt8}(codeunits(name))
            write(output,hton(UInt32(length(data))),tag,data)
            crc = ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),0,tag,length(tag))
            crc = ccall((:crc32,Zlib_jll.libz),Culong,(Culong,Ptr{UInt8},Cuint),crc,data,length(data))
            write(output,hton(UInt32(crc)))
        end
        close(output)
        output = nothing
        close(gzinput)
        gzinput = nothing
        if tiff_output
            tiff_path, temp_tiff_io = mktemp(dirname(destination))
            close(temp_tiff_io)
            # Julia Cmd interpolation passes paths as individual arguments, including spaces.
            # Sequential PNG access avoids an enormous decompressed temporary image.
            run(`$vips_executable tiffsave $(png_path*"[access=sequential]") $tiff_path
                --tile --pyramid --tile-width=512 --tile-height=512
                --compression=deflate --predictor=none --bigtiff --bitdepth=1
                --region-shrink=nearest --vips-concurrency=2`)
            mv(tiff_path,destination;force=true)
            tiff_path = ""
        else
            mv(png_path,destination;force=true)
            png_path = ""
        end
    finally
        gzip_handle==C_NULL || ccall((:gzclose,Zlib_jll.libz),Cint,(Ptr{Cvoid},),gzip_handle)
        output===nothing || close(output)
        gzinput===nothing || close(gzinput)
        isfile(gzip_path) && rm(gzip_path)
        !isempty(png_path) && isfile(png_path) && rm(png_path)
        !isempty(tiff_path) && isfile(tiff_path) && rm(tiff_path)
    end
    return destination
end
