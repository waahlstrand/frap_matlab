using ImageFiltering: imfilter!, imfilter, Kernel, centered
using ImageTransformations: imresize!, interpolate, BSpline, Linear, Constant
using Interpolations: interpolate
using Flux: SamePad, MeanPool
struct ROI
    shape::String
    x::Int64
    y::Int64
    r::Union{Missing,Float64}
    lx::Union{Missing,Float64}
    ly::Union{Missing,Float64}

    function ROI(shape::String,
                 x::Int64,
                 y::Int64,
                 r::Union{Missing,Float64},
                 lx::Union{Missing,Float64},
                 ly::Union{Missing,Float64})

        if shape == "circle"
            new(shape, x, y, r, missing, missing)
        elseif shape == "rectangle"
            new(shape, x, y, missing, lx, ly)
        else
            #throw(ArgumentError(shape, "Field 'shape' must be either 'circle' or 'rectangle'".))
        end
    end # function

end # struct

function bleach!(c::Array{Float64, 2}, masks::Array)
    for mask in masks
        c = c .* mask
    end

    return c
end


function create_imaging_bleach_mask(β::Float64, n_pixels::Int64, n_pad_pixels::Int64)

    mask    = ones((n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels))
    mask[n_pad_pixels+1:end-n_pad_pixels, n_pad_pixels+1:end-n_pad_pixels] .= β

    return mask

end

function create_bleach_mask(α::Float64,
                            γ::Float64,
                            n_pixels::Int64,
                            n_pad_pixels::Int64,
                            bleach_region::ROI)

    upsampling_factor = 15; # Needs to be a multiple of 3 due the 'box' method in imresize.

    lb_x, lb_y, ub_x, ub_y = get_bounds(bleach_region, γ)

    x = range(lb_x-1, stop=ub_x, length=upsampling_factor*(ub_x-lb_x+1))
    y = range(lb_y-1, stop=ub_y, length=upsampling_factor*(ub_y-lb_y+1))

    X, Y = meshgrid(x, y)

    small_mask = ones(size(X))
    inds = get_bleach_region_index(bleach_region, X, Y)

    small_mask[inds] .= α
    if γ > 0.0
        σ = upsampling_factor * γ
        ℓ = convert(Int64, 4 * ceil(2 * upsampling_factor * γ) + 1 )


        #kernel = centered(Kernel.gaussian((σ,), (ℓ,)))

        #small_mask = imfilter(small_mask, kernel, "replicate")
    end

    pool = MeanPool((upsampling_factor, upsampling_factor); pad=0, stride=(upsampling_factor, upsampling_factor))
    resized_mask = pool(reshape(small_mask, (size(small_mask)..., 1, 1)))

    dims = (max(size(resized_mask)...),max(size(resized_mask)...))
    resized_mask = reshape(resized_mask, dims)
    # Use multi-bilinear interpolation
    #small_mask_size = *(size(small_mask)...)
    #interpolation = interpolate(small_mask, BSpline(Linear()))

    # Pre-allocate the new resized mask (necessary due to non-mature package if we want to
    # specify the interpolation ourselves). Currently identical to
    # imresize(small_mask, (ub_x - lb_x + 1, ub_y - lb_y +1))
    #resized_mask = similar(small_mask, (ub_x -  lb_x + 1, ub_y - lb_y +1))
    #resized_mask = imresize!(resized_mask, interpolation)

    mask    = ones((n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels))
    mask[n_pad_pixels+lb_x:n_pad_pixels+ub_x, n_pad_pixels+lb_y:n_pad_pixels+ub_y] .= resized_mask

    return mask

end

function get_bounds(bleach_region::ROI, γ::Float64)

    factor = 8

    if bleach_region.shape == "circle"

        lb_x = floor(0.5 + bleach_region.x - bleach_region.r - factor*γ)
        ub_x = ceil(0.5 + bleach_region.x + bleach_region.r + factor*γ)
        lb_y = floor(0.5 + bleach_region.y - bleach_region.r - factor*γ)
        ub_y = ceil(0.5 + bleach_region.y + bleach_region.r + factor*γ)

    elseif bleach_region.shape == "rectangle"

        lb_x = floor(0.5 + bleach_region.x - 0.5 * bleach_region.lx - factor*γ)
        ub_x = ceil(0.5 + bleach_region.x + 0.5 * bleach_region.lx + factor*γ)
        lb_y = floor(0.5 + bleach_region.y - 0.5 * bleach_region.ly - factor*γ)
        ub_y = ceil(0.5 + bleach_region.y + 0.5 * bleach_region.ly + factor*γ)
    else
        #throw(ArgumentError(bleach_region, "Field 'shape' must be either 'circle' or 'rectangle'".))
    end

    lb_x, lb_y, ub_x, ub_y = map(x -> convert(Int64, x), (lb_x, lb_y, ub_x, ub_y))

    lb_x -= 1
    lb_y -= 1
    ub_x += 1
    ub_y += 1

    return lb_x, lb_y, ub_x, ub_y
end

function get_bleach_region_index(bleach_region::ROI,
                                 X::Array{Float64, 2},
                                 Y::Array{Float64, 2})

    if bleach_region.shape == "circle"

        # Make a circular cutout of a matrix
        idx_bleach = (( X .- bleach_region.x ).^2 + ( Y .- bleach_region.y ).^2
                        .<= bleach_region.r^2)

    elseif bleach_region.shape == "rectangle"

        # Make a rectangular cutout of a matrix
        idx_bleach = ((X .>= bleach_region.x .- 0.5 .* bleach_region.lx) .&
                      (X .<= bleach_region.x .+ 0.5 .* bleach_region.lx) .&
                      (Y .>= bleach_region.y .- 0.5 .* bleach_region.ly) .&
                      (Y .<= bleach_region.y .+ 0.5 .* bleach_region.ly))
    end

    return idx_bleach
end


function meshgrid(x, y)

    # List comprehensions are very quick
    X = [i for i in x, j in y]
    Y = [j for i in x, j in y]

    return X, Y
end # function
