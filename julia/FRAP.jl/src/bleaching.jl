using ImageFiltering: imfilter!, Kernel, centered
using ImageTransformations: imresize!, interpolate, BSpline, Linear
using Interpolations: interpolate

struct ROI{T}
    shape::String
    x::T
    y::T
    r::Union{Missing,T}
    lx::Union{Missing,T}
    ly::Union{Missing,T}

    function ROI(shape::String,
                 x::T,
                 y::T,
                 r::Union{Missing,T},
                 lx::Union{Missing,T},
                 ly::Union{Missing,T}) where T

        if shape == "circle"
            new{T}(shape, x, y, r, missing, missing)
        elseif shape == "rectangle"
            new{T}(shape, x, y, missing, lx, ly)
        else
            throw(ArgumentError(shape, "Field 'shape' must be either 'circle' or 'rectangle'".))
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

    upsampling_factor = 15.0; % Needs to be a multiple of 3 due the 'box' method in imresize.

    if bleach_region.shape == "circle"
        lb_x = floor(0.5 + bleach_region.x - bleach_region.r - 8 * γ)
        ub_x = ceil(0.5 + bleach_region.x + bleach_region.r + 8 * γ)
        lb_y = floor(0.5 + bleach_region.y - bleach_region.r - 8 * γ)
        ub_y = ceil(0.5 + bleach_region.y + bleach_region.r + 8 * γ)
    elseif bleach_region.shape == "rectangle"
        lb_x = floor(0.5 + bleach_region.x - 0.5 * bleach_region.lx - 8 * γ)
        ub_x = ceil(0.5 + bleach_region.x + 0.5 * bleach_region.lx + 8 * γ)
        lb_y = floor(0.5 + bleach_region.y - 0.5 * bleach_region.ly - 8 * γ)
        ub_y = ceil(0.5 + bleach_region.y + 0.5 * bleach_region.ly + 8 * γ)
    else
        throw(ArgumentError(bleach_region, "Field 'shape' must be either 'circle' or 'rectangle'".))
    end

    lb_x = lb_x - 1
    ub_x = ub_x + 1
    lb_y = lb_y - 1
    ub_y = ub_y + 1

    #xx = linspace(lb_x - 1, ub_x, upsampling_factor * (ub_x - lb_x + 1));
    #yy = linspace(lb_y - 1, ub_y, upsampling_factor * (ub_y - lb_y + 1));

    x = range(-n_pixels ÷ 2, length=n_pixels)
    y = range(-n_pixels ÷ 2, length=n_pixels)

    # List comprehensions are very quick
    X = [j for i in x, j in y]
    Y = [i for i in x, j in y]

    small_mask = ones(size(X))

    if γ > 0.0
        σ = upsampling_factor * γ
        ℓ = 4 * ceil(2 * upsampling_factor * γ) + 1


        kernel = centered(Kernel.gaussian(σ, ℓ))


        small_mask = imfilter!(small_mask, kernel, "replicate")
    end

    # Use multi-bilinear interpolation
    interpolation = interpolate(small_mask, BSpline(Linear()))

    # Pre-allocate the new resized mask (necessary due to non-mature package if we want to
    # specify the interpolation ourselves). Currently identical to
    # imresize(small_mask, (ub_x - lb_x + 1, ub_y - lb_y +1))
    resized_mask = similar(small_mask, (ub_x - lb_x + 1, ub_y - lb_y +1))
    resized_mask = imresize!(resized_mask, interpolate(small_mask, BSpline(Linear())))

end
