module Bleach

    include("Utils.jl")

    using FFTW: fft, ifft, ifftshift
    using ImageFiltering: imfilter!, Kernel, centered
    using ImageTransformations: imresize!, interpolate, BSpline, Linear
    using Interpolations: interpolate
    using .Utils: Concentration

    export evolve, evolve!, create_fourier_grid, create_imaging_bleach_mask

    function evolve!(C::Concentration, ξ::Array{Float64, 2}, D::Array{Float64, 1}, δt::Float64, masks::Array, frames::Array{Int64, 1})

        # Use simple Markov memory, one timestep
        c = C[:, :, frames[begin]]

        for frame in frames[begin+1:end]
            
            # Make one step of length δt
            c.mobile = step!(c.mobile, ξ, D, δt)

            # Apply imaging bleach masks
            # c.mobile    = bleach!(c.mobile, masks)
            # c.immobile  = bleach!(c.immobile, masks)

            # Save time evolution
            C[:, :, frame] = c

        end
        
        return C
    end

    function evolve(C::Concentration, ξ::Array{Float64, 2}, D::Array{Float64, 1}, δt::Float64, masks::Array, n_frames::Int64)

        # Use simple Markov memory, one timestep
        c = C[:, :]
        C_to_save = similar(C, (size(C)..., n_frames))

        for frame in 1:n_frames
            
            # Make one step of length δt
            c.mobile = step!(c.mobile, ξ, D, δt)

            # Apply imaging bleach masks
            # c.mobile    = bleach!(c.mobile, masks)
            # c.immobile  = bleach!(c.immobile, masks)

            # Save time evolution
            C_to_save[:, :, frame] = c

        end
        
        return C_to_save
    end


    function bleach!(c::Array{Float64, 2}, masks::Array)
        for mask in masks
            c = c .* mask
        end

        return c
    end

    function step!(c::Array{Float64, 2}, ξ²::Array{Float64, 2}, D::Float64, δt::Float64)

        dims = [1, 2]

        ĉ = fft(c, dims)
        ĉ = exp.(-D .* ξ² .* δt) .* ĉ
        c = abs.(ifft(ĉ, dims))

        return c

    end

    function create_fourier_grid(n_pixels::Int64)
        # TODO: Check dimensions of grid
        # It isn't centralized, but perhaps that is not a problem?

        # Julia has no ndgrid or meshgrid function, out of principle
        # it would seem. Not "Julian" enough.
        x = range(-n_pixels ÷ 2, length=n_pixels)
        y = range(-n_pixels ÷ 2, length=n_pixels)

        # List comprehensions are very quick
        X = [j for i in x, j in y]
        Y = [i for i in x, j in y]

        X *= (2π / n_pixels) 
        Y *= (2π / n_pixels) 
        
        # Centre the transform
        ξ² = ifftshift(X.^2 + Y.^2)

        return ξ²
    end

    function create_imaging_bleach_mask(β::Float64, n_pixels::Int64, n_pad_pixels::Int64)

        mask    = ones((n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels))
        mask[n_pad_pixels+1:end-n_pad_pixels, n_pad_pixels+1:end-n_pad_pixels] .= β

        return mask

    end

    function create_bleach_mask(α::Float64, γ::Float64, n_pixels::Int64, n_pad_pixels::Int64)

        upsampling_factor = 1.0
        lb_x = 0
        lb_y = 0
        ub_x = 0
        ub_y = 0

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


            small_mask = imfilter!(small_mask, kernel, [border="replicate"])
        end

        # Use multi-bilinear interpolation
        interpolation = interpolate(small_mask, BSpline(Linear()))

        # Pre-allocate the new resized mask (necessary due to non-mature package if we want to
        # specify the interpolation ourselves). Currently identical to
        # imresize(small_mask, (ub_x - lb_x + 1, ub_y - lb_y +1))
        resized_mask = similar(small_mask, (ub_x - lb_x + 1, ub_y - lb_y +1))
        resized_mask = imresize!(resized_mask, interpolate(small_mask, BSpline(Linear())))
        
    end
end