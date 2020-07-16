module Bleach

    include("Utils.jl")

    using FFTW: fft, ifft, ifftshift
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

        # Julia has no ndgrid or meshgrid function, out of principle
        # it would seem. Not "Julian" enough.
        x = range(-n_pixels ÷ 2, length=n_pixels+1)
        y = range(-n_pixels ÷ 2, length=n_pixels+1)

        # List comprehensions are very quick
        X = [i for i in x, j in y]
        Y = [j for i in x, j in y]

        # "Fourier transform"
        ξ = 2π ./ X
        η = 2π ./ Y
        
        # Centre the transform
        ξ² = ifftshift(ξ.^2 + η.^2)

        return ξ²
    end

    function create_imaging_bleach_mask(β::Float64, n_pixels::Int64, n_pad_pixels::Int64)

        mask    = ones((n_pixels + 2*n_pad_pixels, n_pixels + 2*n_pad_pixels))
        centre  = [n_pad_pixels+1:end-n_pad_pixels, n_pad_pixels+1:end-n_pad_pixels]

        mask[centre...] = β

        return mask

    end
end