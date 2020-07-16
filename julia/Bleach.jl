module Bleach

    include("Utils.jl")

    using FFTW: fft, ifft
    using .Utils: Concentration

    export evolve

    function evolve(c_init::Concentration, ξ::Array{Float32, 2}, D::Array{Float32, 1}, δt::Float32, masks::Array, n_frames::Integer)

        # Use simple Markov memory, one timestep
        c = c_init[:,:,end]

        for frame = 1:n_frames
            
            # Make one step in time
            c.mobile = step!(c.mobile, ξ, D, δt)

            # Apply imaging bleach masks
            c.mobile    = bleach!(c.mobile, masks)
            c.immobile  = bleach!(c.immobile, masks)

            # Save time evolution

        end
        
        return nothing
    end

    function bleach!(c::Array{Float32, 2}, masks::Array)
        for mask in masks
            c = c .* mask
        end

        return c
    end

    function step!(c::Array{Float32, 2}, ξ::Array{Float32, 2}, D::Float32, δt::Float32)

        dims = [1, 2]

        ĉ = fft(c, dims)
        ĉ = exp.(-D .* ξ .* δt) .* ĉ
        c = abs.(ifft(ĉ, dims))

        return c

    end


end