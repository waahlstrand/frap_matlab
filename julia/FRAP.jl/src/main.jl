using FFTW: fft, ifft, ifftshift

function evolve!(C::Concentration{Float64},
                 ξ²::Array{Float64, 2},
                 D::Float64,
                 δt::Float64,
                 masks::Array{Array{Float64,2},1},
                 frames::UnitRange{Int64})

    # Use simple Markov memory, one timestep
    c = C[:, :, frames[begin]]

    for frame in frames[begin+1:end]

        # Make one step of length δt
        c.mobile .= step!(c.mobile, ξ², D, δt)

        # Apply imaging bleach masks
        c.mobile    .= bleach!(c.mobile, masks)
        c.immobile  .= bleach!(c.immobile, masks)

        # Save time evolution
        C[:, :, frame] = c

    end

    return C
end

function evolve(C::Concentration{Float64},
                ξ²::Array{Float64, 2},
                D::Float64,
                δt::Float64,
                masks::Array{Array{Float64,2},1},
                n_frames::Int64)

    # Use simple Markov memory, one timestep
    c = C[:, :]
    C_to_save = similar(C, (size(C)..., n_frames))

    for frame in 1:n_frames

        # Make one step of length δt
        c.mobile = step!(c.mobile, ξ², D, δt)

        # Apply imaging bleach masks
        # c.mobile    = bleach!(c.mobile, masks)
        # c.immobile  = bleach!(c.immobile, masks)

        # Save time evolution
        C_to_save[:, :, frame] = c

    end

    return C_to_save
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
