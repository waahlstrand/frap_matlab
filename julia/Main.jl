
# Includes
include("./Utils.jl")
include("./Bleach.jl")

# Packages
using .Utils: Concentration
using .Bleach: evolve!, create_fourier_grid, create_imaging_bleach_mask

function main()


    c₀::Float32 = 0.1
    ϕₘ::Float32 = 1.0
    D::Float32  = 0.1
    ξ = Array{Float32, 2}(undef, 2, 2)
    δt::Float32 = 0.1
    masks = nothing

    α = 0.1
    β = 1.0
    γ = 0.0

    const n_pixels     = 256
    const n_pad_pixels = 128
    const n_prebleach_frames    = 10
    const n_bleach_frames       = 1
    const n_postbleach_frames   = 100
    const n_frames              = n_prebleach_frames + n_postbleach_frames 

    dims         = (n_pixels + 2 * n_pad_pixels, n_pixels + 2 * n_pad_pixels, n_frames)

    prebleach   = 1:n_prebleach_frames
    bleach      = n_prebleach_frames:n_prebleach_frames+n_bleach_frames
    postbleach  = n_prebleach_frames+n_bleach_frames:n_frames

    C = Concentration(c₀, ϕₘ, dims)
    
    imaging_mask = create_imaging_bleach_mask(β, n_pixels, n_pad_pixels)

    stages = (prebleach, bleach, postbleach)
    

    for stage in [prebleach, bleach, postbleach]

        C = evolve!(C, ξ, D, masks, stage)

    end

    #C = C |> evolve!(ξ, D, masks, prebleach) |> evolve!(ξ, D, masks, bleach) |> evolve!(ξ, D, masks, postbleach)

end

# Run main script
main()