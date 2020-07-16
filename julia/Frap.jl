
# Includes
include("./Utils.jl")
include("./Bleach.jl")

# Packages
using .Utils: Concentration
using .Bleach: evolve

function main()


    c₀::Float32 = 0.1
    ϕₘ::Float32 = 1.0
    D::Float32  = 0.1
    ξ = Array{Float32, 2}(undef, 2, 2)
    δt::Float32 = 0.1
    masks = nothing


    const n_pixels     = 256
    const n_pad_pixels = 128
    const n_prebleach_frames    = 10
    const n_bleach_frames       = 1
    const n_postbleach_frames   = 100
    const n_frames              = n_prebleach_frames + n_postbleach_frames 

    dims         = (n_pixels + 2 * n_pad_pixels, n_pixels + 2 * n_pad_pixels, n_frames)

    init        = [:, :, 1]
    prebleach   = [:, :, 1:n_prebleach_frames]
    bleach      = [:, :, n_prebleach_frames:n_prebleach_frames+n_bleach_frames]
    postbleach  = [:, :, n_prebleach_frames+n_bleach_frames:end]

    C = Concentration(c₀, ϕₘ, dims)

    # C  = Concentration(Cₜ.mobile[:, :, 1], 
    #                    Cₜ.immobile[:, :, 1])



    # C = Concentration(Cₜ.mobile[:, :, n_prebleach_frames], 
    #                   Cₜ.immobile[:, :, n_prebleach_frames])

    C[prebleach...]    = evolve(C[init...],       ξ, D, masks, n_prebleach_frames)
    C[bleach...]       = evolve(C[prebleach...],  ξ, D, masks, n_bleach_frames)
    C[postbleach...]   = evolve(C[bleach...],     ξ, D, masks, n_bleach_frames)

end

# Run main script
main()