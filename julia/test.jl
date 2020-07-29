# if pwd() in LOAD_PATH
#     println("Local directory in path.")
# else
#     push!(LOAD_PATH, pwd())
#     println("Local directory added to path.")
# end

include("FRAP.jl/src/mask.jl")
include("FRAP.jl/src/concentration.jl")
include("FRAP.jl/src/main.jl")

n_pixels     = 256
n_pad_pixels = 128
n_prebleach_frames    = 10
n_bleach_frames       = 1
n_postbleach_frames   = 100
n_frames              = n_prebleach_frames + n_postbleach_frames
n_elements            = n_pixels + 2*n_pad_pixels

dims = (n_elements, n_elements, n_frames)

c₀ = 1.0
ϕₘ = 1.0
D  = 0.1
ξ² = create_fourier_grid(n_elements)
δt = 0.1
β = 0.9
mask = create_imaging_bleach_mask(β, n_pixels, n_pad_pixels)
masks = [mask]

C = Concentration(c₀, ϕₘ, dims)
C = evolve!(C, ξ², D, δt, masks, 1:n_prebleach_frames)


@gif for i in 1:10
    heatmap(C.mobile[:,:,i], clim=(0,1))
end
