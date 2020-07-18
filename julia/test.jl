include("./julia/Utils.jl")
include("./julia/Bleach.jl")

using Base
using .Utils: Concentration
using .Bleach: create_fourier_grid

n_pixels     = 256
n_pad_pixels = 128
n_prebleach_frames    = 10
n_bleach_frames       = 1
n_postbleach_frames   = 100
n_frames              = n_prebleach_frames + n_postbleach_frames 

dims         = (n_pixels + 2 * n_pad_pixels, n_pixels + 2 * n_pad_pixels, n_frames)

c₀ = 0.1
ϕₘ = 1.0
D  = 0.1
ξ = Array{Float32, 2}(undef, 2, 2)
δt = 0.1
masks = nothing

C2 = Concentration(c₀, ϕₘ, dims)
C1 = Concentration(2*c₀, ϕₘ, dims)

X = create_fourier_grid(20)