
module FRAP
    include("concentration.jl")
    include("bleaching.jl")
    include("main.jl")

    export Concentration, ROI
    export evolve, evolve!, create_fourier_grid, create_imaging_bleach_mask

end
