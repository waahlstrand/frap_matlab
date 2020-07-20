mutable struct Concentration{T<:Real}
    mobile::Array{T}
    immobile::Array{T}

    function Concentration(c₀::T, ϕₘ::T, dims::Tuple{Int,Vararg{Int}}) where {T}

        mobile = c₀ * ϕₘ * ones(T, dims)
        immobile = c₀ * (1 - ϕₘ) * ones(T, dims)

        return new{T}(mobile, immobile)
    end

    function Concentration(mobile::Array{T}, immobile::Array{T}) where {T}
        return new{T}(mobile, immobile)
    end

    function Concentration(mobile::T, immobile::T) where {T}
        return new{T}([mobile], [immobile])
    end
end

function Base.getindex(C::Concentration, inds...)
    return Concentration(C.mobile[inds...], C.immobile[inds...])
end

function Base.setindex!(C::Concentration, v::Concentration, inds...)

    C.mobile[inds...] = v.mobile
    C.immobile[inds...] = v.immobile

    return C
end

function Base.size(C::Concentration)

    return size(C.mobile)

end

function Base.lastindex(C::Concentration, d::Int64)

    return lastindex(C.mobile, d)

end

function Base.firstindex(C::Concentration, d::Int64)

    return firstindex(C.mobile, d)

end

function Base.similar(C::Concentration)

    mobile = similar(C.mobile)
    immobile = similar(C.immobile)

    return Concentration(mobile, immobile)

end

function Base.similar(C::Concentration, dims::Tuple{Int,Vararg{Int}})

    mobile = similar(C.mobile, dims)
    immobile = similar(C.immobile, dims)

    return Concentration(mobile, immobile)

end
