include("TPSA_aux.jl")

mutable struct Tpsa{V, N, T}
    c::Vector{T}

    function Tpsa{V, N, T}(a::X, v::Union{Int, UInt}) where {V, N, T, X<:Number}
        _c = zeros(T, binomial(V+N, N))
        _c[1] = T(a)
        if (1<=v<=V) 
            _c[v+1] = T(1) 
        end
        if !haskey(tpsa_global_params, (V,N))
            init_tpsa_params(V, N)
        end
        new{V, N, T}(_c)
    end

    function Tpsa{V, N}(a::X, v::Union{Int, UInt}) where {V, N, X <: Number}
        return Tpsa{V, N, Float64}(a, v)
    end

    function Tpsa{V, N, T}() where {V,N,T}
        return Tpsa{V,N,T}(T(0.0), 0)
    end

    function Tpsa{V, N}() where {V, N}
        return Tpsa{V, N, Float64}(0.0, 0)
    end

    function Tpsa{V, N, T}(vec::Union{Array{X}, Vector{X}}) where {V, N, T, X<:Number}
        _c = zeros(T, binomial(V+N, N))
        if length(vec) != length(_c)
            error("Invalid vector")
        end
        _c .= T.(vec)
        if !haskey(tpsa_global_params, (V,N))
            init_tpsa_params(V, N)
        end
        new{V, N, T}(_c)
    end

    function Tpsa{V, N}(vec::Union{Array{X}, Vector{X}}) where {V, N, X<:Number}
        return Tpsa{V, N, Float64}(vec)
    end

end

function Base.show(io::IO, t::Tpsa{V, N, T}) where {V, N, T}
    println(io, "Tpsa{V, N}: $(t.c)")
end

function Base.convert(::Type{Tpsa{V, N, T}}, x::R) where {V, N, T, R<:Number}
    return Tpsa{V, N, T}(T(x), 0) 
end
function Base.convert(::Type{Tpsa{V, N}}, x::R) where {V, N, R<:Number} 
    return Tpsa{V, N, Float64}(Float64(x), 0) 
end

include("TPSA_operators.jl")