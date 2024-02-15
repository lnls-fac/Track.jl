# Pos.jl

using Base
using Printf

export Pos

mutable struct Pos{Float64}
    rx::Float64
    px::Float64
    ry::Float64
    py::Float64
    de::Float64
    dl::Float64
    function Pos(rx::T, px::T, ry::T, py::T, de::T, dl::T) where T
        new{Float64}(Float64(rx), Float64(px), Float64(ry), Float64(py), Float64(de), Float64(dl))
    end
    function Pos(v::Vector{T}) where T
        new{Float64}(Float64(v[1]), Float64(v[2]), Float64(v[3]), Float64(v[4]), Float64(v[5]), Float64(v[6]))
    end
    function Pos(f::T) where T
        new{Float64}(Float64(f), Float64(f), Float64(f), Float64(f), Float64(f), Float64(f))
    end
end

function Base.:+(v1::Pos{T}, v2::Pos{T}) where T
    return Pos(v1.rx + v2.rx, v1.px + v2.px, v1.ry + v2.ry, v1.py + v2.py, v1.de + v2.de, v1.dl + v2.dl)
end

function Base.:+(m1::Vector{Pos{T}}, m2::Vector{Pos{T}}) where T
    return [m1[i] + m2[i] for i in 1:length(m1)]
end

function Base.:-(v1::Pos{T}, v2::Pos{T}) where T
    return Pos(v1.rx - v2.rx, v1.px - v2.px, v1.ry - v2.ry, v1.py - v2.py, v1.de - v2.de, v1.dl - v2.dl)
end

function Base.:-(m1::Vector{Pos{T}}, m2::Vector{Pos{T}}) where T
    return [m1[i] - m2[i] for i in 1:length(m1)]
end

function  Base.:*(v::Pos{T}, scalar::S) where {T, S<:Real}
    return Pos(v.rx * scalar, v.px * scalar, v.ry * scalar, v.py * scalar, v.de * scalar, v.dl * scalar)
end

function  Base.:*(scalar::S, v::Pos{T}) where {T, S<:Real}
    return Pos(v.rx * scalar, v.px * scalar, v.ry * scalar, v.py * scalar, v.de * scalar, v.dl * scalar)
end

function  Base.:*(scalar::S, m1::Vector{Pos{T}}) where {T, S<:Real}
    return [scalar * v for v in m1]
end

function  Base.:/(v::Pos{T}, scalar::S) where {T, S<:Real}
    return (1 / scalar) * v
end

function Base.copy(p::Pos{T}) where T
    return Pos(p.rx, p.px, p.ry, p.py, p.de, p.dl)
end

function Base.show(io::IO, ::MIME"text/plain", p::Pos{T}) where T
    println(io, "[", p.rx, "\n ", p.px, "\n ", p.ry, "\n ", p.py, "\n ", p.de, "\n ", p.dl, "]")
end

function Base.show(io::IO, p::Pos{T}) where T
    println(io, "[", p.rx, "\n ", p.px, "\n ", p.ry, "\n ", p.py, "\n ", p.de, "\n ", p.dl, "]")
end

function Base.show(io::IO, ::MIME"text/plain", p::Vector{Pos{T}}) where T
    println(io, join(["[$(@sprintf("%1.8e  ", v.rx))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.px))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.ry))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.py))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.de))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e] ", v.dl))  " for v in p]))
end

function Base.show(io::IO, p::Vector{Pos{T}}) where T
    println(io, join(["[$(@sprintf("%1.8e  ", v.rx))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.px))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.ry))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.py))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e  ", v.de))  " for v in p]))
    println(io, join([" $(@sprintf("%1.8e] ", v.dl))  " for v in p]))
end

function Base.getindex(p::Pos{T}, index::Int) where T
    if index == 1
        return p.rx
    elseif index == 2
        return p.px
    elseif index == 3
        return p.ry
    elseif index == 4
        return p.py
    elseif index == 5
        return p.de
    elseif index == 6
        return p.dl
    else
        error("invaid Pos index: $index, should stay between 1 and 6")
    end
end

function Base.getindex(p::Pos{T}, ::Colon) where T
   return Float64[p.rx, p.px, p.ry, p.py, p.de, p.dl]
end

function Base.setindex!(p::Pos{T}, value::Real, index::Int) where T
    if index == 1
        p.rx = Float64(value)
    elseif index == 2
        p.px = Float64(value)
    elseif index == 3
        p.ry = Float64(value)
    elseif index == 4
        p.py = Float64(value)
    elseif index == 5
        p.de = Float64(value)
    elseif index == 6
        p.dl = Float64(value)
    else
        error("invaid Pos index: $index, should stay between 1 and 6")
    end
end

function Pos_get_max(v::Pos{T}) where T
    max_val = abs(v.rx)
    max_val = max(abs(v.px), max_val)
    max_val = max(abs(v.ry), max_val)
    max_val = max(abs(v.py), max_val)
    max_val = max(abs(v.de), max_val)
    max_val = max(abs(v.dl), max_val)
    return max_val
end
