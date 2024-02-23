# Pos.jl

using Base
using Printf
using ..TPSA: Tpsa

export Pos

mutable struct Pos{T<:Union{Float64, Tpsa{6, 1, Float64}}}
    rx::T
    px::T
    ry::T
    py::T
    de::T
    dl::T
    function Pos(rx::T, px::T, ry::T, py::T, de::T, dl::T; tpsa::Bool=false) where T#<: Number
        if !tpsa
            new{Float64}(Float64(rx), Float64(px), Float64(ry), Float64(py), Float64(de), Float64(dl))
        else
            new{Tpsa{6, 1, Float64}}(
                Tpsa{6, 1, Float64}(rx, 1),
                Tpsa{6, 1, Float64}(px, 2),
                Tpsa{6, 1, Float64}(ry, 3),
                Tpsa{6, 1, Float64}(py, 4),
                Tpsa{6, 1, Float64}(de, 5),
                Tpsa{6, 1, Float64}(dl, 6),
            )
        end
    end
    function Pos(v::Vector{T}; tpsa::Bool=false) where T#<:Number
        if !(0<length(v)<=6)
            error("Vector agr must be 6-element")
        end
        if !tpsa
            new{Float64}(Float64(v[1]), Float64(v[2]), Float64(v[3]), Float64(v[4]), Float64(v[5]), Float64(v[6]))
        else
            new{Tpsa{6, 1, Float64}}(
                Tpsa{6, 1, Float64}(v[1], 1),
                Tpsa{6, 1, Float64}(v[2], 2),
                Tpsa{6, 1, Float64}(v[3], 3),
                Tpsa{6, 1, Float64}(v[4], 4),
                Tpsa{6, 1, Float64}(v[5], 5),
                Tpsa{6, 1, Float64}(v[6], 6),
            )
        end
    end
    function Pos(f::T; tpsa::Bool=false) where T#<:Number
        if !tpsa
            new{Float64}(Float64(f), Float64(f), Float64(f), Float64(f), Float64(f), Float64(f))
        else
            new{Tpsa{6, 1, Float64}}(
                Tpsa{6, 1, Float64}(f, 1),
                Tpsa{6, 1, Float64}(f, 2),
                Tpsa{6, 1, Float64}(f, 3),
                Tpsa{6, 1, Float64}(f, 4),
                Tpsa{6, 1, Float64}(f, 5),
                Tpsa{6, 1, Float64}(f, 6),
            )
        end
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

function  Base.:*(v::Pos{T}, scalar::T) where T
    return Pos(v.rx * scalar, v.px * scalar, v.ry * scalar, v.py * scalar, v.de * scalar, v.dl * scalar)
end

function  Base.:*(scalar::T, v::Pos{T}) where T
    return Pos(v.rx * scalar, v.px * scalar, v.ry * scalar, v.py * scalar, v.de * scalar, v.dl * scalar)
end

function  Base.:*(scalar::T, m1::Vector{Pos{T}}) where T
    return Vector{Pos{T}}[scalar * v for v in m1]
end

function  Base.:/(v::Pos{T}, scalar::T) where T
    return (1 / scalar) * v
end

function Base.copy(p::Pos{T}) where T
    tpsa=false
    if isa(p.rx, Tpsa)
        tpsa=true
    end
    p_new = Pos(0.0, tpsa=tpsa)
    p_new.rx = p.rx
    p_new.px = p.px
    p_new.ry = p.ry
    p_new.py = p.py
    p_new.de = p.de
    p_new.dl = p.dl
    return p_new
end

function Base.show(io::IO, p::Pos{T}) where T
    println(io, "$(typeof(p)):\n[", p.rx, "\n ", p.px, "\n ", p.ry, "\n ", p.py, "\n ", p.de, "\n ", p.dl, "]")
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
   return T[p.rx, p.px, p.ry, p.py, p.de, p.dl]
end

function Base.setindex!(p::Pos{T}, value::T, index::Int) where T
    if !(0 < index <= 6)
        error("invaid Pos index: $index, should stay between 1 and 6")
    end
    if index == 1
        p.rx = value
    elseif index == 2
        p.px = value
    elseif index == 3
        p.ry = value
    elseif index == 4
        p.py = value
    elseif index == 5
        p.de = value
    elseif index == 6
        p.dl = value
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
