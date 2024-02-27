using Base

###################### auxiliary ########################

# getindex of c
function Base.getindex(t::Tpsa{V,N,T}, index::Union{Int, UInt}) where {V,N,T}
    return t.c[index + 1]
end

# setindex of c
function Base.setindex!(t::Tpsa{V,N,T}, value::T, index::Union{Int, UInt}) where {V,N,T}
    t.c[index + 1] = value
end

# basic copy
function Base.copy(t::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = Tpsa{V,N,T}(t.c)
    return r
end

# length / size
function Base.length(t::Tpsa{V,N,T}) where {V,N,T}
    return length(t.c)
end

function Base.size(t::Tpsa{V,N,T}) where {V,N,T}
    return length(t.c)
end

########################## sum #############################

# sum constant
function Base.:+(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    r::Tpsa{V,N,T} = copy(t)
    @inbounds r[0] += T(a)
    return r
end

function Base.:+(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    r::Tpsa{V,N,T} = copy(t)
    return r + a
end

# sum "tpsa"s
function Base.:+(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = copy(t1)
    for i in 1:length(r)
        @inbounds r.c[i] += t2.c[i]
    end
    return r
end

########################## sub #############################

# subtract constant
function Base.:-(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    return t + (-a)
end

# subtract "tpsa"s
function Base.:-(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = copy(t1)
    for i in 1:length(r)
        @inbounds r.c[i] -= t2.c[i]
    end
    return r
end

# operator (-)
function Base.:-(t1::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = copy(t1)
    for i in 1:length(r)
        @inbounds r.c[i] = -(r.c[i])
    end
    return r
end

# subtract "tpsa"
function Base.:-(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return a + (-t)
end

########################## mult ############################

# multiply constant
function Base.:*(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    r::Tpsa{V,N,T} = copy(t)
    for i in 1:length(r)
        @inbounds r.c[i] *= T(a)
    end
    return r
end

function Base.:*(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return t * a
end

# multiply "tpsa"s
function Base.:*(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = Tpsa{V,N,T}()
    counter::ui = 1
    for n1::ui in 0:N
        for i1::ui in tpsa_global_params[V,N].first_at_order(n1):tpsa_global_params[V,N].last_at_order(n1)-1
            n2::ui = 0
            while n1+n2 <= N
                for i2::ui in tpsa_global_params[V,N].first_at_order(n2):tpsa_global_params[V,N].last_at_order(n2)-1
                    @inbounds r[tpsa_global_params[V,N].osip[counter]] += t1[i1] * t2[i2]
                    counter += 1
                end
                n2 += 1
            end
        end
    end
    return r
end

########################## inv ############################

# invert tpsa
function inverse(t::Tpsa{V,N,T}) where {V, N, T}
    r::Tpsa{V,N,T} = Tpsa{V,N,T}()
    @inbounds a::T = t[0]
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    x /= a
    p::Tpsa{V,N,T} = Tpsa{V,N,T}(1, 0)
    for i in 0:N
        r += p * (Bool(i&1) ? -1 : 1)
        p *= x
    end
    return r / a
end

########################## div ############################

# divide constant
function Base.:/(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    return t * (1/a)
end

function Base.:/(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return a * inverse(t)
end

# divide "tpsa"s
function Base.:/(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    return t1 * inverse(t2)
end

###################### comparators ########################

# equal constant
function Base.:(==)(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    @inbounds return (t[0] == T(a))
end

function Base.:(==)(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    @inbounds return (t[0] == T(a))
end

# equal "tpsa"s
function Base.:(==)(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    return ((t1 - t2) == T(0))
end

# not equal -> not necessary: the function "==" do the job

# is less constant
function Base.:(<)(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    @inbounds return (T(a) < t[0])
end

function Base.:(<)(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    @inbounds return (t[0] < T(a))
end

# is less "tpsa"
function Base.:(<)(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T}
    return ((t1 - t2) < T(0))
end

# is greater -> not necessary: the function "<" do the job

#################### math functions ######################

# pow
function Base.:(^)(t::Tpsa{V,N,T}, a::Union{Int,UInt}) where {V, N, T}
    r::Tpsa{V,N,T} = copy(t)
    for _ in 1:a-1
        r *= t
    end
    return r
end

# abs
function Base.abs(t::Tpsa{V,N,T}) where {V, N, T}
    if t >= T(0)
        return t
    else
        return -t
    end
end

# sqrt
function Base.sqrt(t::Tpsa{V,N,T}) where {V, N, T}
    r::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    @inbounds x /= t[0]
    p::Tpsa{V,N,T} = Tpsa{V, N, T}(1.0, 0)
    f::T = T(1)
    for i in 0:N
        r += p * f
        f *= (0.5 - i)/(i+1)
        p *= x
    end
    @inbounds r *= sqrt(t[0])
    return r
end

# log
function Base.log(t::Tpsa{V,N,T}) where {V, N, T}
    r::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    @inbounds x /= t[0]
    p::Tpsa{V,N,T} = copy(x)
    for i in 1:N
        r += p / ((Bool(i&1) ? 1.0 : -1.0) * i)
        p *= x
    end
    @inbounds r += log(t[0])
    return r
end

# cos
function Base.cos(t::Tpsa{V,N,T}) where {V, N, T}
    rc::Tpsa{V,N,T} = Tpsa{V, N, T}(1.0, 0)
    rs::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    p::Tpsa{V,N,T} = copy(x)
    fac::ui = 1
    for i in 1:N
        if Bool(i&1)
            rs += (Bool(i&2) ? -1 : 1) * p / fac
        else
            rc += (Bool(i&2) ? -1 : 1) * p / fac
        end
        p *= x
        fac *= i + 1
    end
    @inbounds return (cos(t[0]) * rc) - (sin(t[0]) * rs)
end

# sin
function Base.sin(t::Tpsa{V,N,T}) where {V, N, T}
    rc::Tpsa{V,N,T} = Tpsa{V, N, T}(1.0, 0)
    rs::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    p::Tpsa{V,N,T} = copy(x)
    fac::ui = 1
    for i in 1:N
        if Bool(i&1)
            rs += (Bool(i&2) ? -1 : 1) * p / fac
        else
            rc += (Bool(i&2) ? -1 : 1) * p / fac
        end
        p *= x
        fac *= i + 1
    end
    @inbounds return (sin(t[0]) * rc) + (cos(t[0]) * rs)
end

# tan
function Base.tan(t::Tpsa{V,N,T}) where {V, N, T}
    return sin(t) / cos(t)
end

# cosh
function Base.cosh(t::Tpsa{V,N,T}) where {V, N, T}
    rc::Tpsa{V,N,T} = Tpsa{V, N, T}(1.0, 0)
    rs::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    p::Tpsa{V,N,T} = copy(x)
    fac::ui = 1
    for i in 1:N
        if Bool(i&1)
            rs += p / fac
        else
            rc += p / fac
        end
        p *= x
        fac *= i + 1
    end
    @inbounds return (cosh(t[0]) * rc) + (sinh(t[0]) * rs)
end

# sinh
function Base.sinh(t::Tpsa{V,N,T}) where {V, N, T}
    rc::Tpsa{V,N,T} = Tpsa{V, N, T}(1.0, 0)
    rs::Tpsa{V,N,T} = Tpsa{V, N, T}()
    x::Tpsa{V,N,T} = copy(t)
    @inbounds x[0] = T(0)
    p::Tpsa{V,N,T} = copy(x)
    fac::ui = 1
    for i in 1:N
        if Bool(i&1)
            rs += p / fac
        else
            rc += p / fac
        end
        p *= x
        fac *= i + 1
    end
    @inbounds return (sinh(t[0]) * rc) - (cosh(t[0]) * rs)
end

# tanh
function Base.tanh(t::Tpsa{V,N,T}) where {V, N, T}
    return sinh(t) / cosh(t)
end

# isinfinite
function Base.isfinite(t::Tpsa{V,N,T}) where {V, N, T}
    @inbounds return isfinite(t[0])
end