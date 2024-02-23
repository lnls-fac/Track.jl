using Base

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
    r[0] += T(a)
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
        r.c[i] += t2.c[i]
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
        r.c[i] -= t2.c[i]
    end
    return r
end

# operator (-)
function Base.:-(t1::Tpsa{V,N,T}) where {V,N,T}
    r::Tpsa{V,N,T} = copy(t1)
    for i in 1:length(r)
        r.c[i] = -(r.c[i])
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
        r.c[i] *= T(a)
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
        for i1::ui in r.first_at_order(n1):r.last_at_order(n1)-1
            n2::ui = 0
            while n1+n2 <= N
                for i2::ui in r.first_at_order(n2):r.last_at_order(n2)-1
                    r[tpsa_global_params[V,N].osip[counter]] += t1[i1] * t2[i2]
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
    a::T = t[0]
    x::Tpsa{V,N,T} = copy(t)
    x[0] = T(0)
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
    return (t[0] == T(a))
end

function Base.:(==)(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return (t[0] == T(a))
end

# equal "tpsa"s
function Base.:(==)(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return ((t1 - t2) == T(0))
end

# is less constant
function Base.:(<)(a::R, t::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return (T(a) < t[0])
end

function Base.:(<)(t::Tpsa{V,N,T}, a::R) where {V,N,T,R<:Real}
    return (t[0] < T(a))
end

# is less "tpsa"
function Base.:(<)(t1::Tpsa{V,N,T}, t2::Tpsa{V,N,T}) where {V,N,T,R<:Real}
    return ((t1 - t2) < T(0))
end

# is greater constant -> not necessary: the function "<" do the job