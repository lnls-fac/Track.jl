using Base
const ui = UInt32
function et_osip(V::I, N::I, n::I) where I<:Union{UInt, Int}
    if n == 0
        return I(binomial(V, 0)*binomial(V+N, N))
    else
        return I((V * binomial(V+n, n) * binomial(V+N-n, N-n))/(V+n) + et_osip(V, N, n-1))
    end
end
struct TpsaParams
    binomials::Vector{ui}
    powers::Matrix{ui}
    osip::Vector{ui}
    function TpsaParams(V,N)
        return new(calc_tpsa_params(V,N))
    end
end
global tpsa_global_params::Dict{Tuple{ui}, TpsaParams} = Dict()
function init_tpsa_params(V, N)
    push!(tpsa_global_params, TpsaParams(V, N))
end
mutable struct Tpsa{V, N, T}
    c::Vector{T}

    function Tpsa{V, N, T}(a::X, v::Union{Int, UInt}) where {V, N, T, X}
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

    function Tpsa{V, N}(a::X, v::Union{Int, UInt}) where {V, N, X}
        return Tpsa{V, N, Float64}(a, v)
    end

    function Tpsa{V, N, T}() where {V,N,T}
            return Tpsa{V,N,T}(0.0, 0)
    end

    function Tpsa{V, N}() where {V, N}
        return Tpsa{V, N, Float64}(0.0, 0)
    end

end

C(v, n) = ui((((v+n-1)*(n+v))>>1) + n) + 1
function calc_tpsa_params(V, N)
    
    osip_val_ = et_osip(V, N, N)
    osip_ = zeros(ui, _osip_val)
    powers_ = zeros(ui, (_osip_val, V))
    binomials_ = zeros(ui, ((N+V+1)*(N+V+2))>>1)

    # set binomials
    for s::ui in 0:N+V
        binomials_[C(s, 0)] = binomials_[C(0, s)] = 1
        for n::ui in 1:(s-1) binomials_[C(s-n, n)] = binomials_[C((s-1)-n,n)] + binomials_[C((s-1)-(n-1), (n-1))] end
    end

    # set OSIP and powers
    counter::ui = 1
    power ::Vector{ui} = zeros(ui, V)
    power1::Vector{ui} = zeros(ui, V)
    power2::Vector{ui} = zeros(ui, V)
    for n1::ui in 0:N
        for i1::ui in first_at_order(n1):last_at_order(n1)-1
            get_power(i1, power1)
            n2::ui = 0
            while n1+n2 <= N
                for i2::ui in first_at_order(n2):last_at_order(n2)-1
                    get_power(i2, power2)
                    for i::ui in 1:V
                        power[i] = power1[i] + power2[i]
                        powers_[counter][i] = power[i]
                    end
                    osip_[counter] = get_index(power)
                    counter += 1
                end
                n2 += 1
            end
        end
    end

    # auxiliary functions

    first_at_order(order) = (order==0) ? ui(0) : binomials_[C(V, order-1)]

    last_at_order(order)  = (order==0) ? ui(1) : binomials_[C(V, order)]

    function get_power(idx_, pwr::Vector{ui})
        idx = ui(idx_) + 1
        while binomials_[C(V, pwr[1])] < idx
            pwr[1] += 1
        end
        idx = binomials_[C(V, pwr[1])] - idx
        for i in 2:V
            pwr[i] = pwr[i-1]
            while ui((pwr[i]*binomials_[C(V-i, pwr[i])]) / (V+1-i)) > (idx)
                pwr[i] -= 1
            end
            idx -= (pwr[i] * binomials_[C(V-i, pwr[i])]) / (V+1-i)
            pwr[i-1] -= pwr[i]
        end
    end

    function get_index(power_)
        s_i::ui = 0
        idx::ui = 0
        for i::ui in 0:V-2
            s_i += power_[V-i]
            idx += (s_i * binomials_[C(i,s_i)]) / (i+1)
        end
        s_i += power_[1]
        idx = binomials_[C(V, s_i)] - idx
        return idx - 1
    end

    return binomials_, powers_, osip_

end