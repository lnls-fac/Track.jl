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
first_at_order(V, order) = (order==0) ? ui(0) : _binomials[C(V, order-1)]
last_at_order(V, order)  = (order==0) ? ui(1) : _binomials[C(V, order)]
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
    counter::ui = 0
    power ::Vector{ui} = zeros(ui, V)
    power1::Vector{ui} = zeros(ui, V)
    power2::Vector{ui} = zeros(ui, V)
    for n1::ui in 0:N
        for i1::ui in first_at_order(V,n1):last_at_order(V,n1)-1
            println("n1 = $n1, i1 = $i1")
            # get_power(i1, power1)
            # n2::ui = 0
            # while n1+n2 <= N
            #     for i2::ui in first_at_order(n2):last_at_order(n2)-1
            #         get_power(i2, power2)
            #         for i::ui in 0 +1:V-1 +1
            #             power[i] = power1[i] + power2[i]
            #             powers_[counter][i] = power[i]
            #         end
            #         osip_[counter+1] = get_index(power)
            #         counter += 1
            #     end
            #     n2 += 1
            # end
        end
    end

    return binomials_, powers_, osip_
end
V=2; N=2

_osip_val = et_osip(V, N, N)
_osip = zeros(ui, _osip_val)
_powers = zeros(ui, (_osip_val, V))
_binomials = zeros(ui, ((N+V+1)*(N+V+2))>>1)

# set binomials
for s::ui in 0:N+V
    _binomials[C(s, 0)] = _binomials[C(0, s)] = 1
    for n::ui in 1:(s-1) _binomials[C(s-n, n)] = _binomials[C((s-1)-n,n)] + _binomials[C((s-1)-(n-1), (n-1))] end
end
_powers
function get_power(idx_, pwr::Vector{ui})
    idx = ui(idx_) + 1
    while _binomials[C(V, pwr[1])] < idx
        pwr[1] += 1
    end
    idx = _binomials[C(V, pwr[1])] - idx
    for i in 2:V
        pwr[i] = pwr[i-1]
        while ui((pwr[i]*_binomials[C(V-i, pwr[i])]) / (V+1-i)) > (idx)
            pwr[i] -= 1
        end
        idx -= (pwr[i] * _binomials[C(V-i, pwr[i])]) / (V+1-i)
        pwr[i-1] -= pwr[i]
    end
end
pwr1 = zeros(ui, V)
get_power(0, pwr1)
println("$(pwr1[1]) $(pwr1[2])")
pwr1 = zeros(ui, V)
get_power(1, pwr1)
println("$(pwr1[1]) $(pwr1[2])")
pwr1 = zeros(ui, V)
get_power(2, pwr1)
println("$(pwr1[1]) $(pwr1[2])")
pwr1 = zeros(ui, V)
get_power(3, pwr1)
println("$(pwr1[1]) $(pwr1[2])")
pwr1 = zeros(ui, V)
get_power(4, pwr1)
println("$(pwr1[1]) $(pwr1[2])")
pwr1 = zeros(ui, V)
get_power(5, pwr1)
println("$(pwr1[1]) $(pwr1[2])")