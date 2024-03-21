const ui = Int

function _calc_osip(V::I, N::I, n::I) where I<:ui
    if n == 0
        return I(binomial(V, 0)*binomial(V+N, N))
    else
        return I((V * binomial(V+n, n) * binomial(V+N-n, N-n))/(V+n) + _calc_osip(V, N, n-1))
    end
end

_C(v, n) = ui((((v+n-1)*(n+v))>>1) + n) + 1

function _calc_tpsa_params(V, N)
    V = ui(V)
    N = ui(N)
    osip_val_ = _calc_osip(V, N, N)
    osip_ = zeros(ui, osip_val_)
    powers_ = zeros(ui, (osip_val_, V))
    binomials_ = zeros(ui, ((N+V+1)*(N+V+2))>>1)

    # set binomials
    for s::ui in 0:N+V
        binomials_[_C(s, 0)] = binomials_[_C(0, s)] = 1
        for n::ui in 1:(s-1) binomials_[_C(s-n, n)] = binomials_[_C((s-1)-n,n)] + binomials_[_C((s-1)-(n-1), (n-1))] end
    end

     # auxiliary functions

    _first_at_order(order) = (order==0) ? ui(0) : binomials_[_C(V, order-1)]

    _last_at_order(order)  = (order==0) ? ui(1) : binomials_[_C(V, order)]
 
    function get_power(idx_)
        pwr::Vector{ui} = zeros(ui, V)
        idx = ui(idx_) + 1
        while binomials_[_C(V, pwr[1])] < idx
            pwr[1] += 1
        end
        idx = binomials_[_C(V, pwr[1])] - idx
        for i in 2:V
            pwr[i] = pwr[i-1]
            while ui((pwr[i]*binomials_[_C(V-i, pwr[i])]) / (V+1-i)) > (idx)
                pwr[i] -= 1
            end
            idx -= (pwr[i] * binomials_[_C(V-i, pwr[i])]) / (V+1-i)
            pwr[i-1] -= pwr[i]
        end
        return pwr
    end
 
    function get_index(power_)
        s_i::ui = 0
        idx::ui = 0
        for i::ui in 0:V-2
            s_i += power_[V-i]
            idx += (s_i * binomials_[_C(i,s_i)]) / (i+1)
        end
        s_i += power_[1]
        idx = binomials_[_C(V, s_i)] - idx
        return idx - 1
    end

    # set OSIP and powers
    counter::ui = 1
    power ::Vector{ui} = zeros(ui, V)
    power1::Vector{ui} = zeros(ui, V)
    power2::Vector{ui} = zeros(ui, V)
    for n1::ui in 0:N
        for i1::ui in _first_at_order(n1):_last_at_order(n1)-1
            power1 = get_power(i1)
            n2::ui = 0
            while n1+n2 <= N
                for i2::ui in _first_at_order(n2):_last_at_order(n2)-1
                    power2 = get_power(i2)
                    for i::ui in 1:V
                        power[i] = power1[i] + power2[i]
                        powers_[counter, i] = power[i]
                    end
                    osip_[counter] = get_index(power)
                    counter += 1
                end
                n2 += 1
            end
        end
    end
    return binomials_, powers_, osip_
end

struct TpsaParams
    # base vectors
    binomials::Vector{ui}
    powers::Matrix{ui}
    osip::Vector{ui}
    # functions
    first_at_order ::Function
    last_at_order  ::Function
    get_power      ::Function
    get_index      ::Function
    function TpsaParams(V, N)
        b_, p_, o_ = _calc_tpsa_params(V, N)
        first_at_order::Function = order -> (order==0) ? ui(0) : b_[_C(V, order-1)]  
        last_at_order ::Function = order -> (order==0) ? ui(1) : b_[_C(V, order)] 
        get_power     ::Function = idx_ -> begin
            pwr::Vector{ui} = zeros(ui, V)
            idx = ui(idx_) + 1
            while b_[_C(V, pwr[1])] < idx
                pwr[1] += 1
            end
            idx = b_[_C(V, pwr[1])] - idx
            for i in 2:V
                pwr[i] = pwr[i-1]
                while ui((pwr[i]*b_[_C(V-i, pwr[i])]) / (V+1-i)) > (idx)
                    pwr[i] -= 1
                end
                idx -= (pwr[i] * b_[_C(V-i, pwr[i])]) / (V+1-i)
                pwr[i-1] -= pwr[i]
            end
            return pwr
        end
        get_index::Function = power_ -> begin
            s_i::ui = 0
            idx::ui = 0
            for i::ui in 0:V-2
                s_i += power_[V-i]
                idx += (s_i * b_[_C(i,s_i)]) / (i+1)
            end
            s_i += power_[1]
            idx = b_[_C(V, s_i)] - idx
            return idx - 1
        end 
        return new(b_, p_, o_, first_at_order, last_at_order, get_power, get_index)
    end
end

global tpsa_global_params::Dict{Tuple, TpsaParams} = Dict()

function init_tpsa_params(V, N)
    push!(tpsa_global_params, (V, N) => TpsaParams(V, N))
end

# Base.@kwdef struct TPSA_Aux_Funcs
#     v::Union{Int, UInt}
#     n::Union{Int, UInt}
#     first_at_order::Function = order -> (order==0) ? ui(0) : tpsa_global_params[v,n].binomials[_C(v, order-1)]  
#     last_at_order ::Function = order -> (order==0) ? ui(1) : tpsa_global_params[v,n].binomials[_C(v, order)]  
#     get_power     ::Function = idx_ -> begin
#         pwr::Vector{ui} = zeros(ui, v)
#         idx = ui(idx_) + 1
#         while tpsa_global_params[v,n].binomials[_C(v, pwr[1])] < idx
#             pwr[1] += 1
#         end
#         idx = tpsa_global_params[v,n].binomials[_C(v, pwr[1])] - idx
#         for i in 2:v
#             pwr[i] = pwr[i-1]
#             while ui((pwr[i]*tpsa_global_params[v,n].binomials[_C(v-i, pwr[i])]) / (v+1-i)) > (idx)
#                 pwr[i] -= 1
#             end
#             idx -= (pwr[i] * tpsa_global_params[v,n].binomials[_C(v-i, pwr[i])]) / (v+1-i)
#             pwr[i-1] -= pwr[i]
#         end
#         return pwr
#     end
#     get_index::Function = power_ -> begin
#         s_i::ui = 0
#         idx::ui = 0
#         for i::ui in 0:v-2
#             s_i += power_[v-i]
#             idx += (s_i * tpsa_global_params[v,n].binomials[_C(i,s_i)]) / (i+1)
#         end
#         s_i += power_[1]
#         idx = tpsa_global_params[v,n].binomials[_C(v, s_i)] - idx
#         return idx - 1
#     end
# end