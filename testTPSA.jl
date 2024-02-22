
const ui = Int

function et_osip(V::I, N::I, n::I) where I<:Union{UInt, Int}
    if n == 0
        return I(binomial(V, 0)*binomial(V+N, N))
    else
        return I((V * binomial(V+n, n) * binomial(V+N-n, N-n))/(V+n) + et_osip(V, N, n-1))
    end
end

C(v, n) = ui((((v+n-1)*(n+v))>>1) + n) + 1
function calc_tpsa_params(V, N)
    
    osip_val_ = et_osip(V, N, N)
    osip_ = zeros(ui, osip_val_)
    powers_ = zeros(ui, (osip_val_, V))
    binomials_ = zeros(ui, ((N+V+1)*(N+V+2))>>1)

    # set binomials
    for s::ui in 0:N+V
        binomials_[C(s, 0)] = binomials_[C(0, s)] = 1
        for n::ui in 1:(s-1) binomials_[C(s-n, n)] = binomials_[C((s-1)-n,n)] + binomials_[C((s-1)-(n-1), (n-1))] end
    end


     # auxiliary functions

    first_at_order(order) = (order==0) ? ui(0) : binomials_[C(V, order-1)]

    last_at_order(order)  = (order==0) ? ui(1) : binomials_[C(V, order)]
 
    function get_power(idx_)
        pwr::Vector{ui} = zeros(ui, V)
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
        return pwr
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

    # set OSIP and powers
    counter::ui = 1
    power ::Vector{ui} = zeros(ui, V)
    power1::Vector{ui} = zeros(ui, V)
    power2::Vector{ui} = zeros(ui, V)
    for n1::ui in 0:N
        for i1::ui in first_at_order(n1):last_at_order(n1)-1
            power1 = get_power(i1)
            n2::ui = 0
            while n1+n2 <= N
                for i2::ui in first_at_order(n2):last_at_order(n2)-1
                    println("counter = $counter,   n1 = $n1,   i1 = $i1,   n2 = $n2,   i2 = $i2")
                    power2 = get_power(i2)
                    println("pwr1 = $power1,   pwr2 = $power2\n")
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

x, y, z = calc_tpsa_params(3,1);

z

x

y
