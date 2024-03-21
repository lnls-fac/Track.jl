using ..AcceleratorModule: Accelerator
using ..PosModule: Pos
using ..Tracking: element_pass, line_pass
using ..Auxiliary
using ..TPSA: Tpsa
using ..Orbit: find_orbit6


function find_matrix66(accelerator::Accelerator, fixed_point::Union{Pos{T},Nothing}=nothing; indices::String="m66", method::String="tracking") where {T}
    if lowercase(method) == "tracking"
        return _tracking_find_matrix66(accelerator, fixed_point, indices=indices)
    elseif lowercase(method) == "tpsa"
        return _tpsa_find_matrix66(accelerator, fixed_point, indices=indices)
    else
        error("invalid mehtod {String}: should be: \"tracking\" or \"tpsa\".")
    end
end

function _tracking_find_matrix66(accelerator::Accelerator, fixed_point::Union{Pos{T},Nothing}=nothing; indices::String="m66") where {T}
    
    if lowercase(indices) == "m66"
        return_tm_flag = false
        plus_idx = -1
    elseif lowercase(indices) == "open"
        return_tm_flag = true
        plus_idx = -1
    elseif lowercase(indices) == "closed"
        return_tm_flag = true
        plus_idx = 0
    else
        error("invalid indices {String}: should be: \"m66\", \"open\" or \"closed\".")
    end

    radsts::BoolState = accelerator.radiation_on
    if radsts == full
        accelerator.radiation_on = on
    end

    if fixed_point === nothing
        fixed_point, _ = find_orbit6(accelerator)
    end

    # gen and track 12 particles
    scaling = 1e-8*[1, 1, 1, 1, 0, 0] + 1e-6*[0, 0, 0, 0, 1, 1]
    particles_in = [Pos(fixed_point[:]) for _ in 1:12]
    particles_tracked = Vector{Pos{Float64}}[]
    sig = 1
    for i in 1:12
        idx = Int(ceil(i/2))
        particles_in[i][idx] += scaling[idx] * sig / 2.0
        out, _, _ = line_pass(accelerator, particles_in[i], "closed")
        push!(particles_tracked, out)
        sig *= -1
    end

    # construct the 66 matrices from tracked particles
    leng = length(particles_tracked[1])
    mat = zeros(Float64, (6,6,leng))
    for k in 1:leng
        for i in 1:2:12
            pp = particles_tracked[i][k][:]
            pm = particles_tracked[i+1][k][:]
            idx = Int(ceil(i/2))
            v = (pp .- pm) ./ scaling[idx]
            mat[:, idx, k] .= v
        end
    end
    m66 = mat[:, :, end]
    mat = [mat[:, :, j] for j in 1:leng+plus_idx]

    # returns
    if return_tm_flag
        return m66, mat
    end
    return m66

end

function _tpsa_find_matrix66(accelerator::Accelerator, fixed_point::Union{Pos{T}, Nothing}=nothing; indices::String="m66") where {T}
    
    if lowercase(indices) == "m66"
        return_tm_flag = false
        plus_idx = -1
    elseif lowercase(indices) == "open"
        return_tm_flag = true
        plus_idx = -1
    elseif lowercase(indices) == "closed"
        return_tm_flag = true
        plus_idx = 0
    else
        error("invalid indices {String}: should be: \"m66\", \"open\" or \"closed\".")
    end

    radsts::BoolState = accelerator.radiation_on
    if radsts == full
        accelerator.radiation_on = on
    end

    if fixed_point === nothing
        fixed_point, _ = find_orbit6(accelerator)
    end

    map::Pos{Tpsa{6, 1, Float64}} = Pos(fixed_point[:], tpsa=true) 
    
    if return_tm_flag
        tm::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}[]
    end
 
    for elem in accelerator.lattice
        element_pass(elem, map, accelerator)
        if return_tm_flag
            mat = _extract_m66(map)
            push!(tm, mat)
        end
    end

    m66 = _extract_m66(map)
    if plus_idx == 0
        push!(tm, m66)
    end

    accelerator.radiation_on = radsts

    if return_tm_flag
        return m66, tm
    else
        return m66
    end
end


function _extract_m66(p::Pos{Tpsa{6, 1, Float64}})
    return [
        p.rx[1] p.rx[2] p.rx[3] p.rx[4] p.rx[5] p.rx[6];
        p.px[1] p.px[2] p.px[3] p.px[4] p.px[5] p.px[6];
        p.ry[1] p.ry[2] p.ry[3] p.ry[4] p.ry[5] p.ry[6];
        p.py[1] p.py[2] p.py[3] p.py[4] p.py[5] p.py[6];
        p.de[1] p.de[2] p.de[3] p.de[4] p.de[5] p.de[6];
        p.dl[1] p.dl[2] p.dl[3] p.dl[4] p.dl[5] p.dl[6]
    ]    
end