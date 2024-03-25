using ..Track.AcceleratorModule: Accelerator
using ..Track.Elements: Element
using LinearAlgebra

function _flatten(arg::Vector{T}) where T
    flattened = []
    for e in arg
        if isa(e, Vector)
            append!(flattened, _flatten(e))
        else
            push!(flattened, e)
        end
    end
    type = typeof(flattened[1])
    return type.(flattened)
end

function _process_output(values, isflat)
    if isflat
        vals = _flatten(values)
    else
        vals = copy(values)
    end
    if length(vals) == 1
        vals = vals[1]
    end
    return vals
end

function _process_args_errors(acc_len::Int, indices::Union{Int, Vector{Int}, Vector{Vector{Int}}, Nothing}, values::Union{Vector{T}, T, Nothing}) where T<:Real
    isflat::Bool = false
    if indices === nothing
        indcs = Vector(1:acc_len)
    elseif isa(indices, Int)
        indcs = [Int[indices]]
    elseif length(indices) > 0 && isa(indices[1], Int)
        indcs = [Int[ind] for ind in indices]
        isflat = true
    else 
        indcs = copy(indices)
    end
    
    if values === nothing
        vals = [zeros(Float64, length(ind)) for ind in indcs]
    elseif isa(values, Vector{T})
        vals = [fill(values[i], length(indcs[i])) for i in eachindex(indcs)]
    elseif isa(values, T)
        vals = [fill(values, length(ind)) for ind in indcs]
    else
        vals = copy(values)
    end

    if length(indcs) != length(vals)
        error("length of values differs from length of indices.")
    end

    return indcs, vals, isflat
end

function get_error_misalignment_x(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int, Nothing}=nothing)
    idcs, vals, isflat = _process_args_errors(length(acc), indices, nothing)
    for (i,segs) in enumerate(idcs)
        @inbounds misx = -0.5*(acc.lattice[segs[1]].t_in[1] - acc.lattice[segs[end]].t_out[1])
        vals[i] .+= misx
    end
    return _process_output(vals, isflat)
end

function set_error_misalignment_x(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs, vals_) in zip(idcs, vals)
        @inbounds yaw = 0.5*(acc.lattice[segs[1]].t_in[1] + acc.lattice[segs[end]].t_out[1])
        for (ind, val) in zip(segs, vals_)
            @inbounds acc.lattice[ind].t_in[1] = yaw - val
            @inbounds acc.lattice[ind].t_out[1] = yaw + val
        end
    end
end

function add_error_misalignment_x(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs,vals_) in zip(idcs, vals)
        for (ind, val) in zip(segs, vals_)
            @inbounds acc.lattice[ind].t_in[1] -= val
            @inbounds acc.lattice[ind].t_out[1] += val
        end
    end
end

function get_error_misalignment_y(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int, Nothing}=nothing)
    idcs, vals, isflat = _process_args_errors(length(acc), indices, nothing)
    for (i,segs) in enumerate(idcs)
        @inbounds misy = -0.5*(acc.lattice[segs[1]].t_in[3] - acc.lattice[segs[end]].t_out[3])
        vals[i] .+= misy
    end
    return _process_output(vals, isflat)
end

function set_error_misalignment_y(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs, vals_) in zip(idcs, vals)
        @inbounds pitch = 0.5*(acc.lattice[segs[1]].t_in[3] + acc.lattice[segs[end]].t_out[3])
        for (ind, val) in zip(segs, vals_)
            @inbounds acc.lattice[ind].t_in[3] = pitch - val
            @inbounds acc.lattice[ind].t_out[3] = pitch + val
        end
    end
end

function add_error_misalignment_y(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs,vals_) in zip(idcs, vals)
        for (ind, val) in zip(segs, vals_)
            @inbounds acc.lattice[ind].t_in[3] -= val
            @inbounds acc.lattice[ind].t_out[3] += val
        end
    end
end

function get_error_rotation_roll(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int, Nothing}=nothing)
    idcs, vals, isflat = _process_args_errors(length(acc), indices, nothing)
    for (i,segs) in enumerate(idcs)
        @inbounds elem = acc.lattice[segs[1]]
        if elem.angle != 0.0 && elem.length != 0.0
            rho = elem.length / elem.angle
            angle = asin(elem.polynom_a[0] * rho)
        else
            angle = asin(elem.r_in[13])
        end    
        vals[i] .+= angle
    end
    return _process_output(vals, isflat)
end

function set_error_rotation_roll(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs, vals_) in zip(idcs, vals)
        
        C, S = cos(vals_[1]), sin(vals_[1])
        rot = diagm([C, C, C, C, 1.0, 1.0])
        rot[1, 3], rot[2, 4], rot[3, 1], rot[4, 2] = S, S, -S, -S

        for ind in segs
            @inbounds elem = acc.lattice[ind]
            if elem.angle != 0.0 && elem.length != 0.0
                rho = elem.length / elem.angle
                # sin(teta)/rho:
                elem.polynom_a[1] = S/rho
                # (cos(teta)-1)/rho:
                elem.polynom_b[1] = (C - 1.0)/rho
            else
                elem.r_in = Vector{Float64}(reshape(rot, (36)))
                elem.r_out = Vector{Float64}(reshape(transpose(rot), (36)))
            end
        end
    end
end

function add_error_rotation_roll(acc::Accelerator, indices::Union{Vector{Vector{Int}}, Vector{Int}, Int}, values::Union{Vector{T}, T}) where T<:Real
    idcs, vals, _ = _process_args_errors(length(acc), indices, values)
    for (segs, vals_) in zip(idcs, vals)
        
        C, S = cos(vals_[1]), sin(vals_[1])
        rot = diagm([C, C, C, C, 1.0, 1.0])
        rot[1, 3], rot[2, 4], rot[3, 1], rot[4, 2] = S, S, -S, -S

        for ind in segs
            @inbounds elem = acc.lattice[ind]
            if elem.angle != 0.0 && elem.length != 0.0
                rho = elem.length / elem.angle
                orig_S = elem.polynom_a[1] * rho
                # look at bndpolysymplectic4pass:
                orig_C = elem.polynom_b[1] * rho + 1.0
                # sin(teta)/rho:
                elem.polynom_a[1] = (orig_S*C - orig_C*S) / rho
                # (cos(teta)-1)/rho:
                elem.polynom_b[1] = (orig_C*C - orig_S*S - 1.0) / rho
            else
                elem.r_in = Vector(reshape(rot * reshape(copy(elem.r_in), (6,6)), (36)))
                elem.r_out = Vector(reshape(reshape(copy(elem.r_out), (6,6)) * transpose(rot), (36)))
            end
        end
    end
end