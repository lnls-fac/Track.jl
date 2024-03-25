using ..Track.AcceleratorModule: Accelerator
using ..Track.Elements: Element

function get_error_misalignment_x(acc::Accelerator, indices::Union{Vector{Int}, Int, Nothing}=nothing)
    L = length(acc)
    if indices === nothing
        idx = Vector(1:L)
    elseif any([i<1 || i>L for i in indices])
        error("invalid indices")
    else 
        idx = copy(indices)
    end
    vals = Float64[]
    for i in idx
        @inbounds misx = -0.5*(acc.lattice[i].t_in[1] - acc.lattice[i].t_out[1])
        push!(vals, misx)
    end
    return vals
end


function add_error_misalignment_x!(acc::Accelerator, indices::Vector{Int}, values::Union{Vector{T}, T}) where T<:Real
    L = length(acc)
    if any([i<1 || i>L for i in indices])
        error("invalid indices")
    end
    if isa(values, Float64)
        vals = Float64[values for _ in 1:L]
    elseif length(values) == length(indices)
        vals = copy(values)
    else
        error("invalid values")
    end
    for i in eachindex(indices)
        @inbounds acc.lattice[indices[i]].t_in[1] -= vals[i]
        @inbounds acc.lattice[indices[i]].t_out[1] += vals[i]
    end
end