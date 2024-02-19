"""Elements module."""
module Elements

include("Elements.jl")

using ..Auxiliary

export corrector, marker, hcorrector, 
        vcorrector, drift, rbend, 
        quadrupole, sextupole, rfcavity
        # matrix, kickmap

# Define functions for Element
function marker(fam_name::String)
    element = Element(fam_name)
    element.pass_method = pm_identity_pass
    return element
end

function corrector(fam_name::String, hkick::Float64, vkick::Float64)
    element = Element(fam_name)
    element.pass_method = pm_corrector_pass
    element.hkick = hkick
    element.vkick = vkick
    return element
end

function hcorrector(fam_name::String, hkick::Float64)
    element = Element(fam_name)
    element.pass_method = pm_corrector_pass
    element.hkick = hkick
    element.vkick = 0.0
    return element
end

function vcorrector(fam_name::String, vkick::Float64)
    element = Element(fam_name)
    element.pass_method = pm_corrector_pass
    element.hkick = 0.0
    element.vkick = vkick
    return element
end

function drift(fam_name::String, length::Float64)
    element = Element(fam_name)
    element.pass_method = pm_drift_pass
    element.length = length
    return element
end

# # Revisar implementacao
# function matrix(fam_name::String, length::Float64)
#     element = Element(fam_name)
#     element.pass_method = pm_matrix_pass
#     element.length = length
#     return element
# end

function rbend(fam_name::String, length::Float64, angle::Float64, angle_in::Float64, angle_out::Float64,
                gap::Float64, fint_in::Float64, fint_out::Float64, polynom_a::Vector{Float64},
                polynom_b::Vector{Float64}, K::Float64=-999.0, S::Float64=-999.0, nr_steps::Int=20)
    element = Element(fam_name)
    element.pass_method = pm_bnd_mpole_symplectic4_pass
    element.length = length
    element.angle = angle
    element.angle_in = angle_in
    element.angle_out = angle_out
    element.gap = gap
    element.fint_in = fint_in
    element.fint_out = fint_out
    element.polynom_a = polynom_a
    element.polynom_b = polynom_b
    if (K != -999.0) 
        element.polynom_b[2] = K
    end
    if (S != -999.0) 
        element.polynom_b[3] = S
    end
    element.nr_steps = nr_steps
    return element
end

function quadrupole(fam_name::String, length::Float64, K::Float64; nr_steps::Int=10)
    element = Element(fam_name)
    element.pass_method = pm_str_mpole_symplectic4_pass
    element.polynom_a = [0.0, 0.0, 0.0] # copy(default_polynom)
    element.polynom_b = [0.0, 0.0, 0.0] # copy(default_polynom) 
    element.polynom_b[2] = K
    element.nr_steps = nr_steps
    element.length = length
    return element
end

function sextupole(fam_name::String, length::Float64, S::Float64; nr_steps::Int=5)
    element = Element(fam_name)
    element.pass_method = pm_str_mpole_symplectic4_pass
    element.polynom_a = [0.0, 0.0, 0.0] # copy(default_polynom)
    element.polynom_b = [0.0, 0.0, 0.0] # copy(default_polynom) 
    element.polynom_b[3] = S
    element.nr_steps = nr_steps
    element.length = length
    return element
end

function rfcavity(fam_name::String, frequency::Float64, voltage::Float64, phase_lag::Float64, length::Float64=0.0)  
    element = Element(fam_name)
    element.pass_method = pm_cavity_pass
    element.frequency = frequency
    element.voltage = voltage
    element.phase_lag = phase_lag
    element.length = length
    return element
end

# # Revisar implementacao
# function kickmap(fam_name::String, kicktable_idx::Int, nr_steps::Int, rescale_kicks::Float64)
#     element = Element(fam_name)
#     element.pass_method = pm_kickmap_pass
#     element.nr_steps = nr_steps
#     element.kicktable_idx = kicktable_idx
#     element.rescale_kicks = rescale_kicks
#     return element
# end

function flatten_to_element_vector(arg::Any)
    if isa(arg, Vector{Element})
        return arg
    elseif isa(arg, Element)
        return [arg]
    elseif isa(arg, Vector)
        # Recursively flatten each element in the vector
        flattened = []
        for elem in arg
            append!(flattened, flatten_to_element_vector(elem))
        end
        return flattened
    else
        throw(ArgumentError("Unsupported argument type: $(typeof(arg))"))
    end
end

function lattice_flatten!(arg::Any)
    latt = Vector{Element}(flatten_to_element_vector(arg))
    return latt
end

end # module Elements