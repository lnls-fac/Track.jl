# Accelerator.jl

import Base: !=, setproperty!, setfield!, getproperty, show, getindex
using Printf
using ..Auxiliary
using ..Constants: electron_rest_energy_eV, light_speed
using ..Elements: Element

export Accelerator, find_spos, find_indices

mutable struct Accelerator
    energy::Float64
    cavity_state::BoolState
    radiation_state::BoolState
    vchamber_state::BoolState
    harmonic_number::Int
    lattice::Vector{Element}
    lattice_version::String
    length::Float64
    velocity::Float64
    beta_factor::Float64
    gamma_factor::Float64
    brho::Float64

    function Accelerator(energy::Float64)
        if energy < electron_rest_energy_eV
            @warn("The energy set ($energy) below the minimum (electron rest energy) $electron_rest_energy_eV [eV]")
            energy = electron_rest_energy_eV
        end
        gamma = energy/electron_rest_energy_eV
        beta = sqrt(1 - (1 / gamma^2))
        velocity = beta * light_speed
        brho = beta * energy / light_speed
        return new(energy, Auxiliary.off, Auxiliary.off, Auxiliary.off, 0, Element[], "", 0.0, velocity, beta, gamma, brho)
    end
end    

function isequal(acc1::Accelerator, acc2::Accelerator)
    if (acc1.energy != acc2.energy) return false end
    if (acc1.cavity_state != acc2.cavity_state) return false end
    if (acc1.radiation_state != acc2.radiation_state) return false end
    if (acc1.vchamber_state != acc2.vchamber_state) return false end
    if (acc1.harmonic_number != acc2.harmonic_number) return false end
    if (acc1.lattice_version != acc2.lattice_version) return false end
    #if (acc1.lattice != acc2.lattice) return false end revisar comparador de lattices
    return true
end

function Base.:(!=)(acc1::Accelerator, acc2::Accelerator)
    return !(isequal(acc1, acc2))
end

function update_cavity(accelerator::Accelerator)
    cavity_indices = find_cav_indices(accelerator)
    for index in cavity_indices
        cav = accelerator.lattice[index]
        if accelerator.cavity_state == on
            cav.pass_method = Auxiliary.pm_cavity_pass
        elseif accelerator.cavity_state == off && cav.length == 0.0
            cav.pass_method = Auxiliary.pm_identity_pass
        else
            cav.pass_method = Auxiliary.pm_drift_pass
        end
    end
end

function setproperty!(accelerator::Accelerator, symbol::Symbol, value)
    if symbol == :energy
        adjust_beam_parameters(accelerator, :energy, value)    
    
    elseif symbol == :cavity_state
        if isa(value, BoolState) || isa(value, Int) || isa(value, Bool)
            val = Int(value)
            if 0 <= val <= 1
                setfield!(accelerator, :cavity_state, BoolState(val))
                update_cavity(accelerator)
            else
                error("cavity_state should be 0(cavity off) or 1(cavity on)")
            end
        end
    
    elseif symbol == :radiation_state
        if isa(value, BoolState) || isa(value, Int) || isa(value, Bool)
            val = Int(value)
            if 0 <= val <= 2
                setfield!(accelerator, :radiation_state, BoolState(val))
            else
                error("radiation_state should be 0(radiation off), 1(radiation dumping) or 2(radiation full)")
            end
        end
    
    elseif symbol == :vchamber_state
        if isa(value, BoolState) || isa(value, Int) || isa(value, Bool)
            val = Int(value)
            if 0 <= val <= 2
                setfield!(accelerator, :vchamber_state, BoolState(val))
            else
                error("vchamber_state should be: (0, off, false) or (1, on, true)")
            end
        end
    elseif symbol == :harmonic_number
        setfield!(accelerator, :harmonic_number, value)

    elseif symbol == :lattice
        setfield!(accelerator, :lattice, value)
        len = find_spos(accelerator, indices="closed")[end] # closed lattice
        setfield!(accelerator, :length, len)

    elseif symbol == :lattice_version
        @warn("Changing the \"lattice_version\" manually is not recommended... But if its necessary, use: \"setfield!(Accelerator, :latice_version, value)\"")
    
    elseif symbol == :length 
        @warn("Cant manually change the \"length\". Consider changing the accelerator's lattice.")
    
    elseif symbol == :velocity 
        adjust_beam_parameters(accelerator, :velocity, value)

    elseif symbol == :beta_factor 
        adjust_beam_parameters(accelerator, :beta_factor, value)

    elseif symbol == :gamma_factor 
        adjust_beam_parameters(accelerator, :gamma_factor, value)

    elseif symbol == :brho 
        adjust_beam_parameters(accelerator, :brho, value)
    else
        throw(ArgumentError("Field $symbol is not a valid field for Accelerator"))
    end
end

function find_spos(accelerator::Accelerator; indices::T="open") where T<:Union{String, Vector{Int}}
    spos::Vector{Float64} = Float64[0.0]
    s_temp = 0.0
    for elem in accelerator.lattice
        s_temp += elem.length
        push!(spos, s_temp)
    end
    if isa(indices, String)
        if indices == "open"
            return spos[1:end-1]
        elseif indices == "closed"
            return spos
        else
            throw(ArgumentError("invalid indices: should be (String)\"closed\" or \"open\" or (Vector{Int})"))
        end
    elseif isa(indices, Vector{Int})
        if all([(1<=i<=accelerator.length) for i in indices])
            return spos[indices]
        else
            leng = length(accelerator.lattice)
            throw(ArgumentError("invalid index in argument \"indices\": should stay between 1 and $leng"))
        end
    else
        throw(ArgumentError("invalid indices: should be (String)\"closed\" or \"open\" or (Vector{Int})"))
    end
end

function find_indices(accelerator::Accelerator, property::String, value::Union{Real, String, Auxiliary.PassMethod})
    indices::Vector{Int} = Int[]
    for (idx, element) in enumerate(accelerator.lattice)
        if getfield(element, Symbol(property)) == value
            push!(indices, idx)
        end
    end
    return indices
end

function find_cav_indices(accelerator::Accelerator)
    idcs::Vector{Int} = Int[]
    for i in 1:1:length(accelerator.lattice) 
        if accelerator.lattice[i].frequency != 0.0
            push!(idcs, i)
        end
    end
    return idcs
end

function lattice_shift!(accelerator::Accelerator, index::Int)
    lattice::Vector{Element} = accelerator.lattice
    
    # Check if the index is within bounds
    if index < 1 || index > length(lattice)
        throw(ArgumentError("Index out of bounds"))
    end
    
    # Shift the lattice
    new_lattice = vcat(lattice[index:end], lattice[1:index-1])

    accelerator.lattice = new_lattice
end

function Base.show(io::IO, ::MIME"text/plain", accelerator::Accelerator)
    println(io, "------------------ Accelerator -----------------")
    @printf(io, "\tEnergy:            %1.3e  [GeV] \n", (accelerator.energy/1e9))
    println(io, "\tCavity State:      ", accelerator.cavity_state)
    println(io, "\tRadiation State:   ", accelerator.radiation_state)
    println(io, "\tVchamber State:    ", accelerator.vchamber_state)
    println(io, "\tHarmonic Number:   ", accelerator.harmonic_number)
    @printf(io, "\tLength:            %f  [m]\n", (accelerator.length))
    @printf(io, "\tVelocity:          %.8f  [m/s]\n", accelerator.velocity)
    @printf(io, "\tBeta Factor:       %.8f \n", accelerator.beta_factor)
    @printf(io, "\tGamma Factor:      %.8f \n", accelerator.gamma_factor)
    @printf(io, "\tBrho:              %.8f \n", accelerator.brho)
    if (accelerator.lattice_version != "")
        println(io, "\tLattice Version:            $(accelerator.lattice_version)")
    else
        println(io, "\tLattice Version:            (none)")
    end
    println(io, "------------------------------------------------")
end

function Base.show(io::IO, accelerator::Accelerator)
    println(io, "------------------ Accelerator -----------------")
    @printf(io, "\tEnergy:            %1.3e  [GeV] \n", (accelerator.energy/1e9))
    println(io, "\tCavity State:      ", accelerator.cavity_state)
    println(io, "\tRadiation State:   ", accelerator.radiation_state)
    println(io, "\tVchamber State:    ", accelerator.vchamber_state)
    println(io, "\tHarmonic Number:   ", accelerator.harmonic_number)
    @printf(io, "\tLength:            %f  [m]\n", (accelerator.length))
    @printf(io, "\tVelocity:          %.8f  [m/s]\n", accelerator.velocity)
    @printf(io, "\tBeta Factor:       %.8f \n", accelerator.beta_factor)
    @printf(io, "\tGamma Factor:      %.8f \n", accelerator.gamma_factor)
    @printf(io, "\tBrho:              %.8f \n", accelerator.brho)
    if (accelerator.lattice_version != "")
        println(io, "\tLattice Version:            $(accelerator.lattice_version)")
    else
        println(io, "\tLattice Version:            (none)")
    end
    println(io, "------------------------------------------------")
end

function adjust_beam_parameters(accelerator::Accelerator, property::Symbol, value::Float64)
    if property == :energy
        energy = value
        gamma_factor = energy / electron_rest_energy_eV
        beta_factor = sqrt(1.0 - (1.0 / (gamma_factor^2)))
        velocity = beta_factor * light_speed
        brho = gamma_factor * beta_factor * light_speed 
    elseif property == :gamma_factor
        gamma_factor = value
        energy = gamma_factor * electron_rest_energy_eV
        beta_factor = sqrt(1.0 - (1.0 / (gamma_factor^2)))
        velocity = beta_factor * light_speed
        brho = gamma_factor * beta_factor * light_speed 
    elseif property == :beta_factor
        beta_factor = value
        gamma_factor = 1.0 / sqrt(1.0 - (value^2))
        energy = gamma_factor * electron_rest_energy_eV
        velocity = beta_factor * light_speed
        brho = gamma_factor * beta_factor * light_speed 
    elseif property == :velocity
        velocity = value
        beta_factor = value / light_speed
        gamma_factor = 1.0 / sqrt(1.0 - (beta_factor^2))
        energy = gamma_factor * electron_rest_energy_eV
        brho = gamma_factor * beta_factor * light_speed 
    elseif property == :brho
        brho = value
        k = brho / (electron_rest_energy_eV / light_speed)
        gamma_factor = sqrt(1.0 + (k^2))
        energy = gamma_factor * electron_rest_energy_eV
        beta_factor = sqrt(1.0 - (1.0 / (gamma_factor^2)))
        velocity = beta_factor * light_speed
    else
        error("invalid property: $property")
    end
    setfield!(accelerator, :energy, energy)
    setfield!(accelerator, :gamma_factor, gamma_factor)
    setfield!(accelerator, :beta_factor, beta_factor)
    setfield!(accelerator, :velocity, velocity)
    setfield!(accelerator, :brho, brho)
end

function Base.getindex(accelerator::Accelerator, index::Int)
    return accelerator.lattice[index]
end

function Base.getindex(accelerator::Accelerator, ::Colon)
    return accelerator.lattice
end

function Base.length(accelerator::Accelerator)
    return length(accelerator.lattice)
end

function Base.size(accelerator::Accelerator)
    return length(accelerator.lattice)
end

Base.iterate(a::Accelerator, state=1) = state > length(a.lattice) ? nothing : (a.lattice[state], state + 1)