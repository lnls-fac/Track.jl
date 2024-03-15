

export element_pass, line_pass, ring_pass

using ..AcceleratorModule: Accelerator
using ..Auxiliary: PassMethod, Plane, Status, no_plane, on, plane_x, plane_xy, plane_y,
    pm_bnd_mpole_symplectic4_pass, pm_cavity_pass, pm_corrector_pass, pm_drift_pass,
    pm_identity_pass, pm_str_mpole_symplectic4_pass, st_particle_lost, st_success,
    vchamber_ellipse, vchamber_rectangle, vchamber_rhombus, BoolState
using ..Elements: Element
using ..PosModule: Pos
using ..TPSA: Tpsa

function element_pass(
    element::Element,            # element through which to track particle
    particle::Pos{T},      # initial electron coordinates
    accelerator::Accelerator;    # accelerator parameters
    turn_number::Int = 0         # optional turn number parameter
    ) where T
    status::Status = st_success

    pass_method::PassMethod = element.pass_method

    if pass_method == pm_identity_pass
        status = pm_identity_pass!(particle, element)

    elseif pass_method == pm_drift_pass
        status = pm_drift_pass!(particle, element)

    elseif pass_method == pm_str_mpole_symplectic4_pass
        status = pm_str_mpole_symplectic4_pass!(particle, element, accelerator)

    elseif pass_method == pm_bnd_mpole_symplectic4_pass
        status = pm_bnd_mpole_symplectic4_pass!(particle, element, accelerator)

    elseif pass_method == pm_corrector_pass
        status = pm_corrector_pass!(particle, element)

    elseif pass_method == pm_cavity_pass
        status = pm_cavity_pass!(particle, element, accelerator, turn_number)

    else
        return st_passmethod_not_defined
    end

    return status
end

function line_pass(
    accelerator::Accelerator,
    particle::Pos{T},
    indices::Vector{Int};
    element_offset::Int = 1,
    turn_number::Int = 0
    ) where T
    leng::Float64 = length(accelerator.lattice)
    if any([!(1<=i<=leng+1) for i in indices])
        error("invalid indices: outside of lattice bounds. The valid indices should stay between 1 and $(leng+1)")
    end
    if element_offset > leng || element_offset < 1
        error("invalid indices: outside of lattice bounds. The valid indices should stay between 1 and $leng")
    end
    tpsa=false
    if isa(particle.rx, Tpsa)
        tpsa = true
    end

    status::Status = st_success
    lostplane::Plane = no_plane
    lostelement::Int64 = 0
    tracked_pos::Vector{Pos{T}} = Pos{T}[]

    line::Vector{Element} = accelerator.lattice
    nr_elements::Int = length(line)

    # Create vector of booleans to determine when to store position
    indcs::Vector{Bool} = falses(nr_elements + 1)
    indcs[indices] .= true

    for i in 1:nr_elements
        # Read-only access to element object parameters
        element::Element = line[element_offset]

        # Stores trajectory at entrance of each element
        if indcs[i]
            push!(tracked_pos, copy(particle))
        end

        status = element_pass(element, particle, accelerator, turn_number=turn_number)

        status, lostplane = _check_if_lost!(element, particle.rx, particle.ry, accelerator.vchamber_state)

        if status != st_success
            nan_particle::Pos{T} = Pos(NaN,tpsa=tpsa)
            lostelement = i
            # Fill the rest of vector with NaNs
            for j in i+1:length(indcs)
                if indcs[j]
                    push!(tracked_pos, copy(nan_particle))
                end
            end
            particle.rx = nan_particle.rx 
            particle.ry = nan_particle.ry
            particle.px = nan_particle.px
            particle.py = nan_particle.py
            particle.de = nan_particle.de
            particle.dl = nan_particle.dl
            return tracked_pos, status, lostplane, lostelement
        end

        # Moves to the next element index
        element_offset = mod1(element_offset + 1, nr_elements)
    end

    # Stores final particle position at the end of the line
    if indcs[nr_elements+1]
        push!(tracked_pos, copy(particle))
    end
    particle = copy(tracked_pos[end])
    return tracked_pos, status, lostplane, lostelement
end

function line_pass(
    accelerator::Accelerator,
    particle::Pos{T},
    indices::String = "closed";
    element_offset::Int = 1,
    turn_number::Int = 0
    ) where T
    leng::Int = length(accelerator.lattice)
    idcs::Vector{Int} = Int[]
    if indices == "closed"
        idcs = Int[i for i in 1:1:leng+1]
    elseif indices == "open"
        idcs = Int[i for i in 1:1:leng]
    elseif indices == "end"
        idcs = Int[leng+1]
    else
        error("invalid indices: should be a String:(closed, open, end) or a Vector{Int})")
    end
    return line_pass(accelerator, particle, idcs, element_offset=element_offset, turn_number=turn_number)
end

function ring_pass(accelerator::Accelerator, 
    particle::Pos{T}, 
    nr_turns::Int = 1; 
    element_offset::Int = 1,
    turn_by_turn::Bool = false
    ) where T
    if nr_turns<1 
        error("invalid nr_turns: should be >= 1")
    end
    if turn_by_turn
        v = Pos{T}[copy(particle)]
    end
    lostplane = no_plane
    st = st_success
    lostturn = -1
    lostelement = 0
    for turn in 1:nr_turns
        _, st, lostplane, lostelement = line_pass(accelerator, particle, "end", element_offset=element_offset, turn_number=turn-1)
        if st == st_success
            if turn_by_turn
                push!(v, copy(particle))
            end
        else
            lostturn = turn
            if turn_by_turn
                append!(v, [copy(particle) for i in 1:1:(nr_turns-turn+1)])
            else
                v = particle
            end
            return v, st, lostplane, lostturn, lostelement
        end
    end
    if !turn_by_turn
        v = particle
    end
    return v, st, lostplane, lostturn, lostelement
end

function _check_if_lost!(element::Element, x::T, y::T, vchamber_state::BoolState) where T
    status::Status = st_success
    lost_plane::Plane = no_plane
    if !isfinite(x)
        lost_plane = plane_x
        status = st_particle_lost
    end

    if !isfinite(y)
        if status != st_particle_lost
            lost_plane = plane_y
            status = st_particle_lost
        else
            lost_plane = plane_xy
        end
    end

    if (status != st_particle_lost) && (vchamber_state == on)
        status, lost_plane = _aux_check_if_lost!(element, x, y)
    end
    return status, lost_plane
end

function _aux_check_if_lost!(element::Element, x::T, y::T) where T
    status::Status = st_success
    lost_plane::Plane = no_plane
    if element.vchamber == vchamber_rectangle
        if x <= element.hmin || x >= element.hmax
            lost_plane = plane_x
            status = st_particle_lost
        end
        if y <= element.vmin || y >= element.vmax
            if status != st_particle_lost
                lost_plane = plane_y
                status = st_particle_lost
            else
                lost_plane = plane_xy
            end
        end
    else
        lx::Float64 = (element.hmax - element.hmin) / 2
        ly::Float64 = (element.vmax - element.vmin) / 2
        xc::Float64 = (element.hmax + element.hmin) / 2
        yc::Float64 = (element.vmax + element.vmin) / 2
        if lx == 0.0 || ly == 0.0
            println("largura = 0!!!")
        end
        xn::Float64 = abs((x - xc) / lx)
        yn::Float64 = abs((y - yc) / ly)
        amplitude::Float64 = 0.0

        if element.vchamber == vchamber_rhombus
            amplitude = xn + yn
        elseif element.vchamber == vchamber_ellipse
            amplitude = xn^2 + yn^2
        else
            amplitude = xn^Int(element.vchamber) + yn^Int(element.vchamber)
        end

        if amplitude > 1
            status = st_particle_lost
            lost_plane = plane_xy
        else
            status = st_success
            lost_plane = no_plane
        end
    end
    return status, lost_plane
end