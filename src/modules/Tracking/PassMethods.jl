
import Base: !=, getproperty, setfield!, setproperty!

using ..AcceleratorModule: Accelerator
using ..Auxiliary: full, no_plane, off, on, pm_bnd_mpole_symplectic4_pass, pm_cavity_pass,
    pm_corrector_pass, pm_drift_pass, pm_identity_pass, pm_str_mpole_symplectic4_pass,
    st_success
using ..Elements: Element
using ..PosModule: Pos
using ..Constants: reduced_planck_constant, electron_charge, electron_mass, light_speed

const DRIFT1::Float64  =  0.6756035959798286638e00
const DRIFT2::Float64  = -0.1756035959798286639e00
const KICK1::Float64   =  0.1351207191959657328e01
const KICK2::Float64   = -0.1702414383919314656e01

const TWOPI::Float64 = 2*pi               # 2*pi
const CGAMMA::Float64 = 8.846056192e-05   # cgamma, [m]/[GeV^3] Ref[1] (4.1)
const M0C2::Float64 = 5.10999060e5        # Electron rest mass [eV]
const LAMBDABAR::Float64 = 3.86159323e-13 # Compton wavelength/2pi [m]
const CER::Float64 = 2.81794092e-15       # Classical electron radius [m]
const CU::Float64 = 1.323094366892892     # 55/(24*sqrt(3)) factor
const CQEXT::Float64 = sqrt(CU * CER * reduced_planck_constant * electron_charge * light_speed) * electron_charge * electron_charge / ((electron_mass*light_speed*light_speed)^3) #  for quant. diff. kick

function _drift(pos::Pos{Float64}, length::Float64) 
    pnorm::Float64 = 1 / (1 + pos.de)
    norml::Float64 = length * pnorm
    pos.rx += norml * pos.px
    pos.ry += norml * pos.py
    pos.dl += 0.5 * norml * pnorm * (pos.px*pos.px + pos.py*pos.py)
end

function _calcpolykick(pos::Pos{Float64}, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64}) 
    n::Int = min(length(polynom_b), length(polynom_a))
    if n != 0
        real_sum = polynom_b[n]
        imag_sum = polynom_a[n]
        for i = n-1:-1:1
            real_sum_tmp = real_sum * pos.rx - imag_sum * pos.ry + polynom_b[i]
            imag_sum = imag_sum * pos.rx + real_sum * pos.ry + polynom_a[i]
            real_sum = real_sum_tmp
        end
        return real_sum, imag_sum
    end
    return 0.0, 0.0
end

function _b2_perp(bx::Float64, by::Float64, px::Float64, py::Float64, curv::Float64=1.0) 
    curv2::Float64 = curv^2
    v_norm2_inv::Float64 = curv2 + px^2 + py^2
    b2p::Float64 = by^2 + bx^2
    b2p *= curv2
    b2p += (bx * py - by * px)^2
    b2p /= v_norm2_inv
    return b2p
end

function _strthinkick(pos::Pos{Float64}, length::Float64, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64}, rad_const::Float64=0.0, qexcit_const::Float64=0.0) 

    real_sum::Float64 = 0.0
    imag_sum::Float64 = 0.0
    real_sum, imag_sum = _calcpolykick(pos, polynom_a, polynom_b)

    if rad_const != 0.0
        pnorm::Float64 = 1 / (1 + pos.de)
        px::Float64 = pos.px * pnorm
        py::Float64 = pos.py * pnorm
        b2p::Float64 = 0.0
        b2p = _b2_perp(imag_sum, real_sum, px, py)
        delta_factor::Float64  = (1 + pos.de)^2
        dl_ds::Float64  = 1.0 + ((px*px + py*py) / 2)
        pos.de -= rad_const * delta_factor * b2p * dl_ds * length

        if qexcit_const != 0.0
            # quantum excitation kick
            d::Float64 = delta_factor * qexcit_const * sqrt(b2p^1.5 * dl_ds)
            pos.de += d * randn()
        end

        pnorm = 1.0 + pos.de  # actually, this is the inverse of pnorm
        pos.px = px * pnorm
        pos.py = py * pnorm
    end

    pos.px -= length * real_sum
    pos.py += length * imag_sum
end

function _bndthinkick(pos::Pos{Float64}, length::Float64, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64}, irho::Float64, rad_const::Float64=0.0, qexcit_const::Float64=0.0) 
    
    real_sum::Float64 = 0.0
    imag_sum::Float64 = 0.0
    real_sum, imag_sum = _calcpolykick(pos, polynom_a, polynom_b)
    de::Float64 = pos.de

    if rad_const != 0.0
        pnorm::Float64 = 1 / (1 + pos.de)
        px::Float64 = pos.px * pnorm
        py::Float64 = pos.py * pnorm
        curv::Float64 = 1.0 + (irho * pos.rx)
        b2p::Float64 = 0.0
        b2p = _b2_perp(imag_sum, real_sum+irho, px, py, curv)
        delta_factor::Float64 = (1 + pos.de)^2
        dl_ds::Float64 = curv + ((px*px + py*py) / 2)
        pos.de -= rad_const * delta_factor * b2p * dl_ds * length

        if qexcit_const != 0.0
            # quantum excitation kick
            d::Float64 = delta_factor * qexcit_const * sqrt(b2p^1.5 * dl_ds)
            pos.de += d * randn()
        end

        pnorm = 1.0 + pos.de  # actually this is the inverse of pnorm
        pos.px = px * pnorm
        pos.py = py * pnorm
    end
    
    pos.px -= length * (real_sum - (de - pos.rx * irho) * irho)
    pos.py += length * imag_sum
    pos.dl += length * irho * pos.rx
end

function _edge_fringe(pos::Pos{Float64}, inv_rho::Float64, edge_angle::Float64,
    fint::Float64, gap::Float64) 
    de::Float64 = pos.de
    rx::Float64 = pos.rx
    ry::Float64 = pos.ry

    fx::Float64 = inv_rho * tan(edge_angle) / (1.0 + de)

    psi_bar::Float64 = edge_angle - inv_rho * gap * fint * (1 + sin(edge_angle)^2) / cos(edge_angle) / (1.0 + de)
    
    fy::Float64 = inv_rho * tan(psi_bar) / (1.0 + de)
    
    pos.px += rx * fx
    pos.py -= ry * fy
end

function pm_identity_pass!(pos::Pos{Float64}, element::Element) 
    return st_success
end

function pm_drift_pass!(pos::Pos{Float64}, element::Element) 
    _drift(pos, element.length)
    return st_success
end

function pm_str_mpole_symplectic4_pass!(pos::Pos{Float64}, elem::Element, accelerator::Accelerator) 

    steps::Int = elem.nr_steps

    sl::Float64 = elem.length / steps
    l1::Float64 = sl * DRIFT1
    l2::Float64 = sl * DRIFT2
    k1::Float64 = sl * KICK1
    k2::Float64 = sl * KICK2
    polynom_b::Vector{Float64} = elem.polynom_b
    polynom_a::Vector{Float64} = elem.polynom_a
    rad_const::Float64 = 0.0
    qexcit_const::Float64 = 0.0

    if accelerator.radiation_state == on
        rad_const = CGAMMA * (accelerator.energy/1e9)^3 / TWOPI
    end

    if accelerator.radiation_state == full
        qexcit_const = CQEXT * accelerator.energy^2 * sqrt(accelerator.energy * sl)
    end

    for i in 1:steps
        _drift(pos, l1)
        _strthinkick(pos, k1, polynom_a, polynom_b, rad_const, 0.0)
        _drift(pos, l2)
        _strthinkick(pos, k2, polynom_a, polynom_b, rad_const, qexcit_const)
        _drift(pos, l2)
        _strthinkick(pos, k1, polynom_a, polynom_b, rad_const, 0.0)
        _drift(pos, l1)
    end

    return st_success
end

function pm_bnd_mpole_symplectic4_pass!(pos::Pos{Float64}, elem::Element, accelerator::Accelerator) 

    steps::Int = elem.nr_steps

    sl   ::Float64 = elem.length / Float64(steps)
    l1   ::Float64 = sl * DRIFT1
    l2   ::Float64 = sl * DRIFT2
    k1   ::Float64 = sl * KICK1
    k2   ::Float64 = sl * KICK2
    irho ::Float64 = elem.angle / elem.length

    polynom_b::Vector{Float64} = elem.polynom_b
    polynom_a::Vector{Float64} = elem.polynom_a
    
    rad_const::Float64 = 0.0
    qexcit_const::Float64 = 0.0

    if accelerator.radiation_state == on
        rad_const = CGAMMA * (accelerator.energy / 1e9)^3 / (TWOPI)
    end

    if accelerator.radiation_state == full
        qexcit_const = CQEXT * accelerator.energy^2 * sqrt(accelerator.energy * sl)
    end

    ang_in   ::Float64 = elem.angle_in
    ang_out  ::Float64 = elem.angle_out
    fint_in  ::Float64 = elem.fint_in
    fint_out ::Float64 = elem.fint_out
    gap      ::Float64 = elem.gap

    _edge_fringe(pos, irho, ang_in, fint_in, gap)

    for i in 1:1:Int(steps)
        _drift(pos, l1)
        _bndthinkick(pos, k1, polynom_a, polynom_b, irho, rad_const, 0.0)
        _drift(pos, l2)
        _bndthinkick(pos, k2, polynom_a, polynom_b, irho, rad_const, qexcit_const)
        _drift(pos, l2)
        _bndthinkick(pos, k1, polynom_a, polynom_b, irho, rad_const, 0.0)
        _drift(pos, l1)
    end

    _edge_fringe(pos, irho, ang_out, fint_out, gap)

    return st_success
end

function pm_corrector_pass!(pos::Pos{Float64}, elem::Element) 

    xkick::Float64 = elem.hkick
    ykick::Float64 = elem.vkick

    if elem.length == 0
        pos.px += hkick
        pos.py += vkick
    else
        px::Float64 = pos.px
        py::Float64 = pos.py
        de = pos.de
        pnorm::Float64 = 1 / (1 + de)
        norml::Float64 = elem.length * pnorm
        pos.dl += norml * pnorm * 0.5 * (xkick * xkick/3.0 + ykick * ykick/3.0 + px*px + py*py + px * xkick + py * ykick)
        pos.rx += norml * (px + 0.5 * xkick)
        pos.px += xkick
        pos.ry += norml * (py + 0.5 * ykick)
        pos.py += ykick
    end

    return st_success
end

function pm_cavity_pass!(pos::Pos{Float64}, elem::Element, accelerator::Accelerator, turn_number::Int) 
    if accelerator.cavity_state == off
        return pm_drift_pass!(pos, elem)
    end

    nv::Float64 = elem.voltage / accelerator.energy
    philag::Float64 = elem.phase_lag
    frf::Float64 = elem.frequency
    harmonic_number::Int = accelerator.harmonic_number
    velocity::Float64 = accelerator.velocity
    # velocity = light_speed
    L0::Float64 = accelerator.length
    T0::Float64 = L0 / velocity
    #println(stdout,"\nfactor = ",harmonic_number/frf - T0)

    if elem.length == 0
        pos.de += -nv * sin((TWOPI * frf * ((pos.dl/velocity) - ((harmonic_number/frf - T0)*turn_number))) - philag)
        #pos.de += -nv * sin(TWOPI * frf * dl / velocity - philag)
    else
        px::Float64 = pos.px
        py::Float64 = pos.py

        # Drift half length
        pnorm::Float64 = 1 / (1 + pos.de)
        norml::Float64 = (0.5 * elem.length) * pnorm
        pos.rx += norml * px
        pos.ry += norml * py
        pos.dl += 0.5 * norml * pnorm * (px*px + py*py)

        # Longitudinal momentum kick
        pos.de += -nv * sin((TWOPI * frf * ((pos.dl/velocity) - ((harmonic_number/frf - T0)*turn_number))) - philag)
        # pos.de += -nv * sin(TWOPI * frf * dl / velocity - philag)

        # Drift half length
        pnorm = 1.0 / (1.0 + pos.de)
        norml = (0.5 * elem.length) * pnorm
        pos.rx += norml * px
        pos.ry += norml * py
        pos.dl += 0.5 * norml * pnorm * (px*px + py*py)
    end

    return st_success
end