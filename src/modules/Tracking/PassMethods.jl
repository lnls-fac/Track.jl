
import Base: !=, getproperty, setfield!, setproperty!
using Printf
using ..AcceleratorModule: Accelerator
using ..Auxiliary: full, no_plane, off, on, pm_bnd_mpole_symplectic4_pass, pm_cavity_pass,
    pm_corrector_pass, pm_drift_pass, pm_identity_pass, pm_str_mpole_symplectic4_pass,
    st_success
using ..Elements: Element
using ..PosModule: Pos
using ..Constants: reduced_planck_constant, electron_charge, electron_mass, light_speed

function pow3(x::T) where T
    return (x*x)*x
end

function pow2(x::T) where T
    return (x*x)
end

const DRIFT1 ::Float64  = +0.6756035959798286638e00
const DRIFT2 ::Float64  = -0.1756035959798286639e00
const KICK1  ::Float64  = +0.1351207191959657328e01
const KICK2  ::Float64  = -0.1702414383919314656e01

# ATCOMPATIBLE
const TWOPI     ::Float64 = 6.28318530717959e0 # # not ATCOMPATIBLE :::: 2e0*pi
const CGAMMA    ::Float64 = 8.846056192e-05   # cgamma, [m]/[GeV^3] Ref[1] (4.1)
const M0C2      ::Float64 = 5.10999060e5        # Electron rest mass [eV]
const LAMBDABAR ::Float64 = 3.86159323e-13 # Compton wavelength/2pi [m]
const CER       ::Float64 = 2.81794092e-15       # Classical electron radius [m]
const CU        ::Float64 = 1.323094366892892e0     # 55/(24*sqrt(3)) factor
const CQEXT     ::Float64 = sqrt(CU * CER * reduced_planck_constant * electron_charge * light_speed) * electron_charge * electron_charge / pow3(electron_mass*light_speed*light_speed) #  for quant. diff. kick

function _drift(pos::Pos{T}, length::Float64) where T
    pnorm::T = 1 / (1 + pos.de)
    norml::T = length * pnorm
    pos.rx += norml * pos.px
    pos.ry += norml * pos.py
    pos.dl += 0.5 * norml * pnorm * (pos.px*pos.px + pos.py*pos.py)
end

function _calcpolykick(pos::Pos{T}, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64})  where T
    real_sum::T = 0.0
    imag_sum::T = 0.0
    n::Int = min(length(polynom_b), length(polynom_a))
    if n != 0
        real_sum = polynom_b[n]
        imag_sum = polynom_a[n]
        for i = n-1:-1:1
            real_sum_tmp = real_sum * pos.rx - imag_sum * pos.ry + polynom_b[i]
            imag_sum = imag_sum * pos.rx + real_sum * pos.ry + polynom_a[i]
            real_sum = real_sum_tmp
        end
    end
    return real_sum, imag_sum
end

function _b2_perp(bx::T, by::T, px::T, py::T, curv::S=1.0) where {T,S}
    curv2::S = curv*curv
    v_norm2_inv::T = curv2 + (px*px) + (py*py)
    b2p::T = (by*by) + (bx*bx)
    b2p *= curv2
    b2p += pow2((bx * py) - (by * px)) # with tpsa, i didnt checked the speed: pow(x, 2) or x*x
    b2p /= v_norm2_inv
    return b2p
end

function _strthinkick(pos::Pos{T}, length::Float64, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64}, rad_const::Float64=0.0, qexcit_const::Float64=0.0) where T

    real_sum::T = 0.0
    imag_sum::T = 0.0
    real_sum, imag_sum = _calcpolykick(pos, polynom_a, polynom_b)

    if rad_const != 0.0
        pnorm::T = 1.0 / (1.0 + pos.de)
        px::T = pos.px * pnorm
        py::T = pos.py * pnorm
        b2p::T = 0.0
        b2p = _b2_perp(imag_sum, real_sum, px, py)
        delta_factor::T  = pow2(1.0 + pos.de)
        dl_ds::T  = 1.0 + ((px*px + py*py) / 2.0)
        pos.de -= rad_const * delta_factor * b2p * dl_ds * length

        if qexcit_const != 0.0
            # quantum excitation kick
            d::T = delta_factor * qexcit_const * sqrt(sqrt(pow3(b2p)) * dl_ds)
            pos.de += d * randn()
        end

        pnorm = 1.0 + pos.de  # actually, this is the inverse of pnorm
        pos.px = px * pnorm
        pos.py = py * pnorm
    end

    pos.px -= length * real_sum
    pos.py += length * imag_sum
end

function _bndthinkick(pos::Pos{T}, length::Float64, polynom_a::Vector{Float64},
    polynom_b::Vector{Float64}, irho::Float64, rad_const::Float64=0.0, qexcit_const::Float64=0.0) where T
    
    real_sum::T = 0.0
    imag_sum::T = 0.0
    real_sum, imag_sum = _calcpolykick(pos, polynom_a, polynom_b)
    de::T = pos.de

    if rad_const != 0.0
        pnorm::T = 1.0 / (1.0 + pos.de)
        px::T = pos.px * pnorm
        py::T = pos.py * pnorm
        curv::T = 1.0 + (irho * pos.rx)
        b2p::T = 0.0
        b2p = _b2_perp(imag_sum, real_sum+irho, px, py, curv)
        delta_factor::T = pow2(1.0 + pos.de)
        dl_ds::T = curv + ((px*px + py*py) / 2.0)
        pos.de -= rad_const * delta_factor * b2p * dl_ds * length

        if qexcit_const != 0.0
            # quantum excitation kick
            d::T = delta_factor * qexcit_const * sqrt(sqrt(pow3(b2p)) * dl_ds)
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

function _edge_fringe(pos::Pos{T}, inv_rho::Float64, edge_angle::Float64,
    fint::Float64, gap::Float64) where T
    de::T = pos.de
    rx::T = pos.rx
    ry::T = pos.ry

    fx::T = inv_rho * tan(edge_angle) / (1.0 + de)

    psi_bar::T = edge_angle - inv_rho * gap * fint * (1.0 + pow2(sin(edge_angle))) / cos(edge_angle) / (1.0 + de)
    
    fy::T = inv_rho * tan(psi_bar) / (1.0 + de)
    
    pos.px += rx * fx
    pos.py -= ry * fy
end

function pm_identity_pass!(pos::Pos{T}, element::Element) where T
    return st_success
end

function pm_drift_pass!(pos::Pos{T}, element::Element) where T
    _drift(pos, element.length)
    return st_success
end

function pm_str_mpole_symplectic4_pass!(pos::Pos{T}, elem::Element, accelerator::Accelerator) where T

    steps::Int = elem.nr_steps

    sl::Float64 = elem.length / Float64(steps)
    l1::Float64 = sl * DRIFT1
    l2::Float64 = sl * DRIFT2
    k1::Float64 = sl * KICK1
    k2::Float64 = sl * KICK2
    polynom_b::Vector{Float64} = elem.polynom_b
    polynom_a::Vector{Float64} = elem.polynom_a
    rad_const::Float64 = 0.0
    qexcit_const::Float64 = 0.0

    if accelerator.radiation_state == on
        rad_const = CGAMMA * pow3(accelerator.energy/1.0e9) / TWOPI
    end

    if accelerator.radiation_state == full
        qexcit_const = CQEXT * pow2(accelerator.energy) * sqrt(accelerator.energy * sl)
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

function pm_bnd_mpole_symplectic4_pass!(pos::Pos{T}, elem::Element, accelerator::Accelerator) where T

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
        rad_const = CGAMMA * pow3(accelerator.energy / 1e9) / TWOPI
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

    for i in 1:steps
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

function pm_corrector_pass!(pos::Pos{T}, elem::Element) where T

    xkick::Float64 = elem.hkick
    ykick::Float64 = elem.vkick

    if elem.length == 0.0
        pos.px += hkick
        pos.py += vkick
    else
        px::T = pos.px
        py::T = pos.py
        de::T = pos.de
        pnorm::T = 1.0 / (1.0 + de)
        norml::T = elem.length * pnorm
        pos.dl += norml * pnorm * 0.5 * (xkick * xkick/3.0 + ykick * ykick/3.0 + px*px + py*py + px * xkick + py * ykick)
        pos.rx += norml * (px + 0.5 * xkick)
        pos.px += xkick
        pos.ry += norml * (py + 0.5 * ykick)
        pos.py += ykick
    end

    return st_success
end

function pm_cavity_pass!(pos::Pos{T}, elem::Element, accelerator::Accelerator, turn_number::Int) where T
    if accelerator.cavity_state == off
        return pm_drift_pass!(pos, elem)
    end

    nv::Float64 = elem.voltage / accelerator.energy
    philag::Float64 = elem.phase_lag
    frf::Float64 = elem.frequency
    harmonic_number::Int = accelerator.harmonic_number

    # not ATCOMPATIBLE
    #velocity::Float64 = accelerator.velocity / 1e8 # numerical problem
    # ATCOMPATIBLE
    velocity::Float64 = light_speed #/ 1e8

    # L0::Float64 = accelerator.length
    # factor::Float64 = (velocity*harmonic_number/frf*1e8 - L0) / velocity / 1e8

    if elem.length == 0.0
        #pos.de += -nv * sin((TWOPI * frf * ((pos.dl/velocity/1e8) - (factor*turn_number))) - philag)
        fd::Float64 = TWOPI * frf * pos.dl / velocity  - philag
        # @printf("fd = %+.100e\n", fd)
        # @printf("sin(fd) = %+.100e\n", sin(fd))
        # @printf("-nv * sin(fd) = %+.100e\n", -nv*sin(fd))
        # @printf("antes  -- de = %+.100e\n", pos.de)
        pos.de += -nv * mysin(TWOPI * frf * pos.dl / velocity  - philag)
        # @printf("depois -- de = %+.100e\n", pos.de)
    else
        px::T = pos.px
        py::T = pos.py

        # Drift half length
        pnorm::T = 1.0 / (1.0 + pos.de)
        norml::T = (0.5 * elem.length) * pnorm
        pos.rx += norml * px
        pos.ry += norml * py
        pos.dl += 0.5 * norml * pnorm * (px*px + py*py)

        # Longitudinal momentum kick
        # pos.de += -nv * sin((TWOPI * frf * ((pos.dl/velocity/1e8) - (factor*turn_number))) - philag)
        pos.de += -nv * mysin(TWOPI * frf * pos.dl / velocity  - philag)

        # Drift half length
        pnorm = 1.0 / (1.0 + pos.de)
        norml = (0.5 * elem.length) * pnorm
        pos.rx += norml * px
        pos.ry += norml * py
        pos.dl += 0.5 * norml * pnorm * (px*px + py*py)
    end

    return st_success
end

function mysin(x::T) where T
    if isa(x, Number)
        return sinpi(x/pi)
    else
        return sin(x)
    end
end