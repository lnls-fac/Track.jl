"""Track module."""
module Track

include("modules/Constants/constantsModule.jl")
include("modules/Auxiliary/auxiliaryModule.jl")
include("modules/Pos/posModule.jl")
include("modules/Elements/elementsModule.jl")
include("modules/Accelerator/acceleratorModule.jl")
include("modules/Tracking/trackingModule.jl")
include("modules/Orbit/orbitModule.jl")

using .PosModule: Pos
using .Elements: Element
using .AcceleratorModule: Accelerator, find_spos, find_indices
using .Tracking: element_pass, line_pass, ring_pass
using .Orbit: find_orbit4, find_orbit6

export Constants, 
        Auxiliary, 
        Pos, 
        find_spos, find_indices, 
        element_pass, line_pass, ring_pass,
        find_orbit4, find_orbit6

using PrecompileTools

@setup_workload begin
    v::Vector{Float64} = rand(Float64, 6) * 1e-6
    poly::Vector{Float64} = zeros(Float64, 3)
    @compile_workload begin
        p = Pos(v)
        marker = Elements.marker("marker")
        dip = Elements.rbend("dip", 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, copy(poly), copy(poly), -999.0, -999.0, 1)
        quad = Elements.quadrupole("quad", 0.02, 0.0, nr_steps=1)
        sext = Elements.sextupole("sext", 0.02, 0.0, nr_steps=1)
        cav = Elements.rfcavity("cav", 0.0, 1.0, 0.0, 0.0)
        drift = Elements.drift("drift", 0.0)
        lattice = Element[drift, marker, dip, quad, sext, cav]
        m = Accelerator(3e9)
        m.lattice = lattice
        m.radiation_state = Auxiliary.full
        m.cavity_state = Auxiliary.on
        m.vchamber_state = Auxiliary.on
        pf, st, lf = ring_pass(m, p, 1, turn_by_turn=true)
        find_orbit4(m)
        find_orbit6(m)
    end
end

end # module Track
