"""Track module."""
module Track

include("modules/Constants/constantsModule.jl")
include("modules/Auxiliary/auxiliaryModule.jl")
include("modules/Pos/posModule.jl")
include("modules/Elements/elementsModule.jl")
include("modules/Accelerator/acceleratorModule.jl")
include("modules/Tracking/trackingModule.jl")
include("modules/Orbit/orbitModule.jl")
include("modules/Models/modelsModule.jl")

using .PosModule: Pos
using .AcceleratorModule: Accelerator, find_spos, find_indices
using .Tracking: element_pass, line_pass, ring_pass
using .Orbit: find_orbit4, find_orbit6
using .Models: StorageRing

export Constants, 
        Auxiliary, 
        Pos, 
        Elements, 
        Accelerator, find_spos, find_indices, 
        element_pass, line_pass, ring_pass,
        find_orbit4, find_orbit6
        StorageRing

using PrecompileTools

@setup_workload begin
    @compile_workload begin
        v::Vector{Float64} = rand(Float64, 6) * 1e-6
        p = Pos(v)
        m = StorageRing.create_accelerator()
        m.radiation_state = Auxiliary.full
        m.cavity_state = Auxiliary.on
        m.vchamber_state = Auxiliary.on
        pf, st, lf = ring_pass(m, p, 1)
        find_orbit4(m)
        find_orbit6(m)
    end
end

end # module Track
