"""Track module."""
module Track

include("modules/Constants/constantsModule.jl")
include("modules/Auxiliary/auxiliaryModule.jl")
include("modules/TPSA/tpsaModule.jl")
include("modules/Pos/posModule.jl")
include("modules/Elements/elementsModule.jl")
include("modules/Accelerator/acceleratorModule.jl")
include("modules/Tracking/trackingModule.jl")
include("modules/Orbit/orbitModule.jl")
include("modules/Matrix/matrixModule.jl")
include("modules/Optics/opticsModule.jl")

using .PosModule: Pos
using .AcceleratorModule: Accelerator, find_spos, find_indices

# using PrecompileTools

# @setup_workload begin
#     @compile_workload begin
#         F = Elements.quadrupole("F", 0.01, +1.0)
#         D = Elements.quadrupole("D", 0.01, -1.0)
#         O = Elements.drift("O", 1.0)
#         cell = Elements.Element[F, O, D, O]
#         lattice = Elements.Element[]
#         for i=1:10
#             append!(lattice, copy(cell))
#         end
#         acc = Accelerator(3e9)
#         acc.lattice = lattice
#         acc.radiation_state = 1
#         acc.vchamber_state = true
#         acc.cavity_state = 1

#         twiss = Optics.calc_twiss(acc)
#     end
# end

end # module Track
