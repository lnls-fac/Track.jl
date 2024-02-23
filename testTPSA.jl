include("tpsaModule.jl")
using .TPSA

V, N = 1, 2

x = TPSA.Tpsa{V,N}(1, 1)

0 - x
