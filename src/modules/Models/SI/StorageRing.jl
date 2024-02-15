"""SIRIUS StorageRing Model module."""
module StorageRing

    const default_optics_mode::String = "S05.01"

    include("model.jl")

    function create_accelerator(;optics_mode::String=default_optics_mode, simplified::Bool=false, ids=[])
        return _create_lattice(optics_mode, simplified, ids)
    end
    
    export create_accelerator
end # module StorageRing