# Auxiliary.jl

macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end

@exported_enum BoolState off on full

@exported_enum PassMethod pm_identity_pass pm_corrector_pass pm_drift_pass pm_matrix_pass pm_bnd_mpole_symplectic4_pass pm_str_mpole_symplectic4_pass pm_cavity_pass pm_kickmap_pass

@exported_enum VChamberShape vchamber_rectangle vchamber_rhombus vchamber_ellipse

@exported_enum Distributions dist_normal dist_uniform

@exported_enum Status st_success st_particle_lost st_passmethod_not_defined st_passmethod_not_implemented st_inconsistent_dimensions st_uninitialized_memory st_findorbit_not_converged st_findorbit_one_turn_matrix_problem st_file_not_found st_file_not_opened st_kicktable_not_defined st_kicktable_out_of_range st_flat_file_error st_newton_not_converged st_not_implemented

@exported_enum Plane no_plane plane_x plane_y plane_z plane_xy