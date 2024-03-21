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

@exported_enum Status st_success st_particle_lost st_passmethod_not_defined st_passmethod_not_implemented st_inconsistent_dimensions st_uninitialized_memory st_findorbit_not_converged st_findorbit_one_turn_matrix_problem st_file_not_found st_file_not_opened st_kicktable_not_defined st_kicktable_out_of_range st_flat_file_error st_newton_not_converged st_not_implemented st_no_cavities_found

@exported_enum Plane no_plane plane_x plane_y plane_z plane_xy

# baremodule BoolState
#     export on, off, full
#     const on::Int = 0
#     const off::Int = 1
#     const full::Int = 2
# end

# baremodule PassMethod
#     export pm_identity_pass, pm_corrector_pass, pm_drift_pass, pm_matrix_pass,
#            pm_bnd_mpole_symplectic4_pass, pm_str_mpole_symplectic4_pass,
#            pm_cavity_pass, pm_kickmap_pass
#     const pm_identity_pass::Int = 0
#     const pm_corrector_pass::Int = 1
#     const pm_drift_pass::Int = 2
#     const pm_matrix_pass::Int = 3
#     const pm_bnd_mpole_symplectic4_pass::Int = 4
#     const pm_str_mpole_symplectic4_pass::Int = 5
#     const pm_cavity_pass::Int = 6
#     const pm_kickmap_pass::Int = 7
# end

# baremodule VChamberShape
#     export vchamber_rectangle, vchamber_rhombus, vchamber_ellipse
#     const vchamber_rectangle::Int = 0
#     const vchamber_rhombus::Int = 1
#     const vchamber_ellipse::Int = 2
# end

# baremodule Distributions
#     export dist_normal, dist_uniform
#     const dist_normal::Int = 0
#     const dist_uniform::Int = 1
# end

# baremodule Status
#     export st_success, st_particle_lost, st_passmethod_not_defined,
#            st_passmethod_not_implemented, st_inconsistent_dimensions,
#            st_uninitialized_memory, st_findorbit_not_converged,
#            st_findorbit_one_turn_matrix_problem, st_file_not_found,
#            st_file_not_opened, st_kicktable_not_defined,
#            st_kicktable_out_of_range, st_flat_file_error,
#            st_newton_not_converged, st_not_implemented,
#            st_no_cavities_found
#     const st_success::Int = 0
#     const st_particle_lost::Int = 1
#     const st_passmethod_not_defined::Int = 2
#     const st_passmethod_not_implemented::Int = 3
#     const st_inconsistent_dimensions::Int = 4
#     const st_uninitialized_memory::Int = 5
#     const st_findorbit_not_converged::Int = 6
#     const st_findorbit_one_turn_matrix_problem::Int = 7
#     const st_file_not_found::Int = 8
#     const st_file_not_opened::Int = 9
#     const st_kicktable_not_defined::Int = 10
#     const st_kicktable_out_of_range::Int = 11
#     const st_flat_file_error::Int = 12
#     const st_newton_not_converged::Int = 13
#     const st_not_implemented::Int = 14
#     const st_no_cavities_found::Int = 15
# end

# baremodule Plane
#     export no_plane, plane_x, plane_y, plane_z, plane_xy
#     const no_plane::Int = 0
#     const plane_x::Int = 1
#     const plane_y::Int = 2
#     const plane_z::Int = 3
#     const plane_xy::Int = 4
# end
