using ..Auxiliary: BoolState, off, on, full, st_findorbit_not_converged, st_success, st_findorbit_one_turn_matrix_problem, no_plane, pm_bnd_mpole_symplectic4_pass
using ..Tracking: line_pass, CGAMMA
using ..PosModule: Pos, Pos_get_max
using ..AcceleratorModule: Accelerator, find_cav_indices
using ..Elements: Element
using ..Constants: light_speed
using LinearAlgebra

export find_orbit4, find_orbit6

function find_orbit4(accelerator::Accelerator; fixed_point_guess::Pos{Float64} = Pos(0.0), element_offset::Int=1)
    delta = 1e-9              # [m],[rad],[dE/E]
    tolerance = 2.22044604925e-14
    max_nr_iters = 50
    leng = length(accelerator.lattice)
    
    radsts::BoolState = accelerator.radiation_state
    if radsts == full
        accelerator.radiation_state = on
    end
    
    co::Vector{Pos{Float64}} = fill(fixed_point_guess, 7)
    D::Vector{Pos{Float64}} = fill(Pos(0.0), 7)
    M::Vector{Pos{Float64}} = fill(Pos(0.0), 6)
    dco::Pos{Float64} = Pos(1.0, 1.0, 1.0, 1.0, 0.0, 0.0)
    theta::Pos{Float64} = Pos(0.0)
    
    D = matrix6_set_identity_posvec(D, delta=delta)

    nr_iter = 0
    while (Pos_get_max(dco) > tolerance) && (nr_iter <= max_nr_iters)
        co = co + D
        Ri = co[7]
        status = st_success
        co2::Vector{Pos{Float64}} = fill(Pos(0.0), 7)
        for i in [1, 2, 3, 4, 7]
            pf, status, _ = line_pass(accelerator, co[i], [leng+1])
            co2[i] = copy(pf[1])
        end
        if status != st_success
            return Pos{Float64}[], st_findorbit_one_turn_matrix_problem
        end

        Rf = co2[5]
        M[1] = (co2[1] - Rf) / delta
        M[2] = (co2[2] - Rf) / delta
        M[3] = (co2[3] - Rf) / delta
        M[4] = (co2[4] - Rf) / delta
        b = Rf - Ri
        M_1 = fill(Pos(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 6)
        matrix6_set_identity_posvec(M_1)
        M_1 = M_1 - M
        dco = linalg_solve4_posvec(M_1, b)
        co[7] = dco + Ri
        co[1] = co[7]; co[2] = co[7]
        co[3] = co[7]; co[4] = co[7]
        co[5] = co[7]; co[6] = co[7]
        nr_iter += 1
    end
    
    if nr_iter > max_nr_iters
        return Pos{Float64}[], st_findorbit_not_converged
    end
    
    closed_orbit, _, _ = line_pass(accelerator, co[7], "open", element_offset=element_offset)
    accelerator.radiation_state = radsts
    return closed_orbit, st_success
end

function find_orbit6(accelerator::Accelerator; fixed_point_guess::Pos{Float64} = Pos(0.0), element_offset::Int=1)

    delta = 1e-9              # [m],[rad],[dE/E]
    tolerance = 2.22044604925e-14
    max_nr_iters = 50
    leng = length(accelerator.lattice)
    
    radsts::BoolState = accelerator.radiation_state
    if radsts == full
        accelerator.radiation_state = on
    end

    cav::Element = accelerator.lattice[find_cav_indices(accelerator)[1]]

    if true
        u0 = get_U0(accelerator)
        voltage = cav.voltage
        fixed_point_guess[6] = -accelerator.length/(2*pi*accelerator.harmonic_number) * asin(u0/voltage)
    end

    frf::Float64 = cav.frequency
    longitudinal_fixed_point::Float64 = (accelerator.velocity * accelerator.harmonic_number / frf) - accelerator.length

    co::Vector{Pos{Float64}} = fill(fixed_point_guess, 7)
    D::Vector{Pos{Float64}} = fill(Pos(0.0), 7)
    M::Vector{Pos{Float64}} = fill(Pos(0.0), 6)
    dco::Pos{Float64} = Pos(1.0, 1.0, 1.0, 1.0, 0.0, 0.0)
    theta::Pos{Float64} = Pos(0.0)
    theta.dl = longitudinal_fixed_point

    D = matrix6_set_identity_posvec(D, delta=delta)

    nr_iter = 0
    while (Pos_get_max(dco) > tolerance) && (nr_iter <= max_nr_iters)
        co = co + D
        Ri = co[7]
        status = st_success
        co2::Vector{Pos{Float64}} = fill(Pos(0.0), 7)
        for i in [1, 2, 3, 4, 5, 6, 7]
            pf, status, _ = line_pass(accelerator, co[i], [leng+1], turn_number=1)
            co2[i] = copy(pf[1])
        end
        if status != st_success
            return Pos{Float64}[], st_findorbit_one_turn_matrix_problem
        end

        Rf = co2[7]
        M[1] = (co2[1] - Rf) / delta
        M[2] = (co2[2] - Rf) / delta
        M[3] = (co2[3] - Rf) / delta
        M[4] = (co2[4] - Rf) / delta
        M[5] = (co2[5] - Rf) / delta
        M[6] = (co2[6] - Rf) / delta
        b = Rf - Ri - theta
        M_1 = fill(Pos(0.0), 6)
        matrix6_set_identity_posvec(M_1)
        M_1 = M_1 - M
        dco = linalg_solve4_posvec(M_1, b)
        co[7] = dco + Ri
        co[1] = co[7]; co[2] = co[7]
        co[3] = co[7]; co[4] = co[7]
        co[5] = co[7]; co[6] = co[7]
        nr_iter += 1
    end
    
    if nr_iter > max_nr_iters
        return Pos{Float64}[], st_findorbit_not_converged
    end
    
    closed_orbit, _, _ = line_pass(accelerator, co[7], "open", element_offset=element_offset)
    accelerator.radiation_state = radsts
    return closed_orbit, st_success
end

function matrix6_set_identity_posvec(D::Vector{Pos{Float64}}; delta::Float64=1.0)
    M::Vector{Pos{Float64}} = fill(Pos(0.0), length(D))
    for i in 1:1:6
        d = Pos(0.0)
        d[i] = delta
        M[i] = d
    end
    return M
end

function linalg_solve4_posvec(A::Vector{Pos{Float64}}, B::Pos{Float64})
    m = Matrix{Float64}(undef, 4, 4)   # Create a 4x4 matrix
    b = Vector{Float64}(undef, 4)      # Create a vector of size 4
    x = Vector{Float64}(undef, 4)      # Create a solution vector of size 4
    
    # Populate the matrix and vector
    for i in 1:4
        m[:, i] .= [A[j][i] for j in 1:4]
        b[i] = B[i]
    end

    # Solve the system using LU decomposition
    x .= m \ b
    
    # Create the solution vector
    X = Pos(x[1], x[2], x[3], x[4], 0.0, 0.0)
    
    return X
end

function linalg_solve6_posvec(A::Vector{Pos{Float64}}, B::Pos{Float64})
    m = Matrix{Float64}(undef, 6, 6)   # Create a 4x4 matrix
    b = Vector{Float64}(undef, 6)      # Create a vector of size 4
    x = Vector{Float64}(undef, 6)      # Create a solution vector of size 4
    
    # Populate the matrix and vector
    for i in 1:6
        m[:, i] .= [A[j][i] for j in 1:6]
        b[i] = B[i]
    end

    # Solve the system using LU decomposition
    x .= m \ b
    
    # Create the solution vector
    X = Pos(x[1], x[2], x[3], x[4], x[5], x[6])
    
    return X
end

function get_U0(accelerator::Accelerator)
    theta = []
    leng = []
    for e in accelerator.lattice
        if e.pass_method == pm_bnd_mpole_symplectic4_pass
            push!(theta, e.angle)
            push!(leng, e.length)
        end
    end
    coef = CGAMMA*1e9/2/pi * (accelerator.energy/1e9)^4
    u = sum(abs.(theta.*theta./leng)) * coef
    return u
end