using ..Auxiliary: BoolState, off, on, full, st_findorbit_not_converged, st_success, st_findorbit_one_turn_matrix_problem, no_plane, pm_bnd_mpole_symplectic4_pass, st_no_cavities_found
using ..Tracking: line_pass, CGAMMA
using ..PosModule: Pos, Pos_get_max
using ..AcceleratorModule: Accelerator, find_cav_indices
using ..Elements: Element
using ..Constants: light_speed
using ..TPSA: Tpsa
using LinearAlgebra

export find_orbit4, find_orbit6

function find_orbit4(accelerator::Accelerator; energy_offset::Float64=0.0, fixed_point_guess::Pos{T} = Pos(0.0)) where T
    delta = 1e-9              # [m],[rad],[dE/E]
    tolerance = 2.22044604925e-14
    max_nr_iters = 50
    leng = length(accelerator.lattice)
    
    radsts::BoolState = accelerator.radiation_state
    if radsts == full
        accelerator.radiation_state = on
    end
    
    tpsa=false
    if isa(fixed_point_guess.rx, Tpsa)
        tpsa = true
    end
    
    fixed_point_guess.de += energy_offset

    co::Vector{Pos{T}} = fill(fixed_point_guess, 7)
    D::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 7)
    M::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 6)
    dco::Pos{T} = Pos(1.0, 1.0, 1.0, 1.0, 0.0, 0.0, tpsa=tpsa)
    
    D = matrix6_set_identity_posvec(D, delta=delta)

    nr_iter = 0
    while (Pos_get_max(dco) > tolerance) && (nr_iter <= max_nr_iters)
        co = co + D
        Ri = copy(co[7])
        status = st_success
        co2::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 7)
        for i in [1, 2, 3, 4, 7]
            pf, status, _ = line_pass(accelerator, co[i], [leng+1])
            co2[i] = copy(pf[1])
        end
        if status != st_success
            return Pos{T}[], st_findorbit_one_turn_matrix_problem
        end

        Rf = copy(co2[5])
        M[1] = (co2[1] - Rf) / delta
        M[2] = (co2[2] - Rf) / delta
        M[3] = (co2[3] - Rf) / delta
        M[4] = (co2[4] - Rf) / delta
        b = Rf - Ri
        M_1 = fill(Pos(0.0, tpsa=tpsa), 6)
        M_1 = matrix6_set_identity_posvec(M_1)
        M_1 = M_1 - M
        dco = linalg_solve4_posvec(M_1, b)
        co[7] = dco + Ri
        co[1] = copy(co[7]); co[2] = copy(co[7])
        co[3] = copy(co[7]); co[4] = copy(co[7])
        co[5] = copy(co[7]); co[6] = copy(co[7])
        nr_iter += 1
    end
    
    if nr_iter > max_nr_iters
        return Pos{T}[], st_findorbit_not_converged
    end
    
    closed_orbit, _, _ = line_pass(accelerator, co[7], "open")
    accelerator.radiation_state = radsts
    return closed_orbit, st_success
end

function find_orbit6(accelerator::Accelerator; fixed_point_guess::Pos{T} = Pos(0.0)) where T

    xy_delta = 1e-8              # [m],[rad],[dE/E]
    dp_delta = 1e-6
    delta = vcat([xy_delta*ones(Float64, 4)... , dp_delta*ones(Float64, 2)...])
    tolerance = 1e-15
    max_nr_iters = 50
    leng = length(accelerator.lattice)
    
    radsts::BoolState = accelerator.radiation_state
    if radsts == full
        accelerator.radiation_state = on
    end

    cav_indices::Vector{Int} = find_cav_indices(accelerator)
    if isempty(cav_indices)
        return Pos[], st_no_cavities_found
    end
    cavidx::Int = cav_indices[1]
    cav::Element = accelerator.lattice[cavidx]
    tpsa=false

    if isa(fixed_point_guess.rx, Tpsa)
        tpsa = true
    end

    # if true
    #     u0 = get_U0(accelerator)
    #     voltage = cav.voltage
    #     fixed_point_guess[6] = -accelerator.length/(2*pi*accelerator.harmonic_number) * asin(u0/voltage)
    # end

    frf::Float64 = cav.frequency
    longitudinal_fixed_point::Float64 = (accelerator.velocity / 1e8 * accelerator.harmonic_number / frf * 1e8) - accelerator.length
    println(stdout, "lng_fxdpt = $longitudinal_fixed_point")

    co::Vector{Pos{T}} = fill(fixed_point_guess, 7)
    co2::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 7)
    D::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 7)
    M::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), 6)
    dco::Pos{T} = Pos(1.0, tpsa=tpsa)
    theta::Pos{T} = Pos(0.0, tpsa=tpsa)
    theta.dl = longitudinal_fixed_point

    D = matrix6_set_identity_posvec(D, delta=delta)

    nr_iter = 0
    while (Pos_get_max(dco) > tolerance) && (nr_iter <= max_nr_iters)

        co = co + D

        Ri = copy(co[7]) # is *1e8 bigger

        status = st_success
        for i in [1, 2, 3, 4, 5, 6, 7]
            pf, status, _ = line_pass(accelerator, co[i], [leng+1])
            co2[i] = copy(pf[1])
        end

        if status != st_success
            return Pos{T}[], st_findorbit_one_turn_matrix_problem
        end

        Rf = copy(co2[7]) # is *1e8 bigger
        
        M[1] = (co2[1] - Rf) / (delta)
        
        M[2] = (co2[2] - Rf) / (delta)
        M[3] = (co2[3] - Rf) / (delta)
        M[4] = (co2[4] - Rf) / (delta)
        M[5] = (co2[5] - Rf) / (delta)
        M[6] = (co2[6] - Rf) / (delta)
        
        b = (Rf - Ri) - theta
        
        M_1 = fill(Pos(0.0, tpsa=tpsa), 6)
        M_1 = matrix6_set_identity_posvec(M_1)
        
        M_1 = M_1 - M # is *1e8 bigger
        
        dco = Pos(0.0, tpsa=tpsa)
        dco = linalg_solve6_posvec(M_1, b)
        
        co[7] = dco + Ri
        co[1] = copy(co[7]); co[2] = copy(co[7])
        co[3] = copy(co[7]); co[4] = copy(co[7])
        co[5] = copy(co[7]); co[6] = copy(co[7])
        nr_iter += 1
    end
    
    if nr_iter > max_nr_iters
        return Pos{T}[], st_findorbit_not_converged
    end
    
    closed_orbit, _, _ = line_pass(accelerator, co[7], "open")
    accelerator.radiation_state = radsts
    return closed_orbit, st_success
end

function matrix6_set_identity_posvec(D::Vector{Pos{T}}; delta::Union{Float64, Vector{Float64}}=1.0) where T
    tpsa = false

    if isa(D[1].rx, Tpsa)
        tpsa=true
    end
    if isa(delta, Float64)
        delta = ones(Float64, 6) * delta
    end
    M::Vector{Pos{T}} = fill(Pos(0.0, tpsa=tpsa), length(D))
    for i in 1:6
        d = Pos(0.0, tpsa=tpsa)
        d[i] = delta[i]
        M[i] = d
    end
    return M
end

function linalg_solve4_posvec(A::Vector{Pos{T}}, B::Pos{T}) where T

    tpsa = false
    if isa(B.rx, Tpsa)
        tpsa=true
    end

    m = Matrix{T}(undef, 4, 4)   # Create a 4x4 matrix
    b = Vector{T}(undef, 4)      # Create a vector of size 4
    x = Vector{T}(undef, 4)      # Create a solution vector of size 4
    
    # Populate the matrix and vector
    for i in 1:4
        m[:, i] .= [A[i][j] for j in 1:4]
        b[i] = B[i]
    end

    # Solve the system using LU decomposition
    x .= m \ b
    
    # Create the solution vector
    X = Pos(0.0, tpsa=tpsa)
    X[1] = x[1]
    X[2] = x[2]
    X[3] = x[3]
    X[4] = x[4]
    
    return X
end

function linalg_solve6_posvec(A::Vector{Pos{T}}, B::Pos{T}) where T

    tpsa = false
    if isa(B.rx, Tpsa)
        tpsa=true
    end

    m = Matrix{T}(undef, 6, 6)   # Create a 6x6 matrix
    b = Vector{T}(undef, 6)      # Create a vector of size 6
    x = Vector{T}(undef, 6)      # Create a solution vector of size 6
    
    # Populate the matrix and vector
    for i in 1:6
        m[:, i] .= [A[i][j] for j in 1:6]
        b[i] = B[i]
    end

    # Solve the system using LU decomposition
    x .= m \ b
    
    # Create the solution vector
    X = Pos(0.0, tpsa=tpsa)
    X[1] = x[1]
    X[2] = x[2]
    X[3] = x[3]
    X[4] = x[4]
    X[5] = x[5]
    X[6] = x[6]
    
    return X
end

function get_U0(accelerator::Accelerator)
    theta::Vector{Float64} = Float64[]
    leng::Vector{Float64} = Float64[]
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