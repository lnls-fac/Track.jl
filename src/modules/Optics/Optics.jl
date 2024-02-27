using ..MatrixModule: find_matrix66
using ..Orbit: find_orbit6
using ..AcceleratorModule: Accelerator

export calc_twiss

struct Twiss 
    betax
    betay
    alphax
    alphay
    mux
    muy
end

function calc_twiss(accelerator::Accelerator)

    fp, _ = find_orbit6(accelerator)
    fp = fp[end]

    m66, tm = find_matrix66(accelerator, fp, indices="all")

    sin_mux0 = sign(m66[0+1,1+1]) * sqrt(-m66[0+1,1+1]*m66[1+1,0+1] - (m66[0+1,0+1]-m66[1+1,1+1])^2/4.0)
    sin_muy0 = sign(m66[2+1,3+1]) * sqrt(-m66[2+1,3+1]*m66[3+1,2+1] - (m66[2+1,2+1]-m66[3+1,3+1])^2/4.0)
    alphax0 = (m66[0+1,0+1]-m66[1+1,1+1])/2.0/sin_mux0
    alphay0 = (m66[2+1,2+1]-m66[3+1,3+1])/2.0/sin_muy0
    betax0  =  m66[0+1,1+1]/sin_mux0
    betay0  =  m66[2+1,3+1]/sin_muy0

    betax = [betax0]
    betay = [betay0]
    alphax = [alphax0]
    alphay = [alphay0]
    mux = [asin(sin_mux0)]
    muy = [asin(sin_muy0)]

    for m in tm[1:end-1]
        _betax = ((m[0+1, 1+0] * betax0 - m[0+1, 1+1] * alphax0)^2 + (m[0+1, 1+1])^2) / betax0;
        _betay = ((m[2+1, 1+2] * betay0 - m[2+1, 1+3] * alphay0)^2 + (m[2+1, 1+3])^2) / betay0;
        _alphax = -((m[0+1, 1+0] * betax0 - m[0+1, 1+1] * alphax0) * (m[1+1, 1+0] * betax0 - m[1+1, 1+1] * alphax0) + m[0+1, 1+1]*m[1+1, 1+1])/betax0;
        _alphay = -((m[2+1, 1+2] * betay0 - m[2+1, 1+3] * alphay0) * (m[3+1, 1+2] * betay0 - m[3+1, 1+3] * alphay0) + m[2+1, 1+3]*m[3+1, 1+3])/betay0;
        _mux = atan(m[0+1, 1+1]/(m[0+1, 1+0] * betax0 - m[0+1, 1+1] * alphax0));
        _muy = atan(m[2+1, 1+3]/(m[2+1, 1+2] * betay0 - m[2+1, 1+3] * alphay0));

        push!(betax, _betax)
        push!(betay, _betay)
        push!(alphax, _alphax)
        push!(alphay, _alphay)
        push!(mux, _mux)
        push!(muy, _muy)
    end

    twiss = Twiss(betax, betay, alphax, alphay, mux, muy)

    return twiss

end