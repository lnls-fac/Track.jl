using ..AcceleratorModule: Accelerator
using ..PosModule: Pos
using ..Tracking: element_pass
using ..Auxiliary
using ..TPSA: Tpsa

function find_matrix66(accelerator::Accelerator, fixed_point::Pos{T}=Pos(0); indices::String="m66") where {T}
    
    if lowercase(indices) == "m66"
        return_tm_flag = false
    elseif lowercase(indices) == "all"
        return_tm_flag = true
    else
        error("invalid indices {String}: should be: \"m66\" for only the One-Turn Matrix or \"all\" for the cumulative Matrix in all elements.")
    end

    radsts::BoolState = accelerator.radiation_state
    if radsts == full
        accelerator.radiation_state = on
    end

    map::Pos{Tpsa{6, 1, Float64}} = Pos(fixed_point[:], tpsa=true) 
    
    if return_tm_flag
        tm::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}[]
    end
 
    for elem in accelerator.lattice
        element_pass(elem, map, accelerator)
        if return_tm_flag
            mat = extract_m66(map)
            push!(tm, mat)
        end
    end

    m66 = extract_m66(map)

    accelerator.radiation_state = radsts

    if return_tm_flag
        return m66, tm
    else
        return m66
    end
end


function extract_m66(p::Pos{Tpsa{6, 1, Float64}})
    return [
        p.rx[1] p.rx[2] p.rx[3] p.rx[4] p.rx[5] p.rx[6];
        p.px[1] p.px[2] p.px[3] p.px[4] p.px[5] p.px[6];
        p.ry[1] p.ry[2] p.ry[3] p.ry[4] p.ry[5] p.ry[6];
        p.py[1] p.py[2] p.py[3] p.py[4] p.py[5] p.py[6];
        p.de[1] p.de[2] p.de[3] p.de[4] p.de[5] p.de[6];
        p.dl[1] p.dl[2] p.dl[3] p.dl[4] p.dl[5] p.dl[6]
    ]    
end