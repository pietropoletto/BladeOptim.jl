module CCBladeAnalysis

using CCBlade
using Statistics

export analyze_performance, analyze_loads, analyze_tsr_sweep, analyze_cp_robust
export analyze_performance_multipolar, analyze_loads_multipolar
export analyze_tsr_sweep_multipolar, analyze_cp_robust_multipolar

# =============================================================================
#  ORIGINAL: single polar for all sections
# =============================================================================

function analyze_performance(turbine, af_file::String;
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20)
    tsrvec = range(tsr_range[1], tsr_range[2], length=ntsr)
    cpvec = zeros(ntsr); ctvec = zeros(ntsr)
    airfoils = AlphaAF(af_file, radians=false)
    sections = Section.(turbine.r, turbine.chord, turbine.twist, Ref(airfoils))
    for i = 1:ntsr
        Omega = turbine.V_inf * tsrvec[i] / turbine.rotorR
        ops = windturbine_op.(turbine.V_inf, Omega, turbine.pitch, turbine.r,
            turbine.precone, turbine.yaw, turbine.tilt,
            turbine.azimuth', turbine.hubHt, turbine.shearExp, turbine.rho)
        outs = solve.(Ref(turbine.rotor), sections, ops)
        T, Q = thrusttorque(turbine.rotor, sections, outs)
        cpvec[i], ctvec[i], _ = nondim(T, Q, turbine.V_inf, Omega, turbine.rho, turbine.rotor, "windturbine")
    end
    idx_max = argmax(cpvec)
    return tsrvec, cpvec, ctvec, tsrvec[idx_max], cpvec[idx_max]
end

function analyze_loads(turbine, af_file::String)
    airfoils = AlphaAF(af_file, radians=false)
    sections = Section.(turbine.r, turbine.chord, turbine.twist, Ref(airfoils))
    op = windturbine_op.(turbine.V_inf, turbine.omega, turbine.pitch, turbine.r,
        turbine.precone, turbine.yaw, turbine.tilt, 0.0,
        turbine.hubHt, turbine.shearExp, turbine.rho)
    out = solve.(Ref(turbine.rotor), sections, op)
    T, Q = thrusttorque(turbine.rotor, sections, out)
    return T, Q, Q*turbine.omega, getfield.(out, :Np), getfield.(out, :Tp)
end

function analyze_tsr_sweep(turbine, af_file::String, tsr_target::Float64;
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20)
    tsrvec, cpvec, ctvec, tsr_opt, cp_max = analyze_performance(turbine, af_file; tsr_range, ntsr)
    idx = argmin(abs.(tsrvec .- tsr_target))
    return tsrvec, cpvec, ctvec, cpvec[idx], tsr_opt, cp_max
end

function analyze_cp_robust(turbine, af_file::String;
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20, delta_tsr::Float64 = 2.0)
    tsrvec, cpvec, ctvec, tsr_opt, cp_max = analyze_performance(turbine, af_file; tsr_range, ntsr)
    mask = (tsrvec .>= tsr_opt - delta_tsr) .& (tsrvec .<= tsr_opt + delta_tsr)
    cp_robust = sum(mask) == 0 ? cp_max : mean(cpvec[mask])
    return tsrvec, cpvec, ctvec, tsr_opt, cp_max, cp_robust
end

# =============================================================================
#  NEW: multi-polar — one AlphaAF per section
#  af_files : Vector{String} of length n_sections
# =============================================================================

function _build_multipolar_sections(turbine, af_files::Vector{String})
    n_sec = length(turbine.r)
    @assert length(af_files) == n_sec "Need one af_file per section (got $(length(af_files)) for $n_sec sections)"
    af_vec = [AlphaAF(f, radians=false) for f in af_files]
    return [Section(turbine.r[i], turbine.chord[i], turbine.twist[i], af_vec[i]) for i in 1:n_sec]
end

function analyze_performance_multipolar(turbine, af_files::Vector{String};
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20)
    sections = _build_multipolar_sections(turbine, af_files)
    tsrvec = range(tsr_range[1], tsr_range[2], length=ntsr)
    cpvec = zeros(ntsr); ctvec = zeros(ntsr)
    for i = 1:ntsr
        Omega = turbine.V_inf * tsrvec[i] / turbine.rotorR
        ops = windturbine_op.(turbine.V_inf, Omega, turbine.pitch, turbine.r,
            turbine.precone, turbine.yaw, turbine.tilt,
            turbine.azimuth', turbine.hubHt, turbine.shearExp, turbine.rho)
        outs = solve.(Ref(turbine.rotor), sections, ops)
        T, Q = thrusttorque(turbine.rotor, sections, outs)
        cpvec[i], ctvec[i], _ = nondim(T, Q, turbine.V_inf, Omega, turbine.rho, turbine.rotor, "windturbine")
    end
    idx_max = argmax(cpvec)
    return tsrvec, cpvec, ctvec, tsrvec[idx_max], cpvec[idx_max]
end

function analyze_loads_multipolar(turbine, af_files::Vector{String})
    sections = _build_multipolar_sections(turbine, af_files)
    op = windturbine_op.(turbine.V_inf, turbine.omega, turbine.pitch, turbine.r,
        turbine.precone, turbine.yaw, turbine.tilt, 0.0,
        turbine.hubHt, turbine.shearExp, turbine.rho)
    out = solve.(Ref(turbine.rotor), sections, op)
    T, Q = thrusttorque(turbine.rotor, sections, out)
    return T, Q, Q*turbine.omega, getfield.(out, :Np), getfield.(out, :Tp)
end

function analyze_tsr_sweep_multipolar(turbine, af_files::Vector{String}, tsr_target::Float64;
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20)
    tsrvec, cpvec, ctvec, tsr_opt, cp_max = analyze_performance_multipolar(turbine, af_files; tsr_range, ntsr)
    idx = argmin(abs.(tsrvec .- tsr_target))
    return tsrvec, cpvec, ctvec, cpvec[idx], tsr_opt, cp_max
end

function analyze_cp_robust_multipolar(turbine, af_files::Vector{String};
        tsr_range::Tuple{Float64, Float64} = (1.0, 10.0), ntsr::Int = 20, delta_tsr::Float64 = 2.0)
    tsrvec, cpvec, ctvec, tsr_opt, cp_max = analyze_performance_multipolar(turbine, af_files; tsr_range, ntsr)
    mask = (tsrvec .>= tsr_opt - delta_tsr) .& (tsrvec .<= tsr_opt + delta_tsr)
    cp_robust = sum(mask) == 0 ? cp_max : mean(cpvec[mask])
    return tsrvec, cpvec, ctvec, tsr_opt, cp_max, cp_robust
end

end # module
