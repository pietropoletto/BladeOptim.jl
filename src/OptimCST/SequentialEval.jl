module SequentialEval

using DataFrames
using Random
using ..Config
using ..CST
using ..XfoilRunner
using ..AlphaOpt
using ..Viterna
using ..Conversions
using ..CCBladeAnalysis
using ..ObjFunction: OBJECTIVE, _interp_polar, soft_stall_penalty, CM_CAP_TIP, CM_CAP_ROOT

const _EVAL_SPACER = "  -     -     -     -     -     -     -     -"

export run_xfoil_single, evaluate_from_polars
export blend_weight, cl_cd_at_alpha
export check_single_airfoil, has_crossover
export sample_feasible_single, lhs_candidates, param_to_coords

# =============================================================================
# XFOIL SINGLE AIRFOIL
# =============================================================================

"""
    run_xfoil_single(wu, wl, cfg, dt; Re=nothing) → DataFrame | nothing

Generate CST coordinates and run XFOIL.
"""
function run_xfoil_single(wu, wl, cfg, dt; Re::Union{Float64,Nothing}=nothing)
    coord_raw = CST.CST_airfoil(wu, wl, cfg.dz, cfg.N)
    coord     = CST.airfoil_thickness(coord_raw, dt)
    xp        = Config.create_xfoil_params()
    re_val    = Re === nothing ? xp.Reynolds : Re
    df = XfoilRunner.run_xfoil(
        coord[:, 1], coord[:, 2];
        Re=re_val, Mach=xp.Mach,
        alpha_start=xp.alpha_start, alpha_end=xp.alpha_end,
        alpha_step=xp.alpha_step, iter=xp.Iterations,
        ncrit=xp.ncrit, show_airfoil=true)
    return df
end

# =============================================================================
# HELPERS
# =============================================================================

function cl_cd_at_alpha(df::DataFrames.DataFrame, alpha_target::Float64)
    idx = argmin(abs.(df.alpha .- alpha_target))
    cm_val = hasproperty(df, :cm) ? df.cm[idx] : NaN
    return df.cl[idx], df.cd[idx], cm_val
end

function blend_weight(eta::Real; eta_root::Real=0.40, eta_tip::Real=0.90)
    if eta <= eta_root
        return 0.0
    elseif eta >= eta_tip
        return 1.0
    else
        return (eta - eta_root) / (eta_tip - eta_root)
    end
end

function param_to_coords(wu, wl, cfg, dt)
    raw = CST.CST_airfoil(wu, wl, cfg.dz, cfg.N)
    return CST.airfoil_thickness(raw, dt)
end

# =============================================================================
# CORE: evaluate blade from two polar DataFrames
# =============================================================================

"""
    evaluate_from_polars(df_root, df_tip; Re_root, Re_tip, variable_phase, ...)

Evaluate turbine Cp from two polar DataFrames (root and tip).
`variable_phase` controls the Cm hard cap: `:tip`, `:root`, or `:both`.
"""
function evaluate_from_polars(df_root::DataFrames.DataFrame,
                              df_tip::DataFrames.DataFrame;
                              alpha_root_forced::Float64 = NaN,
                              label::String = "",
                              objective::Symbol = OBJECTIVE,
                              Re_root::Union{Float64,Nothing} = nothing,
                              Re_tip::Union{Float64,Nothing}  = nothing,
                              variable_phase::Symbol = :both)

    Rhub = Config.TURBINE.Rhub
    Rtip = Config.TURBINE.Rtip
    n_sections = Config.TURBINE.n_sections

    alpha_opt_tip, cl_opt_tip, cd_opt_tip, cm_opt_tip = AlphaOpt.find_optimal_alpha(df_tip)
    clcd_tip = cl_opt_tip / cd_opt_tip

    alpha_opt_root_nat, cl_opt_root_nat, cd_opt_root_nat, cm_opt_root = AlphaOpt.find_optimal_alpha(df_root)

    if !isnan(alpha_root_forced)
        cl_rf, cd_rf, cm_rf = cl_cd_at_alpha(df_root, alpha_root_forced)
        clcd_root = cl_rf / cd_rf
    else
        clcd_root = cl_opt_root_nat / cd_opt_root_nat
    end

    # ── Cm HARD CAP ─────────────────────────────────────────────────────
    _check_tip  = variable_phase in (:tip, :both)
    _check_root = variable_phase in (:root, :both)

    if _check_tip && !isinf(CM_CAP_TIP) && cm_opt_tip < CM_CAP_TIP
        println("$label  ✘ Cm_tip=$(round(cm_opt_tip, digits=4)) < cap $(CM_CAP_TIP) → REJECTED")
        println(_EVAL_SPACER)
        return -1e6, alpha_opt_tip, clcd_tip, alpha_opt_root_nat, clcd_root, cm_opt_tip, cm_opt_root
    end
    if _check_root && !isinf(CM_CAP_ROOT) && cm_opt_root < CM_CAP_ROOT
        println("$label  ✘ Cm_root=$(round(cm_opt_root, digits=4)) < cap $(CM_CAP_ROOT) → REJECTED")
        println(_EVAL_SPACER)
        return -1e6, alpha_opt_tip, clcd_tip, alpha_opt_root_nat, clcd_root, cm_opt_tip, cm_opt_root
    end

    # ── Soft-stall penalty ──────────────────────────────────────────────
    pen_tip   = soft_stall_penalty(df_tip)
    pen_root  = soft_stall_penalty(df_root)
    pen_max   = max(pen_tip, pen_root)
    stall_mult = 1.0 - pen_max
    if pen_max > 0.0
        println("$label  ⚠ Soft-stall: pen_tip=$(round(pen_tip, digits=3))  pen_root=$(round(pen_root, digits=3))  → obj × $(round(stall_mult, digits=3))")
    end

    # :clcd — quick return
    if objective == :clcd
        clcd_pen = clcd_tip * stall_mult
        println("$label  Cl/Cd_tip=$(round(clcd_tip, digits=1))  Cl/Cd_root=$(round(clcd_root, digits=1))  Cm_tip=$(round(cm_opt_tip, digits=4))  Cm_root=$(round(cm_opt_root, digits=4))$(pen_max > 0 ? "  [penalized=$(round(clcd_pen, digits=1))]" : "")")
        println(_EVAL_SPACER)
        return clcd_pen, alpha_opt_tip, clcd_tip, alpha_opt_root_nat, clcd_root, cm_opt_tip, cm_opt_root
    end

    # Build turbine geometry
    turbine = Config.turbine_parameters(alpha_opt_tip, cl_opt_tip, cd_opt_tip, n_sections;
                                            alpha_opt_root = alpha_opt_root_nat)

    # Per-section blended polars → Viterna → CCBlade
    all_alphas_raw = sort(unique(vcat(df_root.alpha, df_tip.alpha)))
    a_min = max(minimum(df_root.alpha), minimum(df_tip.alpha))
    a_max = min(maximum(df_root.alpha), maximum(df_tip.alpha))
    all_alphas = filter(a -> a_min <= a <= a_max, all_alphas_raw)

    cl_root_g, cd_root_g = _interp_polar(df_root, all_alphas)
    cl_tip_g,  cd_tip_g  = _interp_polar(df_tip,  all_alphas)

    xp       = Config.create_xfoil_params()
    r_vec    = turbine.r
    af_files = String[]

    use_re_interp = (Re_root !== nothing && Re_tip !== nothing)

    for j in 1:n_sections
        eta_j = clamp((r_vec[j] - Rhub) / (Rtip - Rhub), 0.0, 1.0)
        w_tip = blend_weight(eta_j)

        cl_j = (1.0 - w_tip) .* cl_root_g .+ w_tip .* cl_tip_g
        cd_j = (1.0 - w_tip) .* cd_root_g .+ w_tip .* cd_tip_g

        df_j = DataFrames.DataFrame(alpha = all_alphas, cl = cl_j, cd = cd_j)

        df_j_ext, _ = Viterna.apply_viterna_extrapolation(
            df_j, turbine.r, turbine.chord, turbine.Rtip)

        re_j = use_re_interp ?
            Re_root + clamp(eta_j, 0.0, 1.0) * (Re_tip - Re_root) :
            xp.Reynolds

        af_path = tempname() * ".dat"
        Conversions.write_ccblade_airfoil(df_j_ext, af_path, re_j, xp.Mach)
        push!(af_files, af_path)
    end

    if objective == :cp_robust
        _, _, _, tsr_opt, cp_max, cp_robust = CCBladeAnalysis.analyze_cp_robust_multipolar(
            turbine, af_files)
        cp_robust_pen = cp_robust * stall_mult
        println("$label  Cp_robust=$(round(cp_robust, digits=4))  Cp_max=$(round(cp_max, digits=4))  λ=$(round(tsr_opt, digits=2))  Cm_tip=$(round(cm_opt_tip, digits=4))  Cm_root=$(round(cm_opt_root, digits=4))  (Cl/Cd_tip=$(round(clcd_tip, digits=1)), Cl/Cd_root=$(round(clcd_root, digits=1)))$(pen_max > 0 ? "  [penalized=$(round(cp_robust_pen, digits=4))]" : "")")
        println(_EVAL_SPACER)
        return cp_robust_pen, alpha_opt_tip, clcd_tip, alpha_opt_root_nat, clcd_root, cm_opt_tip, cm_opt_root
    end

    _, _, _, _, cp_max = CCBladeAnalysis.analyze_performance_multipolar(
        turbine, af_files)
    cp_max_pen = cp_max * stall_mult

    println("$label  Cp_max=$(round(cp_max, digits=4))  Cm_tip=$(round(cm_opt_tip, digits=4))  Cm_root=$(round(cm_opt_root, digits=4))  (Cl/Cd_tip=$(round(clcd_tip, digits=1)), Cl/Cd_root=$(round(clcd_root, digits=1)))$(pen_max > 0 ? "  [penalized=$(round(cp_max_pen, digits=4))]" : "")")
    println(_EVAL_SPACER)
    return cp_max_pen, alpha_opt_tip, clcd_tip, alpha_opt_root_nat, clcd_root, cm_opt_tip, cm_opt_root
end

# =============================================================================
# GEOMETRY QUALITY CHECKS
# =============================================================================

function _split_upper_lower(coord)
    le_idx = argmin(coord[:, 1])
    return coord[1:le_idx, :], coord[le_idx:end, :]
end

function _local_thickness_at(coord, xc_target; tol=0.03)
    upper, lower = _split_upper_lower(coord)
    yu_pts = upper[abs.(upper[:, 1] .- xc_target) .< tol, 2]
    yl_pts = lower[abs.(lower[:, 1] .- xc_target) .< tol, 2]
    (isempty(yu_pts) || isempty(yl_pts)) && return NaN
    return sum(yu_pts)/length(yu_pts) - sum(yl_pts)/length(yl_pts)
end

function _min_aft_thickness(coord; xc_start=0.50, xc_end=0.92, n_probe=10)
    min_t = Inf
    for xc in range(xc_start, xc_end, length=n_probe)
        t = _local_thickness_at(coord, xc)
        if !isnan(t) && t < min_t; min_t = t; end
    end
    return min_t
end

function has_crossover(coord; xc_start=0.30)
    for xc in range(xc_start, 0.98, length=15)
        t = _local_thickness_at(coord, xc; tol=0.03)
        if !isnan(t) && t < -0.001; return true; end
    end
    return false
end

function check_single_airfoil(coord, label::String;
                               max_tmax, max_camber, min_aft_t)
    reasons = String[]
    t_max   = maximum(coord[:, 2]) - minimum(coord[:, 2])
    camber  = (maximum(coord[:, 2]) + minimum(coord[:, 2])) / 2.0
    aft_min = _min_aft_thickness(coord; xc_start=0.50, xc_end=0.92)

    t_max > max_tmax && push!(reasons, "$label too thick")
    abs(camber) > max_camber && push!(reasons, "$label too cambered")
    (!isinf(aft_min) && aft_min < min_aft_t) && push!(reasons, "$label thin tail")
    has_crossover(coord; xc_start=0.30) && push!(reasons, "$label crossed tail")
    return isempty(reasons), reasons
end

# =============================================================================
# LATIN HYPERCUBE SAMPLING
# =============================================================================

"""
    lhs_candidates(lb, ub, n; oversample=4) → Vector{Vector{Float64}}

Generate `n * oversample` LHS candidates in the [lb, ub] space.
"""
function lhs_candidates(lb::Vector{Float64}, ub::Vector{Float64},
                        n::Int; oversample::Int = 4)
    n_pts  = n * oversample
    n_vars = length(lb)
    pts    = Vector{Vector{Float64}}(undef, n_pts)
    for d in 1:n_vars
        perm = randperm(n_pts)
        for i in 1:n_pts
            u = (perm[i] - 1 + rand()) / n_pts
            v = lb[d] + u * (ub[d] - lb[d])
            if d == 1
                pts[i] = Vector{Float64}(undef, n_vars)
            end
            pts[i][d] = v
        end
    end
    shuffle!(pts)
    return pts
end

# =============================================================================
# CONSTRUCTIVE SAMPLER  (LHS-backed)
# =============================================================================

"""
    sample_feasible_single(lb, ub, cfg, dt, label; ..., lhs_pool=nothing)

Sample a geometrically feasible point.
"""
function sample_feasible_single(lb, ub, cfg, dt, label;
                                 n_wu::Int, n_wl::Int,
                                 max_tmax, max_camber, min_aft_t,
                                 max_attempts::Int = 500,
                                 bias=:none,
                                 lhs_pool::Union{Nothing, Vector{Vector{Float64}}} = nothing,
                                 lhs_idx::Union{Nothing, Ref{Int}} = nothing)

    function next_candidate()
        if lhs_pool !== nothing && lhs_idx !== nothing
            idx = lhs_idx[]
            if idx <= length(lhs_pool)
                lhs_idx[] += 1
                return copy(lhs_pool[idx])
            end
        end
        return [lb[i] + rand() * (ub[i] - lb[i]) for i in eachindex(lb)]
    end

    for attempt in 1:max_attempts
        x = next_candidate()
        x .= clamp.(x, lb, ub)

        if bias == :front_loaded
            wl1_idx = n_wu + 1
            wl2_idx = n_wu + 2
            wl3_idx = n_wu + 3
            r1 = lb[wl1_idx] + rand() * (ub[wl1_idx] - lb[wl1_idx]) * 0.30
            x[wl1_idx] = clamp(r1, lb[wl1_idx], ub[wl1_idx])
            mid2 = (lb[wl2_idx] + ub[wl2_idx]) / 2.0
            x[wl2_idx] = clamp(mid2 + rand() * (ub[wl2_idx] - mid2), lb[wl2_idx], ub[wl2_idx])
            if n_wl >= 3
                mid3 = (lb[wl3_idx] + ub[wl3_idx]) / 2.0
                x[wl3_idx] = clamp(mid3 + rand() * (ub[wl3_idx] - mid3), lb[wl3_idx], ub[wl3_idx])
            end
        end

        x[1] <= abs(x[n_wu+1]) && continue
        n_wl >= 2 && x[min(2, n_wu)] <= abs(x[n_wu + n_wl]) && continue

        wu_v = x[1:n_wu]
        wl_v = x[n_wu+1:end]
        try
            coord = param_to_coords(wu_v, wl_v, cfg, dt)
            pass, _ = check_single_airfoil(coord, label;
                max_tmax=max_tmax, max_camber=max_camber, min_aft_t=min_aft_t)
            if pass && !has_crossover(coord; xc_start=0.30)
                return x
            end
        catch
            continue
        end
    end
    return nothing
end

end # module
