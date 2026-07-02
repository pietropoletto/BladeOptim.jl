# =============================================================================
#   Dual_airfoil_optimizer.jl  (v3 — modular)
#
#   PHASE 1: TIP airfoil optimization  (root polar cached)
#   PHASE 2: ROOT airfoil optimization (tip polar cached)
#
#   Core logic lives in src/OptimCST/SequentialEval.jl
# =============================================================================

include("../src/BladeOptimFull.jl")
using .BladeOptimFull

const PROJECT_ROOT = dirname(@__DIR__)
using DataFrames, Printf, NLopt, Plots, Dates, XLSX, Statistics

Plots.plotlyjs()

const SE = BladeOptimFull.SequentialEval

# =============================================================================
# CONFIGURATION  (from src/input_optimizer.jl)
# =============================================================================
const CG = BladeOptimFull.Config

const N_DOE     = CG.N_DOE
const N_TOP     = CG.N_TOP
const MAX_EVALS = CG.MAX_EVALS
const ETA_ROOT  = CG.ETA_ROOT
const ETA_TIP   = CG.ETA_TIP
const FTOL_REL  = CG.FTOL_REL
const XTOL_REL  = CG.XTOL_REL
const COBYLA_RHOBEG = CG.COBYLA_RHOBEG

const TIP_BIAS_FRACTION = 0.30
const PENALTY_VALUE     = -1e6

# =============================================================================
# SETUP
# =============================================================================

airfoil_root_cfg = BladeOptimFull.Config.airfoil_root_parameters()
airfoil_tip_cfg  = BladeOptimFull.Config.airfoil_tip_parameters()
dt_root = airfoil_root_cfg.dt
dt_tip  = airfoil_tip_cfg.dt

objective_str = string(BladeOptimFull.ObjFunction.OBJECTIVE)
obj_label = Dict("cp_max"=>"Cp_max","clcd"=>"Cl/Cd","cp_robust"=>"Cp_robust")[objective_str]

n_wu = length(airfoil_root_cfg.wu)
n_wl = length(airfoil_root_cfg.wl)
n_single = n_wu + n_wl

Rhub = BladeOptimFull.Config.TURBINE.Rhub
Rtip = BladeOptimFull.Config.TURBINE.Rtip
n_sections = BladeOptimFull.Config.TURBINE.n_sections
delta_alpha_geom = BladeOptimFull.Config.DELTA_ALPHA

# =============================================================================
# REYNOLDS
# =============================================================================
let
    _V   = BladeOptimFull.Config.TURBINE.V_inf
    _TSR = BladeOptimFull.Config.TURBINE.TSR
    _B   = BladeOptimFull.Config.TURBINE.B
    _rho = BladeOptimFull.Config.TURBINE.rho

    _polar(Re, phi) = (0.9, 0.01)

    global Re_root, Re_tip = BladeOptimFull.Config.compute_re_root_tip(
        _V, _TSR, Rtip, Rhub, _B, _polar; rho=_rho
    )
end

# =============================================================================
# OUTPUT DIRECTORY  — Results_cst/<timestamp>/
# =============================================================================

_ts      = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
run_dir  = joinpath(PROJECT_ROOT, "Results_cst", _ts)
mkpath(run_dir)
_t_start = time()


const _sim_log = Vector{NamedTuple}()

function is_sane(val, obj)
    (isnan(val) || isinf(val)) && return false
    obj in ("cp_max","cp_robust") && return -0.5 <= val <= 0.58
    obj == "clcd" && return -50 <= val <= 200
    return abs(val) < 1e4
end

# =============================================================================
# BOUNDS
# =============================================================================

function build_bounds(wu_c, wl_c, delta)
    lb = Float64[]; ub = Float64[]
    for i in eachindex(wu_c)
        push!(lb, max(wu_c[i]-delta, 0.05)); push!(ub, min(wu_c[i]+delta, 0.8))
    end
    for i in eachindex(wl_c)
        push!(lb, max(wl_c[i]-delta, -0.8))
        push!(ub, i==1 ? min(wl_c[i]+delta, -0.02) : min(wl_c[i]+delta, 0.8))
    end
    lb[1] = max(lb[1], 0.10)
    ub[length(wu_c)+1] = min(ub[length(wu_c)+1], -0.10)
    return lb, ub
end

delta = 0.5
wu_center_root = collect(airfoil_root_cfg.wu)
wl_center_root = collect(airfoil_root_cfg.wl)
wu_center_tip  = collect(airfoil_tip_cfg.wu)
wl_center_tip  = collect(airfoil_tip_cfg.wl)

lb_tip,  ub_tip  = build_bounds(wu_center_tip,  wl_center_tip,  delta)
lb_root, ub_root = build_bounds(wu_center_root, wl_center_root, delta)

constraint_le(x::Vector, ::Vector) = abs(x[n_wu+1]) - x[1]
constraint_te(x::Vector, ::Vector) = abs(x[n_wu+n_wl]) - x[2]

# =============================================================================
# HEADER
# =============================================================================

println()
println("=" ^ 72)
println("  DUAL-AIRFOIL SEQUENTIAL OPTIMIZER")
println("-" ^ 72)
println("  Objective    $obj_label")
println("  Turbine      R=$(Rtip)m  B=$(BladeOptimFull.Config.TURBINE.B)  TSR=$(BladeOptimFull.Config.TURBINE.TSR)")
println("  Reynolds     Re_root=$(round(Int, Re_root))  (η=$(ETA_ROOT))   Re_tip=$(round(Int, Re_tip))  (η=$(ETA_TIP))")
println("  Root         $(n_wu)wu + $(n_wl)wl   t/c=$(round(dt_root*100,digits=1))%")
println("  Tip          $(n_wu)wu + $(n_wl)wl   t/c=$(round(dt_tip*100,digits=1))%")
println("  Twist Δα     $(delta_alpha_geom)°  (geometry offset on α_root)")
println("  CST bounds   ±$(delta)")
println("  Cm cap       tip=$(BladeOptimFull.ObjFunction.CM_CAP_TIP)  root=$(BladeOptimFull.ObjFunction.CM_CAP_ROOT)")
println("  DOE          $N_DOE LHS points  →  top $N_TOP  →  COBYLA $MAX_EVALS eval")
println("  Budget/fase  ~$(N_DOE + N_TOP * MAX_EVALS) total evals")
println("  Convergence  ftol=$(FTOL_REL)  xtol=$(XTOL_REL)  rhobeg=$(COBYLA_RHOBEG)")
println("  Results →    $run_dir")
println("=" ^ 72)

# =============================================================================
# SAVE HELPERS  (top-level so Julia 1.12 world-age rules are satisfied)
# =============================================================================

function _interp1_col(df_src, col::Symbol, a_grid)
    vals = zeros(length(a_grid))
    col_data = df_src[!, col]
    for (i, a) in enumerate(a_grid)
        if a <= df_src.alpha[1]
            vals[i] = col_data[1]
        elseif a >= df_src.alpha[end]
            vals[i] = col_data[end]
        else
            idx = findlast(df_src.alpha .<= a)
            t   = (a - df_src.alpha[idx]) / (df_src.alpha[idx+1] - df_src.alpha[idx])
            vals[i] = col_data[idx] + t * (col_data[idx+1] - col_data[idx])
        end
    end
    return vals
end

function _write_df!(sheet, df, row0, col0; header=true)
    col_names = String.(names(df))
    if header
        for (c, nm) in enumerate(col_names)
            sheet[row0, col0 + c - 1] = nm
        end
        row0 += 1
    end
    for (r, row) in enumerate(eachrow(df))
        for (c, val) in enumerate(row)
            v = ismissing(val) ? "" :
                (val isa AbstractFloat && isnan(val)) ? "" : val
            sheet[row0 + r - 1, col0 + c - 1] = v
        end
    end
end

# =============================================================================
#  GENERIC DOE + NLopt PHASE
# =============================================================================


function _write_restart_block!(sheet, phase_df, restart_id, col0, penalty_col_label)
    block_df = filter(r -> r.restart == restart_id, phase_df)
    sheet[2, col0] = "REStart $restart_id"
    headers = ["N SIM", "CP MAX", "CP PENALTY", penalty_col_label, "CM"]
    for (c, h) in enumerate(headers)
        sheet[3, col0 + c - 1] = h
    end
    for (r, row) in enumerate(eachrow(block_df))
        _nan(v) = v isa AbstractFloat && isnan(v) ? "" : v
        sheet[3 + r, col0]     = r
        sheet[3 + r, col0 + 1] = _nan(row.cp_max)
        sheet[3 + r, col0 + 2] = _nan(row.cp_pen)
        sheet[3 + r, col0 + 3] = _nan(row.penalty)
        cm_val = (row.phase == "TIP") ? row.cm_tip : row.cm_root
        sheet[3 + r, col0 + 4] = _nan(cm_val)
    end
end

function run_phase(label_phase, lb, ub, cfg, dt, df_fixed;
                   alpha_forced=NaN, max_tmax, max_camber, min_aft_t,
                   bias_fraction=TIP_BIAS_FRACTION,
                   is_tip_variable::Bool,
                   Re_xfoil::Float64,
                   Re_root_eval::Float64,
                   Re_tip_eval::Float64)

    baseline_x = vcat(collect(cfg.wu), collect(cfg.wl))

    # ── LHS pool: generate N_DOE×5 candidates with uniform coverage ────────────
    # oversample=5 ensures enough candidates even with ~30% geometric rejection.
    _lhs_pool = SE.lhs_candidates(lb, ub, N_DOE; oversample=5)
    _lhs_idx  = Ref(1)

    # ── DOE: baseline + N_DOE-1 LHS feasible points ───────────────────────────
    samples = Vector{Vector{Float64}}()
    push!(samples, baseline_x)   # baseline sempre come primo punto

    n_gen = 1
    for i in 2:N_DOE
        use_bias = rand() < bias_fraction
        sx = SE.sample_feasible_single(lb, ub, cfg, dt, label_phase;
            n_wu=n_wu, n_wl=n_wl,
            max_tmax=max_tmax, max_camber=max_camber, min_aft_t=min_aft_t,
            bias     = use_bias ? :front_loaded : :none,
            lhs_pool = _lhs_pool,
            lhs_idx  = _lhs_idx)
        sx === nothing && continue
        push!(samples, sx)
        n_gen += 1
    end
    println("  DOE (LHS): $n_gen / $N_DOE  (pool: $(_lhs_idx[]-1) / $(length(_lhs_pool)) candidates used)")

    # ── Evaluate all DOE points ──────────────────────────────────────────────
    values = Float64[]
    for (i, x) in enumerate(samples)
        wu_v = x[1:n_wu]; wl_v = x[n_wu+1:end]
        df_new = SE.run_xfoil_single(wu_v, wl_v, cfg, dt; Re=Re_xfoil)
        if df_new === nothing
            push!(values, PENALTY_VALUE); continue
        end
        df_r, df_t = is_tip_variable ? (df_fixed, df_new) : (df_new, df_fixed)
        val, _, _, _, _, _, _ = SE.evaluate_from_polars(df_r, df_t;
            alpha_root_forced=alpha_forced,
            label="[$label_phase DOE $i/$n_gen]",
            Re_root=Re_root_eval, Re_tip=Re_tip_eval,
            variable_phase=is_tip_variable ? :tip : :root)
        push!(values, is_sane(val, objective_str) ? val : PENALTY_VALUE)
    end

    sorted = sortperm(values, rev=true)
    n_show = min(N_TOP, n_gen)
    println("  DOE top-$n_show:")
    for r in 1:n_show
        println("    #$r  $obj_label = $(round(values[sorted[r]], digits=5))")
    end

    # ── COBYLA refinement on top-N_TOP ────────────────────────────────────
    best_val   = -Inf
    best_x     = nothing
    eval_ct    = Ref(0)
    restart_ct = Ref(0)
    nlopt_results = Vector{Tuple{Float64, Vector{Float64}, Int, String}}()

    println()
    for rank in 1:min(N_TOP, n_gen)
        i         = sorted[rank]
        x0        = copy(samples[i])
        start_val = values[i]
        eval_ct[] = 0
        restart_ct[] += 1
        _restart_id = restart_ct[]

        # Initial step proportional to mean bound width
        # → first steps explore ~COBYLA_RHOBEG × range of the space
        _rhobeg = COBYLA_RHOBEG * mean(ub .- lb)

        opt = NLopt.Opt(:LN_COBYLA, n_single)
        NLopt.lower_bounds!(opt, lb); NLopt.upper_bounds!(opt, ub)
        NLopt.initial_step!(opt, _rhobeg)
        NLopt.max_objective!(opt, (x, grad) -> begin
            eval_ct[] += 1
            wu_v = x[1:n_wu]; wl_v = x[n_wu+1:end]
            df_new = SE.run_xfoil_single(wu_v, wl_v, cfg, dt; Re=Re_xfoil)
            df_new === nothing && return PENALTY_VALUE
            df_r, df_t = is_tip_variable ? (df_fixed, df_new) : (df_new, df_fixed)
            v, _, _, _, _, cm_t, cm_r = SE.evaluate_from_polars(df_r, df_t;
                alpha_root_forced=alpha_forced,
                label="[$label_phase #$(eval_ct[])]",
                Re_root=Re_root_eval, Re_tip=Re_tip_eval,
                variable_phase=is_tip_variable ? :tip : :root)
            sane    = is_sane(v, objective_str)
            final_v = sane ? v : PENALTY_VALUE
            push!(_sim_log, (
                phase   = label_phase,
                restart = _restart_id,
                cp_max  = sane ? v : NaN,
                cp_pen  = final_v,
                penalty = sane ? 0.0 : abs(PENALTY_VALUE),
                cm_root = cm_r,
                cm_tip  = cm_t
            ))
            final_v
        end)
        NLopt.maxeval!(opt, MAX_EVALS)
        NLopt.ftol_rel!(opt, FTOL_REL)
        NLopt.xtol_rel!(opt, XTOL_REL)
        NLopt.inequality_constraint!(opt, constraint_le, 1e-8)
        NLopt.inequality_constraint!(opt, constraint_te, 1e-8)

        bv, bx, ret = NLopt.optimize(opt, x0)

        δ = is_sane(bv, objective_str) && is_sane(start_val, objective_str) ?
            (bv - start_val) / abs(start_val) * 100 : NaN
        δ_str = isnan(δ) ? "" : @sprintf("  Δ=%+.2f%%", δ)
        tag   = is_sane(bv, objective_str) && bv > best_val ? "  ★" : ""
        println("  $label_phase $rank/$n_show  $obj_label=$(round(bv,digits=5))  $(eval_ct[]) evals  $ret$δ_str$tag")

        if is_sane(bv, objective_str)
            push!(nlopt_results, (bv, copy(bx), eval_ct[], string(ret)))
            if bv > best_val
                best_val = bv; best_x = copy(bx)
            end
        end
    end

    sort!(nlopt_results, by=first, rev=true)
    total_evals = n_gen + min(N_TOP, n_gen) * MAX_EVALS
    return best_x, best_val, values[1], nlopt_results, total_evals
end

# =============================================================================
#  PHASE 1 — TIP
# =============================================================================

println()
println("=" ^ 72)
println("  PHASE 1 — TIP  (root polar cached)")
println("=" ^ 72)

df_root_cached = SE.run_xfoil_single(wu_center_root, wl_center_root, airfoil_root_cfg, dt_root;
                                     Re=Re_root)
df_root_cached === nothing && error("ROOT baseline XFOIL failed!")
_α_rc, _cl_rc, _cd_rc, _cm_rc = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_root_cached)
println("  Root polar cached  →  α=$(round(_α_rc,digits=1))°  Cl/Cd=$(round(_cl_rc/_cd_rc,digits=1))  Cm=$(round(_cm_rc,digits=4))  Re=$(round(Int,Re_root))")

best_tip_x, best_tip_val, tip_baseline, tip_nlopt_results, tip_doe_n = run_phase(
    "TIP", lb_tip, ub_tip, airfoil_tip_cfg, dt_tip, df_root_cached;
    max_tmax=0.20, max_camber=0.08, min_aft_t=0.01,
    bias_fraction=TIP_BIAS_FRACTION, is_tip_variable=true,
    Re_xfoil=Re_tip, Re_root_eval=Re_root, Re_tip_eval=Re_tip)

best_wu_tip = best_tip_x[1:n_wu]
best_wl_tip = best_tip_x[n_wu+1:end]

df_tip_best = SE.run_xfoil_single(best_wu_tip, best_wl_tip, airfoil_tip_cfg, dt_tip; Re=Re_tip)
alpha_tip_final, cl_tip_f, cd_tip_f, cm_tip_f = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_tip_best)

println("-" ^ 72)
println("  TIP result   $obj_label=$(round(best_tip_val,digits=5))  α=$(round(alpha_tip_final,digits=1))°  Cl/Cd=$(round(cl_tip_f/cd_tip_f,digits=1))  Cm=$(round(cm_tip_f,digits=4))")

# =============================================================================
#  PHASE 2 — ROOT
# =============================================================================

println()
println("=" ^ 72)
println("  PHASE 2 — ROOT  (tip polar cached)")
println("=" ^ 72)

df_tip_cached = df_tip_best
_α_tc, _cl_tc, _cd_tc, _cm_tc = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_tip_cached)
println("  Tip polar cached   →  α=$(round(_α_tc,digits=1))°  Cl/Cd=$(round(_cl_tc/_cd_tc,digits=1))  Cm=$(round(_cm_tc,digits=4))  Re=$(round(Int,Re_tip))")

best_root_x, best_root_val, root_baseline, root_nlopt_results, root_doe_n = run_phase(
    "ROOT", lb_root, ub_root, airfoil_root_cfg, dt_root, df_tip_cached;
    max_tmax=0.30, max_camber=0.070, min_aft_t=0.04,
    is_tip_variable=false,
    Re_xfoil=Re_root, Re_root_eval=Re_root, Re_tip_eval=Re_tip)

best_wu_root = best_root_x[1:n_wu]
best_wl_root = best_root_x[n_wu+1:end]

df_root_best = SE.run_xfoil_single(best_wu_root, best_wl_root, airfoil_root_cfg, dt_root; Re=Re_root)
alpha_root_final, cl_root_f, cd_root_f, cm_root_f = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_root_best)

println("-" ^ 72)
println("  ROOT result  $obj_label=$(round(best_root_val,digits=5))  α=$(round(alpha_root_final,digits=1))°  Cl/Cd=$(round(cl_root_f/cd_root_f,digits=1))  Cm=$(round(cm_root_f,digits=4))")

# =============================================================================
#  FINAL ANALYSIS
# =============================================================================

println()
println("=" ^ 72)
println("  FINAL BLADE ANALYSIS")
println("=" ^ 72)

df_root_final = df_root_best
df_tip_final  = df_tip_cached

if df_root_final === nothing
    println("  ROOT XFOIL failed — skipping.")
else
    a_tip, cl_tip, cd_tip, cm_tip = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_tip_final)
    a_root_nat, cl_root, cd_root, cm_root = BladeOptimFull.AlphaOpt.find_optimal_alpha(df_root_final)

    turbine = BladeOptimFull.Config.turbine_parameters(a_tip, cl_tip, cd_tip, n_sections;
                                                              alpha_opt_root = a_root_nat)

    all_α = let
        raw = sort(unique(vcat(df_root_final.alpha, df_tip_final.alpha)))
        lo = max(minimum(df_root_final.alpha), minimum(df_tip_final.alpha))
        hi = min(maximum(df_root_final.alpha), maximum(df_tip_final.alpha))
        filter(a -> lo <= a <= hi, raw)
    end

    cl_rg, cd_rg = BladeOptimFull.ObjFunction._interp_polar(df_root_final, all_α)
    cl_tg, cd_tg = BladeOptimFull.ObjFunction._interp_polar(df_tip_final,  all_α)

    xp = BladeOptimFull.Config.create_xfoil_params()

    af_files   = String[]
    sec_frames = DataFrame[]

    for j in 1:n_sections
        η = clamp((turbine.r[j]-Rhub)/(Rtip-Rhub), 0.0, 1.0)
        w = SE.blend_weight(η)
        cl_j = (1-w).*cl_rg .+ w.*cl_tg
        cd_j = (1-w).*cd_rg .+ w.*cd_tg
        df_j = DataFrame(alpha=all_α, cl=cl_j, cd=cd_j)
        df_ext, _ = BladeOptimFull.Viterna.apply_viterna_extrapolation(
            df_j, turbine.r, turbine.chord, turbine.Rtip)
        re_j = Re_root + clamp(η, 0.0, 1.0) * (Re_tip - Re_root)
        af = tempname() * ".dat"
        BladeOptimFull.Conversions.write_ccblade_airfoil(df_ext, af, re_j, xp.Mach)
        push!(af_files, af)
        push!(sec_frames, df_ext)
    end

    tsrvec, cpvec, ctvec, tsr_opt, cp_max =
        BladeOptimFull.CCBladeAnalysis.analyze_performance_multipolar(turbine, af_files)
    T, Q, P, Np, Tp =
        BladeOptimFull.CCBladeAnalysis.analyze_loads_multipolar(turbine, af_files)
    ct_opt = ctvec[argmin(abs.(tsrvec .- tsr_opt))]

    Δcp = (cp_max - tip_baseline) / abs(tip_baseline) * 100

    println()
    println("  Performance")
    println("    Cp_max     $(round(cp_max,digits=5))")
    println("    CT         $(round(ct_opt,digits=5))")
    println("    TSR_opt    $(round(tsr_opt,digits=2))")
    println("    Power      $(round(P/1000,digits=2)) kW")
    println("    ΔCp        $(Δcp>=0 ? "+" : "")$(round(Δcp,digits=2))% vs baseline")
    println()
    println("  Aerodynamics")
    println("    α_tip      $(round(a_tip,digits=1))°   Cl/Cd=$(round(cl_tip/cd_tip,digits=1))   Cm=$(round(cm_tip,digits=4))")
    println("    α_root     $(round(a_root_nat,digits=1))°   Cl/Cd=$(round(cl_root/cd_root,digits=1))   Cm=$(round(cm_root,digits=4))")
    da = delta_alpha_geom
    if da != 0.0
        println("    α_root+Δ   $(round(a_root_nat+da,digits=1))°   (twist design, Δα=$(da)°)")
    end

    println()
    println("-" ^ 72)
    println("  COPY-PASTE INTO input_turbine.jl:")
    println("-" ^ 72)
    println("""
const AIRFOIL_ROOT = (
    wu = [$(join(round.(best_wu_root,digits=6), ", "))],
    wl = [$(join(round.(best_wl_root,digits=6), ", "))],
    dz = $(airfoil_root_cfg.dz), dt = $(dt_root), N = $(airfoil_root_cfg.N))

const AIRFOIL_TIP = (
    wu = [$(join(round.(best_wu_tip,digits=6), ", "))],
    wl = [$(join(round.(best_wl_tip,digits=6), ", "))],
    dz = $(airfoil_tip_cfg.dz), dt = $(dt_tip), N = $(airfoil_tip_cfg.N))
""")

    # =========================================================================
    #  SAVE — results.xlsx  +  cst_coefficients.txt
    # =========================================================================

    _has_cm_root = hasproperty(df_root_final, :cm)
    _has_cm_tip  = hasproperty(df_tip_final,  :cm)

    sec_polar_frames = DataFrame[]
    for j in 1:n_sections
        η      = clamp((turbine.r[j]-Rhub)/(Rtip-Rhub), 0.0, 1.0)
        w      = SE.blend_weight(η)
        df_sec = copy(sec_frames[j])
        if _has_cm_root && _has_cm_tip
            cm_r_g = _interp1_col(df_root_final, :cm, df_sec.alpha)
            cm_t_g = _interp1_col(df_tip_final,  :cm, df_sec.alpha)
            df_sec[!, :cm] = (1-w) .* cm_r_g .+ w .* cm_t_g
        else
            df_sec[!, :cm] = fill(NaN, nrow(df_sec))
        end
        push!(sec_polar_frames, df_sec)
    end

    _omega_sv  = tsr_opt * turbine.V_inf / turbine.Rtip
    _lam_sv    = collect(tsrvec)
    _cp_sv     = collect(cpvec)
    _U_sv      = _omega_sv .* turbine.Rtip ./ _lam_sv
    _Pavail_sv = 0.5 .* turbine.rho .* _U_sv.^3 .* (π * turbine.Rtip^2)
    _P_sv      = _cp_sv .* _Pavail_sv

    XLSX.openxlsx(joinpath(run_dir, "results.xlsx"), mode="w") do xf

        sh1 = xf[1]; XLSX.rename!(sh1, "Cp_Ct")
        _write_df!(sh1, DataFrame(
            lambda = collect(tsrvec),
            Cp     = collect(cpvec),
            Ct     = collect(ctvec)), 1, 1)

        XLSX.addsheet!(xf, "Power_Curve")
        _write_df!(xf["Power_Curve"], DataFrame(
            lambda    = _lam_sv,
            U_inf_ms  = _U_sv,
            P_avail_W = _Pavail_sv,
            P_W       = _P_sv,
            P_kW      = _P_sv ./ 1e3), 1, 1)

        XLSX.addsheet!(xf, "Chord_Twist")
        _write_df!(xf["Chord_Twist"], DataFrame(
            r         = turbine.r,
            chord     = turbine.chord,
            twist_deg = rad2deg.(turbine.twist)), 1, 1)

        XLSX.addsheet!(xf, "Polars")
        sh4 = xf["Polars"]
        for j in 1:n_sections
            col0 = 1 + (j-1) * 5
            sh4[1, col0] = "Section $j  (r = $(round(turbine.r[j], digits=3)) m)"
            _write_df!(sh4, sec_polar_frames[j], 2, col0)
        end

        # ── Sheet 5: SimLog ─────────────────────────────────────────────────
        # Layout (mirrors the screenshot):
        #
        #   col:  A          ...        |  (gap)  |  L         ...
        #   row1:     TIP                              ROOT
        #   row2: REStart 1  ... | gap | REStart 2     REStart 1 ... | gap | REStart 2
        #   row3: N SIM | CP MAX | CP PENALTY | PENALTY TIP | CM   (repeated per block)
        #   row4+: data
        #
        # Each restart block is 5 cols wide + 1 col gap between restarts.
        # Gap of 1 col between TIP and ROOT macro-blocks.

        XLSX.addsheet!(xf, "SimLog")
        sh5 = xf["SimLog"]

        if !isempty(_sim_log)
            df_log = DataFrame(_sim_log)

            # Separate TIP and ROOT entries
            tip_data  = filter(r -> r.phase == "TIP",  df_log)
            root_data = filter(r -> r.phase == "ROOT", df_log)

            tip_restarts  = sort(unique(tip_data.restart))
            root_restarts = sort(unique(root_data.restart))

            _COLS_PER_RESTART = 5   # N SIM | CP MAX | CP PENALTY | PENALTY | CM
            _GAP = 1                # empty col between restart blocks
            _BLOCK_GAP = 2          # empty cols between TIP and ROOT macro-blocks

            # Number of restart blocks on left (TIP) side
            n_tip_blocks = length(tip_restarts)

            # Starting column for ROOT macro-block
            root_start_col = 1 + n_tip_blocks * (_COLS_PER_RESTART + _GAP) + _BLOCK_GAP - _GAP

            # ── TIP macro-block ──
            sh5[1, 1] = "TIP"
            for (k, rid) in enumerate(tip_restarts)
                col0 = 1 + (k-1) * (_COLS_PER_RESTART + _GAP)
                _write_restart_block!(sh5, tip_data, rid, col0, "PENALTY TIP")
            end

            # ── ROOT macro-block ──
            sh5[1, root_start_col] = "ROOT"
            for (k, rid) in enumerate(root_restarts)
                col0 = root_start_col + (k-1) * (_COLS_PER_RESTART + _GAP)
                _write_restart_block!(sh5, root_data, rid, col0, "PENALTY ROOT")
            end
        else
            sh5[1, 1] = "No simulation data logged."
        end
    end

    open(joinpath(run_dir, "cst_coefficients.txt"), "w") do io
        for (s_idx, (lbl, wu, wl, dt_v, cfg)) in enumerate([
                ("ROOT", best_wu_root, best_wl_root, dt_root, airfoil_root_cfg),
                ("TIP",  best_wu_tip,  best_wl_tip,  dt_tip,  airfoil_tip_cfg)])
            println(io, "Section_$(s_idx) ($lbl):")
            println(io, "    wu = [$(join(round.(wu, digits=7), ", "))],")
            println(io, "    wl = [$(join(round.(wl, digits=7), ", "))],")
            println(io, "    dz = $(cfg.dz),")
            println(io, "    dt = $(round(dt_v, digits=6)),")
            println(io, "    N  = $(cfg.N)")
            println(io)
        end
    end

    coord_base    = SE.param_to_coords(wu_center_root, wl_center_root, airfoil_root_cfg, dt_root)
    coord_root_sc = SE.param_to_coords(best_wu_root, best_wl_root, airfoil_root_cfg, dt_root)
    coord_tip_sc  = SE.param_to_coords(best_wu_tip,  best_wl_tip,  airfoil_tip_cfg,  dt_tip)
    r_root_eval   = Rhub + ETA_ROOT*(Rtip-Rhub)
    r_tip_eval    = Rhub + ETA_TIP*(Rtip-Rhub)

    BladeOptimFull.OptimizerPlots.plot_optimizer_results(
        coord_base, coord_root_sc, coord_tip_sc,
        dt_root, dt_tip,
        df_root_final, df_tip_final,
        a_root_nat, a_tip,
        cl_root, cd_root, cl_tip, cd_tip,
        tsrvec, cpvec, ctvec, tsr_opt, cp_max, ct_opt,
        turbine, Np, Tp,
        r_root_eval, r_tip_eval,
        tip_baseline, run_dir;
        delta_alpha = delta_alpha_geom)

    # ── Report.txt ────────────────────────────────────────────────────────────
    _t_elapsed = time() - _t_start
    _t_h = floor(Int, _t_elapsed / 3600)
    _t_m = floor(Int, (_t_elapsed % 3600) / 60)
    _t_s = floor(Int, _t_elapsed % 60)
    _fmt_time = @sprintf("%02dh %02dm %02ds", _t_h, _t_m, _t_s)

    open(joinpath(run_dir, "report.txt"), "w") do io
        _line  = "=" ^ 80
        _dash  = "-" ^ 80
        _fv(x, d=5) = string(round(x, digits=d))

        println(io, _line)
        println(io, "  OPTIMIZATION REPORT  —  $(_ts)")
        println(io, _line)
        println(io)
        println(io, "  TURBINE CONFIG")
        println(io, "    Rotor radius      $(Rtip) m              Blades           $(BladeOptimFull.Config.TURBINE.B)")
        println(io, "    Design TSR        $(BladeOptimFull.Config.TURBINE.TSR)                Wind speed       $(BladeOptimFull.Config.TURBINE.V_inf) m/s")
        println(io, "    Air density       $(BladeOptimFull.Config.TURBINE.rho) kg/m³          Reynolds root    $(round(Int,Re_root))  tip $(round(Int,Re_tip))")
        println(io, "    Hub radius        $(Rhub) m              Sections         $(n_sections)")
        println(io)
        println(io, _dash)
        println(io, "  PERFORMANCE")
        println(io, "    Cp max            $(_fv(cp_max))         $(Δcp>=0 ? "+" : "")$(_fv(Δcp,2))% vs baseline")
        println(io, "    Ct @ TSR_opt      $(_fv(ct_opt))")
        println(io, "    TSR_opt           $(_fv(tsr_opt,2))")
        println(io, "    Power             $(round(P/1000,digits=2)) kW")
        println(io)
        println(io, _dash)
        println(io, "  AERODYNAMICS")
        println(io, "                        ROOT              TIP")
        println(io, "    α opt             $(lpad(round(a_root_nat,digits=1),8))°         $(lpad(round(a_tip,digits=1),8))°")
        println(io, "    Cl                $(lpad(round(cl_root,digits=4),8))          $(lpad(round(cl_tip,digits=4),8))")
        println(io, "    Cd                $(lpad(round(cd_root,digits=5),8))         $(lpad(round(cd_tip,digits=5),8))")
        println(io, "    Cl/Cd             $(lpad(round(cl_root/cd_root,digits=1),8))          $(lpad(round(cl_tip/cd_tip,digits=1),8))")
        println(io, "    Cm                $(lpad(round(cm_root,digits=4),8))         $(lpad(round(cm_tip,digits=4),8))")
        println(io)
        println(io, _dash)
        println(io, "  CST COEFFICIENTS")
        println(io, "    ROOT  wu = [$(join(round.(best_wu_root,digits=6), ", "))]")
        println(io, "          wl = [$(join(round.(best_wl_root,digits=6), ", "))]")
        println(io, "          dt = $(dt_root)   dz = $(airfoil_root_cfg.dz)   N = $(airfoil_root_cfg.N)")
        println(io)
        println(io, "    TIP   wu = [$(join(round.(best_wu_tip,digits=6), ", "))]")
        println(io, "          wl = [$(join(round.(best_wl_tip,digits=6), ", "))]")
        println(io, "          dt = $(dt_tip)   dz = $(airfoil_tip_cfg.dz)   N = $(airfoil_tip_cfg.N)")
        println(io)
        println(io, _dash)
        println(io, "  OPTIMIZER STATS")
        println(io, "    Total sim time    $(_fmt_time)")
        println(io)
        println(io, "  Strategy     DOE (LHS) + COBYLA refinement")
        println(io, "  Config       DOE=$N_DOE (LHS)  top=$N_TOP  evals/restart=$MAX_EVALS")
        println(io)
        println(io, "    Phase TIP    total evals=$(tip_doe_n)")
        isempty(tip_nlopt_results) || begin
            _best_tip_cp = maximum(r[1] for r in tip_nlopt_results)
            for (k, res) in enumerate(tip_nlopt_results)
                cp_r, _, evals_r, ret_r = res
                best_tag = (cp_r == _best_tip_cp) ? "  ★" : ""
                println(io, "      Restart $k   $(lpad(evals_r,3)) evals   result  $(_fv(cp_r))  $(ret_r)$best_tag")
            end
        end
        println(io)
        println(io, "    Phase ROOT   total evals=$(root_doe_n)")
        isempty(root_nlopt_results) || begin
            _best_root_cp = maximum(r[1] for r in root_nlopt_results)
            for (k, res) in enumerate(root_nlopt_results)
                cp_r, _, evals_r, ret_r = res
                best_tag = (cp_r == _best_root_cp) ? "  ★" : ""
                println(io, "      Restart $k   $(lpad(evals_r,3)) evals   result  $(_fv(cp_r))  $(ret_r)$best_tag")
            end
        end
        println(io)
        println(io, _line)
    end

    println("  Results saved → $run_dir")
end

println()
println("=" ^ 72)
println("  DONE")
println("=" ^ 72)

# =============================================================================
#  3D BLADE VISUALIZATION
# =============================================================================

if @isdefined(turbine)
    p_blade, blade_sections = BladeOptimFull.BladeVisual.plot_blade_3d(
        best_wu_root, best_wl_root, airfoil_root_cfg,
        best_wu_tip,  best_wl_tip,  airfoil_tip_cfg,
        turbine;
        n_pts     = 80,
        title_str = "Optimized Blade — Cp=$(round(cp_max, digits=4))",
        save_path = joinpath(run_dir, "blade_3d.html")
    )
else
    println("  Skipping 3D visualization (turbine not defined)")
end






