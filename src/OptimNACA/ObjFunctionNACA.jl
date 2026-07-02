module ObjFunctionNACA

# =============================================================================
#  Standalone objective function for NACA Kriging.
#  No penalty, no dependency on SequentialEval or CST.
#  Depends on: ConfigGeom, naca_cord, XfoilBuiltin, AlphaOpt,
#              Viterna, conversions, CCBladeAnalysis
# =============================================================================

using Surrogates
using Statistics
using ..NACACord
using ..Config
using ..XfoilRunner
using ..Conversions
using ..AlphaOpt
using ..Viterna
using ..CCBladeAnalysis

using DataFrames

export obj_func_load, krig_error_minim, OBJECTIVE

const _EVAL_SPACER = "  -     -     -     -     -     -     -     -"

const OBJECTIVE = Config.OBJECTIVE_NACA


# =============================================================================
# NACA OBJECTIVE FUNCTION  (no penalty)
# =============================================================================

function obj_func_load(x; label="", objective=OBJECTIVE)
    n_points = 25

    fname         = NACACord.naca_shape_func(x, n_points)
    naca_dat_path = fname

    if !isfile(naca_dat_path)
        error("NACA airfoil file not created: $naca_dat_path")
    end

    lines    = readlines(naca_dat_path)
    coords   = [parse.(Float64, split(strip(l))) for l in lines[2:end] if !isempty(strip(l))]
    x_coords = [c[1] for c in coords]
    y_coords = [c[2] for c in coords]

    rm(naca_dat_path)

    xfoil_params = Config.create_xfoil_params()

    df_polar = XfoilRunner.run_xfoil(
        x_coords, y_coords;
        Re          = xfoil_params.Reynolds,
        Mach        = xfoil_params.Mach,
        alpha_start = xfoil_params.alpha_start,
        alpha_end   = xfoil_params.alpha_end,
        alpha_step  = xfoil_params.alpha_step,
        iter        = xfoil_params.Iterations,
        ncrit       = xfoil_params.ncrit
    )

    if df_polar === nothing
        println("$label NACA  →  XFOIL failed — returning 0.0")
        return 0.0
    end

    alpha_opt, cl_opt, cd_opt, cm_opt = AlphaOpt.find_optimal_alpha(df_polar)
    clcd = cl_opt / cd_opt

    if objective == :clcd
        println("$label NACA$(round(Int, x[1]))$(round(Int, x[2]))12  →  Cl/Cd = $(round(clcd, digits=1))  Cm=$(round(cm_opt, digits=4))")
        println(_EVAL_SPACER)
        return clcd
    end

    paths      = Config.create_paths(tempdir())
    n_sections = Config.TURBINE.n_sections
    turbine    = Config.turbine_parameters(alpha_opt, cl_opt, cd_opt, n_sections)

    df_ext, _ = Viterna.apply_viterna_extrapolation(
        df_polar, turbine.r, turbine.chord, turbine.Rtip,
        paths.xfoildata_extended_csv
    )
    Conversions.write_ccblade_airfoil(
        df_ext, paths.ccblade_af_extended,
        xfoil_params.Reynolds, xfoil_params.Mach
    )

    if objective == :cp_robust
        _, _, _, tsr_opt, cp_max, cp_robust = CCBladeAnalysis.analyze_cp_robust(
            turbine, paths.ccblade_af_extended
        )
        println("$label NACA  →  Cp_robust = $(round(cp_robust, digits=4))  Cp_max = $(round(cp_max, digits=4))  Cm=$(round(cm_opt, digits=4))  (Cl/Cd = $(round(clcd, digits=1)))")
        println(_EVAL_SPACER)
        return cp_robust
    end

    # :cp_max
    _, _, _, _, cp_max = CCBladeAnalysis.analyze_performance(
        turbine, paths.ccblade_af_extended
    )
    println("$label NACA  →  Cp_max = $(round(cp_max, digits=4))  Cm=$(round(cm_opt, digits=4))  (Cl/Cd = $(round(clcd, digits=1)))")
    println(_EVAL_SPACER)
    return cp_max
end


# =============================================================================
# KRIGING ERROR MINIMIZATION
# =============================================================================

function krig_error_minim(x_fine, krig_model::Kriging, obj_func,
                           x_samples, y_samples,
                           lower_bound::AbstractVector,
                           upper_bound::AbstractVector,
                           maxtol;
                           max_iter::Int = 100)

    iter            = 0
    outlier_indices = Int[]

    while true
        std_error              = Surrogates.std_error_at_point.(krig_model, x_fine)
        sqr_std_error          = abs2.(std_error)
        maximum_sqr_std_error, max_err_index = findmax(sqr_std_error)

        println("  └─ Max error: $(round(maximum_sqr_std_error, digits=6))")

        if maximum_sqr_std_error < maxtol
            println("Error below threshold — converged")
            break
        end

        if iter >= max_iter
            println("Reached maximum of $max_iter iterations")
            break
        end

        x_new_point = x_fine[max_err_index]

        if x_new_point in x_samples
            println("Point already sampled — stopping")
            break
        end

        y_new_point = obj_func(x_new_point; label="[Iteration n. $(iter+1)]")
        iter += 1

        push!(x_samples, x_new_point)
        push!(y_samples, y_new_point)

        y_clean = [y_samples[i] for i in eachindex(y_samples) if i ∉ outlier_indices]
        if length(y_clean) >= 4
            y_mean = mean(y_clean)
            y_std  = std(y_clean)
            if y_std > 0 && abs(y_new_point - y_mean) > 4.0 * y_std
                println("  ⚠ Outlier detected — included in Kriging but excluded from best")
                push!(outlier_indices, length(y_samples))
            end
        end

        try
            krig_model = Kriging(x_samples, y_samples, lower_bound, upper_bound,
                                 p = [2.0, 2.0])
        catch e
            println("Kriging failed: $e — stopping")
            pop!(x_samples)
            pop!(y_samples)
            break
        end
    end

    println("Minimization completed in $iter iterations")
    return krig_model, x_samples, y_samples, outlier_indices
end

end # module
