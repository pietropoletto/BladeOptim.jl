module OptimizerPlots

using Plots
using Printf

export plot_optimizer_results, plot_top_result

"""
    plot_optimizer_results(...; delta_alpha=0.0)

`delta_alpha` è l'offset geometrico [deg] aggiunto a α_root nel twist.
Deve corrispondere a DELTA_ALPHA in Config.
"""
function plot_optimizer_results(
    coord_base_sc, coord_root_sc, coord_tip_sc,
    dt_root_eval, dt_tip_eval,
    df_root_final, df_tip_final,
    alpha_opt_r, alpha_opt_f,
    cl_opt_r, cd_opt_r, cl_opt_f, cd_opt_f,
    tsrvec_f, cpvec_f, ctvec_f, tsr_opt_f, cp_max_f, ct_at_opt_f,
    turbine_f, Np_f, Tp_f,
    r_root_eval, r_tip_eval,
    baseline_val,
    plots_dir;
    delta_alpha::Real = 0.0
)

# P1 — Airfoil shapes (root + tip overlay)
    pb1 = plot(coord_base_sc[:, 1], coord_base_sc[:, 2],
        label="",
        color=:gray, linestyle=:dash, lw=1.5,
        aspect_ratio=:equal, xlabel="x/c", ylabel="y/c",
        title="Airfoil Shapes — Root & Tip")
    plot!(pb1, coord_root_sc[:, 1], coord_root_sc[:, 2],
        label="",
         color=:red, lw=2.5)
    plot!(pb1, coord_tip_sc[:, 1], coord_tip_sc[:, 2],
        label="",
         color=:blue, lw=2.5)

    # P2 — Cl/Cd vs alpha (root and tip)
    clcd_root_f = df_root_final.cl ./ df_root_final.cd
    clcd_tip_f  = df_tip_final.cl  ./ df_tip_final.cd
    pb2 = plot(df_root_final.alpha, clcd_root_f,
        label="Root", lw=2, color=:red, legend=:topright,
        xlabel="α [deg]", ylabel="Cl/Cd", title="Aerodynamic Efficiency")
    plot!(pb2, df_tip_final.alpha, clcd_tip_f,
        label="Tip", lw=2, color=:blue)
    scatter!(pb2, [alpha_opt_r], [cl_opt_r/cd_opt_r],
        markersize=4, color=:red,  markershape=:circle, label="")
    scatter!(pb2, [alpha_opt_f], [cl_opt_f/cd_opt_f],
        markersize=4, color=:blue, markershape=:circle, label="") 

    # P3 — Cp and CT vs TSR
    pb3 = plot(tsrvec_f, cpvec_f,
        label="Cp", lw=2.5, color=:blue,
        xlabel="Tip Speed Ratio (λ)", ylabel="Cp, CT",
        title="Power and Thrust Coefficients", legend=:topleft)
    plot!(pb3, tsrvec_f, ctvec_f, label="CT", lw=2.5, color=:green)
    scatter!(pb3, [tsr_opt_f], [cp_max_f],
        label="Cp_max (λ=$(round(tsr_opt_f,digits=2)))",
        markersize=6, color=:blue, markershape=:circle)
    scatter!(pb3, [tsr_opt_f], [ct_at_opt_f],
        label="CT@λ_opt", markersize=6, color=:orange, markershape=:circle)

    # P4 — Chord distribution (STANDALONE)
    pb4 = plot(turbine_f.r, turbine_f.chord,
        label="", 
        lw=2.5, color=RGB(0.88, 0.39, 0.32), legend=false,
        xlabel="Radius [m]", ylabel="Chord [m]",
        grid=false,
        guidefontsize=18, tickfontsize=13,
        size=(1000, 600), margin=10Plots.mm)
    vline!(pb4, [r_root_eval, r_tip_eval], linestyle=:dot, color=:black,
        lw=1.5, label="")

    # P5 — Twist distribution (STANDALONE)
    _p5_Rhub   = turbine_f.Rhub
    _p5_Rtip   = turbine_f.Rtip
    _p5_r      = turbine_f.r
    _p5_lam    = turbine_f.TSR
    _p5_n      = length(_p5_r)
    _p5_phi    = zeros(_p5_n)
    for _p5_i in 1:_p5_n
        _p5_x = _p5_lam * _p5_r[_p5_i] / _p5_Rtip
        _p5_a = 1/3
        for _ in 1:50
            _p5_fp  = 16_p5_a^3 - 24_p5_a^2 + _p5_a*(9 - 3_p5_x^2) - 1 + _p5_x^2
            _p5_fpp = 48_p5_a^2 - 48_p5_a + (9 - 3_p5_x^2)
            abs(_p5_fpp) < 1e-10 && break
            _p5_an = clamp(_p5_a - _p5_fp/_p5_fpp, 0.01, 0.95)
            abs(_p5_an - _p5_a) < 1e-6 && (_p5_a = _p5_an; break)
            _p5_a = _p5_an
        end
        _p5_ap       = (1 - 3_p5_a)/(4_p5_a - 1)
        _p5_phi[_p5_i] = atan((1 - _p5_a)/((1 + _p5_ap) * _p5_x))
    end
    _p5_phi_deg      = rad2deg.(_p5_phi)
    _p5_twist_betz   = _p5_phi_deg .- alpha_opt_f

    _p5_alpha_rd   = alpha_opt_r + delta_alpha
    _p5_eta        = clamp.((_p5_r .- _p5_Rhub) ./ (_p5_Rtip - _p5_Rhub), 0.0, 1.0)
    _p5_alpha_des  = _p5_alpha_rd .+ _p5_eta .* (alpha_opt_f - _p5_alpha_rd)
    _p5_twist_dual = _p5_phi_deg .- _p5_alpha_des

    pb5 = plot(_p5_r, _p5_twist_dual;
        label = "",
        lw = 2.5, color = RGB(0.27, 0.60, 0.81),
        xlabel = "Radius [m]", ylabel = "Twist [deg]",
        legend = false, grid = false,
        guidefontsize = 18, tickfontsize = 13,
        size = (1000, 600), margin = 10Plots.mm)
    hline!(pb5, [0]; linestyle=:dot, color=:gray, lw=1, label="")
    vline!(pb5, [r_root_eval, r_tip_eval]; linestyle=:dot, color=:black, lw=1.5, label="")

    # P6 — Blade planform
    le_f = zeros(length(turbine_f.r))
    te_f = turbine_f.chord
    pb6  = plot(turbine_f.r, le_f, lw=2, color=:black, linestyle=:solid, label="")
    plot!(pb6, turbine_f.r, te_f, lw=2, color=:black, linestyle=:solid, label="")
    plot!(pb6, turbine_f.r, le_f, fillrange=te_f, fillalpha=0.3, fillcolor=:lightblue, label="")
    plot!(pb6, [turbine_f.r[1],   turbine_f.r[1]],   [le_f[1],   te_f[1]],   lw=2, color=:black, label="")
    plot!(pb6, [turbine_f.r[end], turbine_f.r[end]], [le_f[end], te_f[end]], lw=2, color=:black, label="")
    xlabel!(pb6, "Radius [m]"); ylabel!(pb6, "Chord [m]"); title!(pb6, "Blade Planform (Top View)")
    plot!(pb6, aspect_ratio=:equal)

    # P7 — Cl curves
    pb7 = plot(df_root_final.alpha, df_root_final.cl,
        label="",
         lw=2, color=:red, legend=:topleft,
        xlabel="α [deg]", ylabel="Cl", title="Lift Coefficient")
    plot!(pb7, df_tip_final.alpha, df_tip_final.cl,
        label="",
         lw=2, color=:blue)
    scatter!(pb7, [alpha_opt_r], [cl_opt_r], markersize=4, color=:red,  markershape=:circle, label="")
    scatter!(pb7, [alpha_opt_f], [cl_opt_f], markersize=4, color=:blue, markershape=:circle, label="")

    # P8 — Cd curves
    pb8 = plot(df_root_final.alpha, df_root_final.cd,
        label="",
         lw=2, color=:red, legend=:topleft,
        xlabel="α [deg]", ylabel="Cd", title="Drag Coefficient")
    plot!(pb8, df_tip_final.alpha, df_tip_final.cd,
        label="",
        lw=2, color=:blue)
    scatter!(pb8, [alpha_opt_r], [cd_opt_r], markersize=4, color=:red,  markershape=:circle, label="")
    scatter!(pb8, [alpha_opt_f], [cd_opt_f], markersize=4, color=:blue, markershape=:circle, label="")

    # P9 — Blade loads
    pb9 = plot(turbine_f.r, Np_f,
        label="Flapwise (Np)", lw=2.5, color=:blue, legend=:topleft,
        xlabel="Radius [m]", ylabel="Load [N/m]", title="Blade Loads")
    plot!(pb9, turbine_f.r, Tp_f, label="Lead-lag (Tp)", lw=2.5, color=:orange)

    # P10 — Cm curves
    _has_cm_root = hasproperty(df_root_final, :cm)
    _has_cm_tip  = hasproperty(df_tip_final,  :cm)
    pb10 = plot(xlabel="α [deg]", ylabel="Cm", title="Pitching Moment Coefficient", legend=:topright)
    if _has_cm_root
        plot!(pb10, df_root_final.alpha, df_root_final.cm,
           label="",
            lw=2, color=:red)
        _cm_r_opt = df_root_final.cm[argmin(abs.(df_root_final.alpha .- alpha_opt_r))]
        scatter!(pb10, [alpha_opt_r], [_cm_r_opt], markersize=4, color=:red, markershape=:circle, label="")
    end
    if _has_cm_tip
        plot!(pb10, df_tip_final.alpha, df_tip_final.cm,
            label="", 
            lw=2, color=:blue)
        _cm_f_opt = df_tip_final.cm[argmin(abs.(df_tip_final.alpha .- alpha_opt_f))]
        scatter!(pb10, [alpha_opt_f], [_cm_f_opt], markersize=4, color=:blue, markershape=:circle, label="")
    end
    hline!(pb10, [0], linestyle=:dot, color=:gray, lw=1, label="")

    # --- Compose figures ---
    # Chord and Twist are now standalone — compose remaining panels
    fig_blade = plot(pb1, pb6,
        layout=(1, 2), size=(1200, 500), margin=8Plots.mm,
        plot_title="Dual-Airfoil Blade — Geometry",
        plot_titlefontsize=16)

    fig_aero = plot(pb7, pb8, pb10, pb9, pb2,
        layout=(3, 2), size=(1200, 1200), margin=8Plots.mm,
        plot_title="Dual-Airfoil Blade — Aerodynamic Polars & Loads",
        plot_titlefontsize=16)

    fig_cp = plot(pb3,
        layout=(1, 1), size=(1200, 800), margin=8Plots.mm,
        plot_title="Dual-Airfoil Blade — Power & Thrust",
        plot_titlefontsize=16)

    display(fig_blade)
    display(fig_aero)
    display(fig_cp)
    display(pb4)
    display(pb5)

    savefig(fig_blade, joinpath(plots_dir, "dual_blade_geometry.png"))
    savefig(fig_aero,  joinpath(plots_dir, "dual_blade_aero.png"))
    savefig(fig_cp,    joinpath(plots_dir, "dual_blade_cp.png"))
    savefig(pb4,       joinpath(plots_dir, "dual_blade_chord.png"))
    savefig(pb5,       joinpath(plots_dir, "dual_blade_twist.png"))

end

"""
    plot_top_result(rank, phase, obj_val, coord, df_root, df_tip,
                    tsrvec, cpvec, ctvec, tsr_opt, cp_max, ct_opt, plots_dir)

Generate a 2×2 page for result #rank of phase `phase`:
  top-left  : Cp & Ct vs TSR
  top-right : Airfoil shape
  bottom-left : Cm vs α  (root + tip)
  bottom-right: Cl vs α  (root + tip)
"""
function plot_top_result(
    rank::Int, phase::String, obj_val::Float64,
    coord, df_root, df_tip,
    tsrvec, cpvec, ctvec, tsr_opt, cp_max, ct_opt,
    plots_dir::String)

    tag = "$phase #$rank — Cp=$(round(cp_max, digits=4))"

    # ── top-right: Airfoil shape ──
    p_af = plot(coord[:, 1], coord[:, 2],
        label="$phase #$rank", lw=2.5, color=:black,
        aspect_ratio=:equal, xlabel="x/c", ylabel="y/c",
        title="Airfoil — $tag")

    # ── top-left: Cp & Ct vs TSR ──
    p_cpct = plot(tsrvec, cpvec,
        label="Cp", lw=2.5, color=:blue,
        xlabel="TSR (λ)", ylabel="Cp, Ct",
        title="Cp & Ct — $tag", legend=:topleft)
    plot!(p_cpct, tsrvec, ctvec, label="Ct", lw=2.5, color=:green)
    scatter!(p_cpct, [tsr_opt], [cp_max],
        label="Cp_max=$(round(cp_max,digits=4))  λ=$(round(tsr_opt,digits=2))",
        markersize=6, color=:blue, markershape=:circle)
    scatter!(p_cpct, [tsr_opt], [ct_opt],
        label="Ct@λ_opt", markersize=6, color=:orange, markershape=:circle)

    # ── bottom-right: Cl vs α ──
    p_cl = plot(xlabel="α [deg]", ylabel="Cl", title="Cl — $tag", legend=:topleft)
    plot!(p_cl, df_root.alpha, df_root.cl, label="ROOT", lw=2, color=:red)
    plot!(p_cl, df_tip.alpha,  df_tip.cl,  label="TIP",  lw=2, color=:blue)

    # ── bottom-left: Cm vs α ──
    p_cm = plot(xlabel="α [deg]", ylabel="Cm", title="Cm — $tag", legend=:topright)
    if hasproperty(df_root, :cm)
        plot!(p_cm, df_root.alpha, df_root.cm, label="ROOT", lw=2, color=:red)
    end
    if hasproperty(df_tip, :cm)
        plot!(p_cm, df_tip.alpha, df_tip.cm, label="TIP", lw=2, color=:blue)
    end
    hline!(p_cm, [0], linestyle=:dot, color=:gray, lw=1, label="")

    # ── Compose: layout (2,2) ──
    fig = plot(p_cpct, p_af, p_cm, p_cl,
        layout=(2, 2), size=(1200, 800), margin=8Plots.mm,
        plot_title="Top Result — $tag",
        plot_titlefontsize=14)

    display(fig)

    fname = "top5_$(lowercase(phase))_$(rank).png"
    savefig(fig, joinpath(plots_dir, fname))
    println("    #$rank  $tag  → $fname")
end

end