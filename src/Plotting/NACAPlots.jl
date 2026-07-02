module NACAPlots

using Plots
using Plots: RGB

export plot_naca_airfoil

function plot_naca_airfoil(
    naca_surrogate,
    naca_surrogate_initial,
    x1_range, x2_range,
    x1s_doe, x2s_doe,
    sample_x1, sample_x2,
    xvals, yvals,
    best_x, best_value,
    obj_label,
    output_dir,
    x_true_range,
    y_true
)

    # --- Custom color palettes ---
    true_palette   = cgrad([:white, RGB(0.47, 0.55, 0.66)])   # grigio blu
    before_palette = cgrad([:white, RGB(0.36, 0.73, 0.47)])   # verde salvia
    after_palette  = cgrad([:white, RGB(0.88, 0.39, 0.32)])   # corallo

    # --- True NACA map (3D) ---

    if y_true !== nothing
        p1 = Plots.surface(x1_range, x2_range, y_true,
                title="",
                xlabel="Camber", ylabel="Camber position", zlabel=obj_label,
                size=(900, 700), c=true_palette)
        display(p1)
    end

    # --- True NACA map (2D contour) ---

    if y_true !== nothing
        p1c = Plots.contourf(x1_range, x2_range, y_true,
                title="",
                xlabel="Camber", ylabel="Camber position",
                colorbar_title=obj_label,
                size=(900, 700), c=true_palette, levels=20, linewidth=0.5,
                guidefontsize=14, titlefontsize=16, tickfontsize=11,
                colorbar_titlefontsize=13,
                legendfontsize = 12,
                legend=(0.555, 0.95),
                xlims=(x1_range[1], x1_range[end]), ylims=(x2_range[1], x2_range[end]))
        # Find and mark the maximum on the true map
        max_idx = argmax(y_true)
        max_x1 = x1_range[max_idx[2]]
        max_x2 = x2_range[max_idx[1]]
        Plots.scatter!(p1c, [max_x1], [max_x2],
                markersize=5, markercolor=:black, markershape=:star, label="Max $obj_label (true)")
        display(p1c)
    end

    # --- Kriging BEFORE error minimization (3D) ---

    y_kriging_before = [naca_surrogate_initial([x1, x2]) for x2 in x2_range, x1 in x1_range]

    p2 = Plots.surface(x1_range, x2_range, y_kriging_before,
            title="",
            xlabel="Camber", ylabel="Camber position", zlabel=obj_label,
            size=(900, 700), c=before_palette)
    y_samples_initial = [naca_surrogate_initial([x1s_doe[i], x2s_doe[i]]) for i in 1:length(x1s_doe)]
    Plots.scatter!(p2, x1s_doe, x2s_doe, y_samples_initial,
            markersize=1.5, markercolor=:black, label="DOE samples")
    display(p2)

    # --- Kriging BEFORE error minimization (2D contour) ---

    p2c = Plots.contourf(x1_range, x2_range, y_kriging_before,
            title="",
            xlabel="Camber", ylabel="Camber position",
            colorbar_title=obj_label,
            size=(900, 700), c=before_palette, levels=20, linewidth=0.5,
            guidefontsize=14, titlefontsize=16, tickfontsize=11, legendfontsize = 12,
            xlims=(x1_range[1], x1_range[end]), ylims=(x2_range[1], x2_range[end]),
            legend=(0.575, 0.95))
    Plots.scatter!(p2c, x1s_doe, x2s_doe,
            markersize=2, markercolor=:black, label="DOE samples")
    display(p2c)

    # --- Kriging AFTER error minimization (3D) ---

    y_kriging_after = [naca_surrogate([x1, x2]) for x2 in yvals, x1 in xvals]

    p3 = Plots.surface(xvals, yvals, y_kriging_after,
            title="",
            xlabel="Camber", ylabel="Camber position", zlabel=obj_label,
            size=(900, 700), c=after_palette)
    y_all_samples = [naca_surrogate([sample_x1[i], sample_x2[i]]) for i in 1:length(sample_x1)]
    Plots.scatter!(p3, sample_x1, sample_x2, y_all_samples,
            markersize=1.5, markercolor=:black, label="Error samples")
    Plots.scatter!(p3, [best_x[1]], [best_x[2]], [best_value],
            markersize=4, markercolor=:red, markershape=:diamond, label="Best airfoil (kriging)")
    display(p3)

    # --- Kriging AFTER error minimization (2D contour) ---

    p3c = Plots.contourf(xvals, yvals, y_kriging_after,
            title="",
            xlabel="Camber", ylabel="Camber position",
            colorbar_title=obj_label,
            size=(900, 700), c=after_palette, levels=20, linewidth=0.5,
            guidefontsize=14, titlefontsize=16, tickfontsize=11, legendfontsize = 12,
            xlims=(xvals[1], xvals[end]), ylims=(yvals[1], yvals[end]),
            legend=(0.53, 0.95))
    Plots.scatter!(p3c, sample_x1, sample_x2,
            markersize=2, markercolor=:black, label="Error samples")
    Plots.scatter!(p3c, [best_x[1]], [best_x[2]],
            markersize=5, markercolor=:black, markershape=:star, label="Best airfoil (kriging)")
    display(p3c)

    # --- True NACA map vs Kriging overlay (3D) ---

    if y_true !== nothing
        true_x1 = [x[1] for x in x_true_range]
        true_x2 = [x[2] for x in x_true_range]
        p4 = Plots.surface(x1_range, x2_range, y_true,
                title="",
                xlabel="Camber", ylabel="Camber position", zlabel=obj_label,
                size=(900, 700), c=true_palette, alpha=1.0)
        Plots.surface!(p4, xvals, yvals, y_kriging_after, c=after_palette, alpha=0.8)
        Plots.scatter!(p4, true_x1, true_x2, vec(y_true),
                markersize=1.5, markercolor=:black, label="NACA true samples")
        display(p4)
    end

end

end