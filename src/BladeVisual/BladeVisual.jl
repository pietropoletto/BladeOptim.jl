module BladeVisual

# =============================================================================
#   blade_visual.jl  —  GLMakie 3D blade visualization
#
#   LOGICA SPAN:
#     η ∈ [0.00, 0.20]  →  blend cerchio → root airfoil
#     η ∈ [0.20, 0.40]  →  pure root airfoil
#     η ∈ [0.40, 0.90]  →  blend root → tip
#     η ∈ [0.90, 1.00]  →  pure tip airfoil
#
#   SISTEMA COORDINATE (pala):
#     X = span (radiale, hub→tip)
#     Y = chord (LE→TE nel piano del disco)
#     Z = spessore (out-of-plane)
#
#   STACKING: quarter chord (25%) sull'asse X (span).
#   TWIST:    rotazione in piano YZ attorno al quarter chord.
#
#   UTILIZZO:
#     fig, sections = BladeVisual.plot_blade_3d(
#         best_wu_root, best_wl_root, airfoil_root_cfg,
#         best_wu_tip,  best_wl_tip,  airfoil_tip_cfg,
#         turbine;
#         n_pts         = 100,
#         n_sec_shown   = 10,
#         color_hub     = RGB(0.88, 0.39, 0.32),
#         color_tip     = RGB(0.47, 0.55, 0.66),
#         save_path     = joinpath(plots_dir, "blade_3d.png")
#     )
#
#   DIPENDENZE: GLMakie (finestra interattiva desktop), Colors
# =============================================================================

export plot_blade_3d, build_blade_sections

using GLMakie
using LinearAlgebra
using Colors

const ETA_CYL_END   = 0.20
const ETA_ROOT_PURE = 0.40
const ETA_TIP_START = 0.90

@inline _smoothstep(t::Real) = t * t * (3 - 2t)

# =============================================================================
#  CST — internal replica of CST.jl
# =============================================================================
function _class_shape(w, x, N1, N2, dz)
    n = length(w) - 1
    K = [binomial(n, j-1) for j in 1:(n+1)]
    y = zeros(length(x))
    for i in eachindex(x)
        xi = x[i]
        C  = xi^N1 * (1 - xi)^N2
        S  = sum(w[j] * K[j] * xi^(j-1) * (1-xi)^(n-(j-1)) for j in 1:(n+1))
        y[i] = C*S + xi*dz
    end
    return y
end

function _airfoil_norm(wu, wl, dz, dt, n_pts)
    zeta  = range(0, 2π, length = n_pts + 1)
    x_all = @. 0.5*(cos(zeta) + 1)
    le_i  = argmin(x_all)
    xu = x_all[1:le_i-1]
    xl = x_all[le_i:end]
    yu = _class_shape(wu, xu, 0.5, 1.0,  dz)
    yl = _class_shape(wl, xl, 0.5, 1.0, -dz)
    xc = vcat(xu, xl)
    yc = vcat(yu, yl)
    lei      = argmin(abs.(xc))
    L        = min(lei, length(xc)-lei+1)
    t_actual = maximum(yc[1:L] .- yc[end-L+1:end])
    t_actual > 1e-8 && (yc .*= dt/t_actual)
    return hcat(xc, yc)
end

function _circle_physical(r_phys, r_span, n_pts)
    theta = range(0, 2π, length = n_pts + 1)
    pts   = zeros(n_pts + 1, 3)
    for i in eachindex(theta)
        pts[i, 1] = r_span
        pts[i, 2] = r_phys * cos(theta[i])
        pts[i, 3] = r_phys * sin(theta[i])
    end
    return pts
end

# =============================================================================
#  PLACE SECTION — coordinate normalizzate → 3D fisico
# =============================================================================
function _place_section(cn, chord, r_span, twist_rad)
    n    = size(cn, 1)
    pts  = zeros(n, 3)
    cosT = cos(twist_rad)
    sinT = sin(twist_rad)
    for i in 1:n
        yp = (cn[i,1] - 0.25) * chord
        zp =  cn[i,2]         * chord
        pts[i,1] = r_span
        pts[i,2] = yp*cosT - zp*sinT
        pts[i,3] = yp*sinT + zp*cosT
    end
    return pts
end

# =============================================================================
#  BUILD BLADE SECTIONS
# =============================================================================
"""
    build_blade_sections(wu_root, wl_root, cfg_root,
                         wu_tip,  wl_tip,  cfg_tip,
                         turbine; n_pts=100)

Returns Vector{Matrix{Float64}}: each section is (n_pts+1)×3,
colonne = [X=span, Y=chord, Z=thickness].
"""
function build_blade_sections(
        wu_root, wl_root, cfg_root,
        wu_tip,  wl_tip,  cfg_tip,
        turbine; n_pts=100)

    Rhub    = turbine.Rhub
    Rtip    = turbine.Rtip
    r_v     = turbine.r
    chord_v = turbine.chord
    twist_v = turbine.twist

    af_root = _airfoil_norm(wu_root, wl_root, cfg_root.dz, cfg_root.dt, n_pts)
    af_tip  = _airfoil_norm(wu_tip,  wl_tip,  cfg_tip.dz,  cfg_tip.dt,  n_pts)

    r_hub_phys        = chord_v[1] / 2.0
    r_transition_phys = chord_v[argmin(abs.(r_v .- (Rhub + ETA_CYL_END*(Rtip-Rhub))))] * cfg_root.dt / 2.0

    sections = Vector{Matrix{Float64}}()
    for i in eachindex(r_v)
        eta = clamp((r_v[i] - Rhub)/(Rtip - Rhub), 0.0, 1.0)

        if eta <= ETA_CYL_END
            w       = _smoothstep(eta / ETA_CYL_END)
            r_cyl   = r_hub_phys + w * (r_transition_phys - r_hub_phys)
            sec_cyl = _circle_physical(r_cyl, r_v[i], n_pts)
            sec_af  = _place_section(af_root, chord_v[i], r_v[i], twist_v[i])
            sec = (1-w) .* sec_cyl .+ w .* sec_af
            push!(sections, sec)

        elseif eta <= ETA_ROOT_PURE
            push!(sections, _place_section(af_root, chord_v[i], r_v[i], twist_v[i]))

        elseif eta <= ETA_TIP_START
            w  = _smoothstep((eta - ETA_ROOT_PURE)/(ETA_TIP_START - ETA_ROOT_PURE))
            cn = (1-w) .* af_root .+ w .* af_tip
            push!(sections, _place_section(cn, chord_v[i], r_v[i], twist_v[i]))

        else
            push!(sections, _place_section(af_tip, chord_v[i], r_v[i], twist_v[i]))
        end
    end
    return sections
end

# =============================================================================
#  TRIANGOLAZIONE
# =============================================================================
function _triangulate(sections)
    n_sec     = length(sections)
    n_pts_sec = size(sections[1], 1)

    verts = zeros(Float32, n_sec * n_pts_sec, 3)
    for i in 1:n_sec
        for j in 1:n_pts_sec
            idx = (i-1)*n_pts_sec + j
            verts[idx, :] = sections[i][j, :]
        end
    end

    faces = Vector{GLMakie.GeometryBasics.TriangleFace{Int}}()
    for i in 1:(n_sec-1)
        for j in 1:(n_pts_sec-1)
            v00 = (i-1)*n_pts_sec + j
            v10 = (i  )*n_pts_sec + j
            v01 = (i-1)*n_pts_sec + j+1
            v11 = (i  )*n_pts_sec + j+1
            push!(faces, GLMakie.GeometryBasics.TriangleFace(v00, v10, v01))
            push!(faces, GLMakie.GeometryBasics.TriangleFace(v10, v11, v01))
        end
    end

    return verts, faces
end

# =============================================================================
#  _make_colormap — generate a linear two-color colormap
# =============================================================================
function _make_colormap(c1, c2; n_steps::Int=256)
    col1 = parse(Colorant, c1)
    col2 = parse(Colorant, c2)
    return cgrad([col1, col2], n_steps)
end

# =============================================================================
#  PLOT BLADE 3D  —  GLMakie
# =============================================================================
function plot_blade_3d(
        wu_root, wl_root, cfg_root,
        wu_tip,  wl_tip,  cfg_tip,
        turbine;
        n_pts::Int        = 100,
        n_sec_shown::Int  = 10,
        color_hub         = RGB(0.88, 0.39, 0.32),   # corallo
        color_tip         = RGB(0.47, 0.55, 0.66),   # grigio blu
        title_str::String = "Optimized Blade",
        save_path::Union{String,Nothing} = nothing)

    println("  [BladeVisual] Building sections...")
    secs      = build_blade_sections(wu_root, wl_root, cfg_root,
                                     wu_tip,  wl_tip,  cfg_tip,
                                     turbine; n_pts)
    n_sec     = length(secs)
    n_pts_sec = size(secs[1], 1)
    println("  [BladeVisual] $n_sec sections × $n_pts_sec pts")

    span_min = turbine.Rhub
    span_max = turbine.Rtip

    blade_cmap = _make_colormap(color_hub, color_tip)

    verts, faces = _triangulate(secs)
    span_vals = verts[:, 1]
    colors    = (span_vals .- span_min) ./ (span_max - span_min)

    bg_color   = :white
    text_color = :black
    sec_color  = (:gray40, 0.5)

    # ── figure ───────────────────────────────────────────────────────────────
    fig = Figure(size = (1600, 900), backgroundcolor = bg_color)
    ax  = Axis3(fig[1,1];
        title              = title_str,
        aspect             = :data,
        titlesize          = 20,
        titlecolor         = text_color,
        backgroundcolor    = bg_color,
        xlabelvisible      = false,
        ylabelvisible      = false,
        zlabelvisible      = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        zticklabelsvisible = false,
        xticksvisible      = false,
        yticksvisible      = false,
        zticksvisible      = false,
        xgridvisible       = false,
        ygridvisible       = false,
        zgridvisible       = false,
        xspinesvisible     = false,
        yspinesvisible     = false,
        zspinesvisible     = false,
        xypanelvisible     = false,
        xzpanelvisible     = false,
        yzpanelvisible     = false,
        azimuth            = 1.3π,
        elevation          = 0.15π,
        perspectiveness    = 0.3)

    # ── surface mesh ──────────────────────────────────────────────────────
    pts_makie = [GLMakie.Point3f(verts[i,1], verts[i,2], verts[i,3])
                 for i in 1:size(verts,1)]

    mesh!(ax, pts_makie, faces;
        color        = colors,
        colormap     = blade_cmap,
        shading      = NoShading,
        transparency = false)

    # ── leading edge / trailing edge ─────────────────────────────────────────
    le_pts = [GLMakie.Point3f(secs[i][argmin(secs[i][:,2]), :]...)
              for i in 1:n_sec]
    te_pts = [GLMakie.Point3f(secs[i][argmax(secs[i][:,2]), :]...)
              for i in 1:n_sec]

    lines!(ax, le_pts; color = RGB(0.36, 0.73, 0.47),  linewidth = 3.0)
    lines!(ax, te_pts; color = RGB(0.58, 0.46, 0.78),  linewidth = 3.0)

    # ── distributed cross-sections ──────────────────────────────────────
    idxs = unique(clamp.(round.(Int, range(2, n_sec-1, length=n_sec_shown)), 1, n_sec))
    for idx in idxs
        sec = secs[idx]
        pts = [GLMakie.Point3f(sec[j,1], sec[j,2], sec[j,3]) for j in 1:n_pts_sec]
        lines!(ax, pts; color = sec_color, linewidth = 1.0)
    end

    # ── sezioni riferimento η=40% e η=90% ───────────────────────────────────
    for (eta_mark, col, _lbl) in [
            (ETA_ROOT_PURE, :black, "Root η=40%"),
            (ETA_TIP_START, :black, "Tip η=90%")]
        r_m = turbine.Rhub + eta_mark*(turbine.Rtip - turbine.Rhub)
        idx = argmin(abs.(turbine.r .- r_m))
        sec = secs[idx]
        pts = [GLMakie.Point3f(sec[j,1], sec[j,2], sec[j,3]) for j in 1:n_pts_sec]
        lines!(ax, pts; color = col, linewidth = 2.5)
    end

    # ── legend ──────────────────────────────────────────────────────────────
    Legend(fig[1,2],
        [LineElement(color=RGB(0.36, 0.73, 0.47), linewidth=3),
         LineElement(color=RGB(0.58, 0.46, 0.78), linewidth=3),
         LineElement(color=sec_color,             linewidth=1),
         LineElement(color=:black,                linewidth=3),
         LineElement(color=:black,                linewidth=3)],
        ["Leading Edge", "Trailing Edge", "Sections",
         "Root airfoil η=40%", "Tip airfoil η=90%"],
        backgroundcolor = (:white, 0.9),
        labelcolor      = text_color,
        framecolor      = :gray70,
        padding         = (8, 8, 6, 6))

    # ── colorbar span ────────────────────────────────────────────────────────
    Colorbar(fig[1,3];
        colormap       = blade_cmap,
        limits         = (span_min, span_max),
        label          = "Span [m]",
        labelcolor     = text_color,
        ticklabelcolor = text_color)

    colsize!(fig.layout, 1, Relative(0.78))
    colsize!(fig.layout, 2, Auto())
    colsize!(fig.layout, 3, Auto())

    println("  [BladeVisual] Plot ready.")

    if save_path !== nothing
        save(save_path, fig)
        println("  [BladeVisual] Saved → $save_path")
    end

    display(fig)
    return fig, secs
end

end # module BladeVisual