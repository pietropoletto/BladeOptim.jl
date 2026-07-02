module BladeGeom

export generate_sections, calculate_chord_twist
export apply_root_taper!, thickness_distribution

# =============================================================================
#  COSINE-CLUSTERED SECTIONS
# =============================================================================
function generate_sections(Rhub::Real, Rtip::Real, n_sections::Int)
    beta = range(0, π, length=n_sections)
    return @. Rhub + (Rtip - Rhub) * (1 - cos(beta)) / 2
end

# =============================================================================
#  NEWTON-RAPHSON FOR 'a' (Betz axial induction factor)
# =============================================================================
@inline function _solve_a(x::Real; a_init=1/3, max_iter=50, tol=1e-6)
    a  = a_init
    x2 = x^2
    for _ in 1:max_iter
        f  = 16a^3 - 24a^2 + a*(9 - 3x2) - 1 + x2
        fp = 48a^2 - 48a  + (9 - 3x2)
        abs(fp) < 1e-10 && break
        a_new = clamp(a - f/fp, 0.01, 0.95)
        abs(a_new - a) < tol && (a = a_new; break)
        a = a_new
    end
    return a
end

# =============================================================================
#  CHORD AND TWIST COMPUTATION — linear interpolation with optional delta
#
#  SCHEMA:
#    φ(r)        = angolo di flusso Betz
#    η(r)        = (r - Rhub) / (Rtip - Rhub)   ∈ [0, 1]
#    α_root_des  = α_root + delta_alpha
#    α_design(r) = α_root_des + η × (α_tip - α_root_des)
#    θ(r)        = φ(r) - α_design(r)
#
#    delta_alpha > 0 → increases root α → less twist → earlier root stall
#    delta_alpha < 0 → decreases root α → more twist
#    delta_alpha = 0 → pure Betz twist with natural interpolation
#
#  If alpha_opt_root = NaN → α_design = α_tip uniform (single-airfoil)
# =============================================================================
"""
    calculate_chord_twist(alpha_opt_tip, cl_opt_tip, cd_opt_tip,
                          r, Rtip, lambda, B, V_inf; kwargs...)

Compute Betz chord and twist for each radial station.
α_design linearly interpolated from `alpha_opt_root + delta_alpha` (hub)
a `alpha_opt_tip` (tip).

Keyword arguments:
- `alpha_opt_root` : α ottimale root [deg] (NaN = uniforme α_tip)
- `delta_alpha`    : Δ [deg] — offset added to α_root in twist design (default 0.0)
- `Rhub`           : raggio hub [m]
"""
function calculate_chord_twist(
        alpha_opt_tip::Real,
        cl_opt_tip::Real,
        cd_opt_tip::Real,
        r::AbstractVector{<:Real},
        Rtip::Real,
        lambda::Real,
        B::Int,
        V_inf::Real;
        a_init::Real            = 1/3,
        max_iter::Int           = 50,
        tol::Real               = 1e-6,
        chord_min::Real         = 0.0,
        # ── Root taper ──────────────────────────────────────────────────────
        Rhub::Real              = NaN,
        root_taper::Bool        = true,
        chord_root_target::Real = 0.30,
        eta_taper::Real         = 0.20,
        # ── Twist ───────────────────────────────────────────────────────────
        alpha_opt_root::Real    = NaN,
        delta_alpha::Real       = 0.0)

    n     = length(r)
    chord = zeros(n)
    twist = zeros(n)

    use_dual = !isnan(alpha_opt_root) && !isnan(Rhub)

    alpha_tip_rad  = deg2rad(alpha_opt_tip)
    alpha_root_rad = use_dual ? deg2rad(alpha_opt_root + delta_alpha) : alpha_tip_rad

    for i in 1:n
        x  = lambda * r[i] / Rtip
        a  = _solve_a(x; a_init, max_iter, tol)
        ap = (1 - 3a)/(4a - 1)
        phi = atan((1 - a)/((1 + ap)*x))

        # ── α_design(η): linear interpolation root → tip ──────────────────
        if use_dual
            eta = clamp((r[i] - Rhub)/(Rtip - Rhub), 0.0, 1.0)
            alpha_design_rad = alpha_root_rad + eta * (alpha_tip_rad - alpha_root_rad)
        else
            alpha_design_rad = alpha_tip_rad
        end

        twist[i] = phi - alpha_design_rad

        # ── Prandtl tip loss + chord Betz ────────────────────────────────────
        f_loss   = (B*(Rtip - r[i]))/(2*r[i]*sin(phi))
        F        = (2/π)*acos(exp(-f_loss))
        cn       = cl_opt_tip*cos(phi) + cd_opt_tip*sin(phi)
        chord[i] = (8π*F*a*x*sin(phi)^2*Rtip)/((1-a)*B*cn*lambda)
        chord[i] = max(chord[i], chord_min)
    end

    # ── Root taper ───────────────────────────────────────────────────────────
    if root_taper && !isnan(Rhub)
        apply_root_taper!(chord, r, Rhub, Rtip;
                          chord_root_target = chord_root_target,
                          eta_taper         = eta_taper)
    end

    return chord, twist
end

# =============================================================================
#  ROOT TAPER  (smoothstep C1)
# =============================================================================
function apply_root_taper!(chord::AbstractVector{<:Real},
                           r::AbstractVector{<:Real},
                           Rhub::Real, Rtip::Real;
                           chord_root_target::Real = 0.30,
                           eta_taper::Real         = 0.20)
    for i in eachindex(r)
        eta = (r[i] - Rhub)/(Rtip - Rhub)
        if eta < eta_taper
            w        = eta/eta_taper
            t        = w*w*(3 - 2w)
            chord[i] = chord_root_target + t*(chord[i] - chord_root_target)
        end
    end
    return chord
end

# =============================================================================
#  THICKNESS DISTRIBUTION  (legge potenza)
# =============================================================================
function thickness_distribution(r::AbstractVector{<:Real},
                                Rhub::Real, Rtip::Real,
                                dt_root::Real, dt_tip::Real;
                                exponent::Real = 1.5)
    return [dt_tip + (dt_root - dt_tip)*(1 - clamp((r[i]-Rhub)/(Rtip-Rhub), 0, 1))^exponent
            for i in eachindex(r)]
end

end # module BladeGeom