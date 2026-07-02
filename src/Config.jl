module Config

using CCBlade

using ..BladeGeom: generate_sections, calculate_chord_twist, apply_root_taper!, thickness_distribution

# =============================================================================
# INPUT FILES
# =============================================================================
include("input_turbine.jl")
include("input_optimizer.jl")

# =============================================================================
# NACA BOUNDS (aliases for backward compatibility with NACA_kriging.jl)
# =============================================================================
const nDOE = NACA_NDOE
const lb1 = NACA_LB1
const ub1 = NACA_UB1
const lb2 = NACA_LB2
const ub2 = NACA_UB2

# =============================================================================
# EXPORTS
# =============================================================================
export AirfoilParameters, TurbineParameters, XFOILParams, Paths
export airfoil_root_parameters, airfoil_tip_parameters,
       airfoil_parameters, airfoil2_parameters,
       turbine_parameters, create_xfoil_params, create_paths
export DT_ROOT, DT_TIP, DT_EXPONENT
export ROOT_TAPER_ENABLED, CHORD_ROOT_TARGET, ROOT_TAPER_ETA
export DELTA_ALPHA
export get_thickness_distribution
export compute_re_root_tip
export AIR_RHO, AIR_MU
export OBJECTIVE_CST, OBJECTIVE_NACA
export N_DOE, N_TOP, MAX_EVALS, ETA_ROOT, ETA_TIP
export FTOL_REL, XTOL_REL, COBYLA_RHOBEG
export STALL_PENALTY_K, CM_CAP_TIP, CM_CAP_ROOT
export NACA_NDOE, NACA_LB1, NACA_UB1, NACA_LB2, NACA_UB2
export NACA_THICKNESS, NACA_TOL_PERCENTAGE

# =============================================================================
# STRUCTS
# =============================================================================
struct AirfoilParameters
    wu::Vector{Float64}
    wl::Vector{Float64}
    dz::Float64
    dt::Float64
    N::Int
end

struct TurbineParameters
    Rhub::Float64
    Rtip::Float64
    B::Int
    rotorR::Float64
    precone::Float64
    yaw::Float64
    tilt::Float64
    pitch::Float64
    V_inf::Float64
    TSR::Float64
    omega::Float64
    rho::Float64
    hubHt::Float64
    shearExp::Float64
    azimuth::Vector{Float64}
    r::Vector{Float64}
    chord::Vector{Float64}
    twist::Vector{Float64}
    rotor::Rotor
end

struct XFOILParams
    Reynolds::Float64
    Mach::Float64
    Iterations::Int
    alpha_start::Float64
    alpha_end::Float64
    alpha_step::Float64
    ncrit::Float64
end

struct Paths
    dat_path::String
    run_dir::String
    polar_file::String
    xfoildata_csv::String
    xfoildata_extended_csv::String
    ccblade_af_extended::String
end

# =============================================================================
# FUNCTIONS
# =============================================================================

function airfoil_root_parameters()
    return AirfoilParameters(
        collect(AIRFOIL_ROOT.wu),
        collect(AIRFOIL_ROOT.wl),
        AIRFOIL_ROOT.dz,
        AIRFOIL_ROOT.dt,
        AIRFOIL_ROOT.N
    )
end

airfoil_parameters() = airfoil_root_parameters()

function airfoil_tip_parameters()
    return AirfoilParameters(
        collect(AIRFOIL_TIP.wu),
        collect(AIRFOIL_TIP.wl),
        AIRFOIL_TIP.dz,
        AIRFOIL_TIP.dt,
        AIRFOIL_TIP.N
    )
end

airfoil2_parameters() = airfoil_tip_parameters()

function create_xfoil_params(;
    Re::Union{Float64, Nothing}          = nothing,
    Mach::Union{Float64, Nothing}        = nothing,
    Iterations::Union{Int, Nothing}      = nothing,
    alpha_start::Union{Float64, Nothing} = nothing,
    alpha_end::Union{Float64, Nothing}   = nothing,
    alpha_step::Union{Float64, Nothing}  = nothing,
    ncrit::Union{Float64, Nothing}       = nothing)

    return XFOILParams(
        isnothing(Re)          ? XFOIL.Reynolds    : Re,
        isnothing(Mach)        ? XFOIL.Mach        : Mach,
        isnothing(Iterations)  ? XFOIL.Iterations  : Iterations,
        isnothing(alpha_start) ? XFOIL.alpha_start : alpha_start,
        isnothing(alpha_end)   ? XFOIL.alpha_end   : alpha_end,
        isnothing(alpha_step)  ? XFOIL.alpha_step  : alpha_step,
        isnothing(ncrit)       ? XFOIL.ncrit       : ncrit
    )
end

"""
    turbine_parameters(alpha_opt_tip, cl_opt_tip, cd_opt_tip, n_sections; alpha_opt_root=NaN)

Build the turbine geometry (chord, twist, rotor).
"""
function turbine_parameters(alpha_opt_tip::Real, cl_opt_tip::Real, cd_opt_tip::Real, n_sections::Int;
                            alpha_opt_root::Real = NaN)
    Rhub     = TURBINE.Rhub
    Rtip     = TURBINE.Rtip
    B        = TURBINE.B
    precone  = TURBINE.precone
    yaw      = TURBINE.yaw
    tilt     = TURBINE.tilt
    pitch    = TURBINE.pitch
    V_inf    = TURBINE.V_inf
    TSR      = TURBINE.TSR
    rho      = TURBINE.rho
    hubHt    = TURBINE.hubHt
    shearExp = TURBINE.shearExp
    azimuth  = TURBINE.azimuth

    rotorR = Rtip * cos(precone)
    omega  = TSR * V_inf / rotorR

    r = generate_sections(Rhub, Rtip, n_sections)
    chord, twist = calculate_chord_twist(alpha_opt_tip, cl_opt_tip, cd_opt_tip, r, Rtip, TSR, B, V_inf;
                                         Rhub              = Rhub,
                                         root_taper        = ROOT_TAPER_ENABLED,
                                         chord_root_target = CHORD_ROOT_TARGET,
                                         eta_taper         = ROOT_TAPER_ETA,
                                         alpha_opt_root    = alpha_opt_root,
                                         delta_alpha       = DELTA_ALPHA)
    rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)

    return TurbineParameters(
        Rhub, Rtip, B, rotorR,
        precone, yaw, tilt, pitch,
        V_inf, TSR, omega, rho,
        hubHt, shearExp, azimuth,
        r, chord, twist, rotor
    )
end

function create_paths(base_dir::String=tempdir())
    run_dir = joinpath(base_dir, "results_temp")
    mkpath(run_dir)

    return Paths(
        joinpath(run_dir, "airfoil.dat"),
        run_dir,
        joinpath(run_dir, "polar.out"),
        joinpath(run_dir, "xfoildata.csv"),
        joinpath(run_dir, "xfoildata_extended.csv"),
        joinpath(run_dir, "xfoil_airfoil_extended.dat")
    )
end

# =============================================================================
# THICKNESS DISTRIBUTION HELPER
# =============================================================================
function get_thickness_distribution(r::AbstractVector{<:Real})
    return thickness_distribution(r, TURBINE.Rhub, TURBINE.Rtip,
                                  DT_ROOT, DT_TIP; exponent=DT_EXPONENT)
end

# =============================================================================
# REYNOLDS COMPUTATION
# =============================================================================
"""
    compute_re_root_tip(V_inf, TSR, Rtip, Rhub, B, polar_func; ...) → (Re_root, Re_tip)

Iteratively compute Reynolds number at root and tip.
"""
function compute_re_root_tip(V_inf::Real, TSR::Real,
                             Rtip::Real, Rhub::Real, B::Int,
                             polar_func;
                             rho::Real     = AIR_RHO,
                             mu::Real      = AIR_MU,
                             Re_init::Real = 1e6,
                             tol::Real     = 1e-3,
                             maxiter::Int  = 30)

    Omega = TSR * V_inf / Rtip

    function _solve_re(r_eval)
        x = TSR * r_eval / Rtip

        a = 1/3
        for _ in 1:50
            x2 = x^2
            f  = 16a^3 - 24a^2 + a*(9 - 3x2) - 1 + x2
            fp = 48a^2 - 48a  + (9 - 3x2)
            abs(fp) < 1e-10 && break
            a_new = clamp(a - f/fp, 0.01, 0.95)
            abs(a_new - a) < 1e-6 && (a = a_new; break)
            a = a_new
        end

        ap  = (1 - 3a) / (4a - 1)
        phi = atan((1 - a) / ((1 + ap) * x))

        V_ax  = V_inf * (1 - a)
        V_tan = Omega * r_eval * (1 + ap)
        V_rel = sqrt(V_ax^2 + V_tan^2)

        Re = Re_init
        for _ in 1:maxiter
            Cl, Cd = polar_func(Re, phi)
            cn     = Cl * cos(phi) + Cd * sin(phi)
            chord  = (8π * a * x * sin(phi)^2 * Rtip) /
                     ((1 - a) * B * cn * TSR)
            Re_new = rho * V_rel * chord / mu
            abs(Re_new - Re) / Re < tol && return Re_new
            Re = Re_new
        end
        return Re
    end

    Re_root = _solve_re(Rhub)
    Re_tip  = _solve_re(Rtip)

    return Re_root, Re_tip
end

end
