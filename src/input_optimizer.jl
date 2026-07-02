# =============================================================================
#   INPUT — OPTIMIZER SETTINGS
#
#   All optimization parameters, split into two sections:
#     1. NACA Kriging optimizer
#     2. CST dual-airfoil optimizer
#
#   Read by Config_geom.jl, obj_function.jl, obj_function_naca.jl,
#   and Optimizer.jl.
# =============================================================================

# =============================================================================
# XFOIL PARAMETERS (shared by NACA and CST pipelines)
# =============================================================================
const XFOIL = (
    Reynolds    = 620000,                         # Reynolds number for polar computation
    Mach        = 0.0,                             # freestream Mach number
    Iterations  = 500,                             # max XFOIL viscous iterations per α
    alpha_start = -4.0,                            # sweep start angle [deg]
    alpha_end   = 24,                              # sweep end angle [deg]
    alpha_step  = 0.5,                             # angle increment [deg]
    ncrit       = 6.0                              # e^N transition criterion
)

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 1 — NACA KRIGING                                                 ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Objective function for NACA optimization
const OBJECTIVE_NACA = :clcd                      # maximize Cl/Cd at best α

# NACA 4-digit search bounds
const NACA_NDOE = 20                               # number of initial Kriging DOE samples
const NACA_LB1  = 2                                # lower bound — max camber digit
const NACA_UB1  = 9                                # upper bound — max camber digit
const NACA_LB2  = 2                                # lower bound — camber position digit
const NACA_UB2  = 8                                # upper bound — camber position digit

# NACA thickness
const NACA_THICKNESS = 0.234                       # t/c for the NACA XX-YY-ZZ family (last two digits)

# Kriging error convergence
const NACA_TOL_PERCENTAGE = 0.5                    # max squared-error tolerance [% of mean objective]

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SECTION 2 — CST OPTIMIZER                                                ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Objective function for CST optimization
 const OBJECTIVE_CST = :cp_max                     # maximize power coefficient
#const OBJECTIVE_CST = :clcd                       # maximize Cl/Cd at best α

# ── DOE + COBYLA settings ─────────────────────────────────────────────────
const N_DOE     = 60                               # LHS feasible samples per phase
const N_TOP     = 4                                # top-N DOE points refined by COBYLA
const MAX_EVALS = 100                              # max COBYLA evaluations per restart

const ETA_ROOT  = 0.40                             # spanwise station for root evaluation (0 = hub, 1 = tip)
const ETA_TIP   = 0.90                             # spanwise station for tip evaluation

const FTOL_REL      = 1e-3                         # COBYLA relative function tolerance
const XTOL_REL      = 1e-3                         # COBYLA relative parameter tolerance
const COBYLA_RHOBEG = 0.25                         # initial step size as fraction of mean bound width
                                                   #   0.05 = conservative, 0.15 = balanced, 0.30 = aggressive

# ── Penalties & caps ──────────────────────────────────────────────────────
# Soft-stall penalty — penalizes non-monotonic Cl before stall
const STALL_PENALTY_K = 94.0                       # sensitivity (0 = disabled, 20 = lenient, 50 = balanced, 100 = strict)

# Cm hard cap — reject profiles with pitching moment too negative (too much camber)
const CM_CAP_TIP  = -0.19                          # tip  Cm threshold (set to -Inf to disable)
const CM_CAP_ROOT = -0.22                          # root Cm threshold (set to -Inf to disable)
