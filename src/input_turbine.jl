# =============================================================================
#   INPUT — TURBINE PARAMETERS
#
#   All geometric and operating parameters for the turbine.
#   This file is read by Config.jl — no need to import it elsewhere.
# =============================================================================

# =============================================================================
# THICKNESS
# =============================================================================
const DT_TIP          = 0.13                       # tip airfoil t/c ratio
const DT_ROOT_FACTOR  = 1.8                        # root = tip × factor
const DT_ROOT         = DT_TIP * DT_ROOT_FACTOR    # root airfoil t/c ratio (computed)
const DT_EXPONENT     = 1.5                        # spanwise t/c distribution exponent (1 = linear, 1.5 = rapid at root)

# =============================================================================
# CST AIRFOIL COEFFICIENTS
# =============================================================================

const AIRFOIL_ROOT = (
    wu = [0.446, 0.43, 0.701],
    wl = [-0.233, -0.18, 0.035, 0.084],
    dz = 0.004,                                    # trailing edge gap half-thickness
    dt = DT_ROOT,                                  # design t/c ratio
    N  = 220                                       # number of coordinate points
)

const AIRFOIL_TIP = (
    wu = [0.261, 0.339, 0.629],
    wl = [-0.127, 0.006, 0.117, 0.337],
    dz = 0.004,                                    # trailing edge gap half-thickness
    dt = DT_TIP,                                   # design t/c ratio
    N  = 220                                       # number of coordinate points
)

# =============================================================================
# ROOT CHORD TAPER
# =============================================================================
const ROOT_TAPER_ENABLED = true                    # enable Betz chord tapering near the hub
const CHORD_ROOT_TARGET  = 0.26                    # target chord [m] at the hub section
const ROOT_TAPER_ETA     = 0.07                    # spanwise fraction where taper blends out (0 = hub, 1 = tip)

# =============================================================================
# TWIST OFFSET
# =============================================================================
#   DELTA_ALPHA = 0.0 → pure Betz twist (pitch-controlled)
#   DELTA_ALPHA > 0   → increases root α → less twist → earlier root stall (stall-regulated)
#   DELTA_ALPHA < 0   → decreases root α → more twist
const DELTA_ALPHA = 0.0                            # geometric offset on α_root [deg]

# =============================================================================
# TURBINE OPERATING PARAMETERS
# =============================================================================
const TURBINE = (
    Rhub      = 0.35,                               # hub radius [m]
    Rtip      = 6.0,                                # tip radius [m]
    B         = 3,                                   # number of blades
    n_sections = 40,                                 # radial BEM sections
    precone   = 0.0 * π/180,                         # precone angle [rad]
    yaw       = 0.0,                                 # yaw misalignment [rad]
    tilt      = 0.0 * π/180,                         # shaft tilt [rad]
    pitch     = 0.0 * π/180,                         # blade pitch [rad]
    V_inf     = 10.0,                                # design wind speed [m/s]
    TSR       = 6,                                   # design tip-speed ratio
    rho       = 1.225,                               # air density [kg/m³]
    hubHt     = 12.0,                                # hub height [m]
    shearExp  = 0.2,                                 # wind shear exponent (power law)
    azimuth   = [0.0, 90.0, 180.0, 270.0] * π/180   # azimuthal positions for loads [rad]
)

# =============================================================================
# FLUID PROPERTIES (ISA standard atmosphere)
# =============================================================================
const AIR_RHO = 1.225                              # air density @ 15 °C, 1 atm [kg/m³]
const AIR_MU  = 1.789e-5                           # dynamic viscosity @ 15 °C [Pa·s]
