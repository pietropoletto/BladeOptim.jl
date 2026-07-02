# BladeOptim.jl

Wind turbine blade optimization framework in Julia. Combines CST airfoil parameterization, XFOIL aerodynamic analysis, Kriging surrogate modeling, and BEM theory (via CCBlade.jl) to design high-performance blade geometries.

## Features

- **Dual-airfoil CST optimizer** — separate root and tip airfoil optimization with spanwise blending
- **NACA 4-digit Kriging** — surrogate-based exploration of the NACA design space with adaptive error minimization
- **Full BEM analysis** — Betz-optimal chord and twist distributions, Prandtl tip-loss correction, Viterna extrapolation for post-stall
- **3D blade visualization** — interactive GLMakie rendering with CST-blended sections, root cylinder transition, and spanwise color mapping
- **Constraint handling** — Cm hard caps, soft-stall penalty, geometric feasibility checks (thickness, camber, tail crossover)

## Project Structure

```
BladeOptim.jl/
├── src/
│   ├── BladeOptimFull.jl      # Full optimizer entry point
│   ├── BladeOptimNACA.jl      # NACA Kriging standalone entry point
│   ├── Config.jl              # Structs, parameters, Reynolds computation
│   ├── input_turbine.jl       # Turbine geometry & operating conditions
│   ├── input_optimizer.jl     # Optimizer settings (NACA + CST)
│   ├── BladeGeom/             # Chord, twist, taper, thickness distribution
│   ├── CST/                   # CST airfoil parameterization
│   ├── XFOIL/                 # XfoilRunner, Conversions, Viterna, AlphaOpt
│   ├── CCBlade/               # BEM performance & loads analysis
│   ├── OptimCST/                 # ObjFunction, SequentialEval (DOE + COBYLA)
│   ├── OptimNACA/             # NACACord, ObjFunctionNACA (Kriging)
│   ├── Plotting/              # NACAPlots, OptimizerPlots
│   └── BladeVisual/           # GLMakie 3D blade rendering
├── examples/
│   ├── naca_kriging.jl        # NACA design space exploration
│   └── optimizer.jl           # Full dual-airfoil optimization
├── tools/
│   └── Coord_to_cst.jl        # Fit CST coefficients to a .dat airfoil
│   └── install_packages.jl    # Pre install all the packages needed 
├── Results_cst/               # CST optimizer outputs
└── Results_naca/              # NACA Kriging outputs
```



## Quick Start

### 1. NACA Kriging

Explore the NACA 4-digit design space to find the best camber/position combination:

```julia
cd("examples")
include("naca_kriging.jl")
```

Configure search bounds and objective in `src/input_optimizer.jl`:
```julia
const OBJECTIVE_NACA = :clcd      # or :cp_max, :cp_robust
const NACA_LB1 = 2                # camber lower bound
const NACA_UB1 = 9                # camber upper bound
```

### 2. Dual-Airfoil Optimizer

Run the full two-phase optimization (tip → root):

```julia
cd("examples")
include("optimizer.jl")
```

This generates optimized CST coefficients, performance plots, an Excel report, and a 3D blade visualization in `Results_cst/`.

### 3. Airfoil Fitting

Convert a `.dat` airfoil to CST coefficients:

```bash
julia tools/Coord_to_cst.jl path/to/airfoil.dat
```

## Configuration

All parameters are in two input files:

- **`src/input_turbine.jl`** — rotor geometry (radius, blades, TSR), airfoil CST coefficients, thickness distribution, twist offset
- **`src/input_optimizer.jl`** — objective function, XFOIL settings, DOE/COBYLA parameters, Kriging bounds, penalty tuning

## How It Works

The optimizer runs in two sequential phases:

1. **Phase 1 (Tip)** — optimizes the tip airfoil CST coefficients while holding the root polar fixed. Generates N DOE samples via Latin Hypercube, evaluates each through XFOIL → Viterna → CCBlade, then refines the top candidates with COBYLA.

2. **Phase 2 (Root)** — optimizes the root airfoil using the best tip polar from Phase 1, same DOE + COBYLA strategy.

Spanwise blending interpolates between root (η ≤ 0.40) and tip (η ≥ 0.90) polars at each BEM section. Chord follows the Betz optimum with Prandtl tip-loss correction and optional hub taper.

