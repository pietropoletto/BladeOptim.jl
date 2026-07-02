module BladeOptimNACA

# =============================================================================
#  Standalone module for NACA Kriging.
#  Does NOT include: SequentialEval, OptimizerPlots, BladeVisual, CST.
#  Loads only dependencies needed for NACA 4-digit kriging.
# =============================================================================

using Plots
using DataFrames
using CSV
using CCBlade
using Printf
using FLOWMath
using DelimitedFiles
using Surrogates
using PlotlyJS

gr()

# ── Core ──
include("BladeGeom/BladeGeom.jl")
include("Config.jl")

# ── XFOIL ──
include("XFOIL/XfoilRunner.jl")
include("XFOIL/Conversions.jl")
include("XFOIL/Viterna.jl")
include("XFOIL/AlphaOpt.jl")

# ── CCBlade ──
include("CCBlade/CCBladeAnalysis.jl")

# ── Plotting NACA ──
include("Plotting/NACAPlots.jl")

# ── NACA Kriging modules (standalone, no penalty) ──
include("OptimNACA/NACACord.jl")
include("OptimNACA/ObjFunctionNACA.jl")

# ── Using ──
using .Config
using .BladeGeom
using .XfoilRunner
using .Conversions
using .Viterna
using .AlphaOpt
using .CCBladeAnalysis
using .NACAPlots
using .NACACord
using .ObjFunctionNACA

# ── Exports ──
export Config, BladeGeom
export XfoilRunner, Conversions, Viterna, AlphaOpt
export CCBladeAnalysis
export NACAPlots
export NACACord, ObjFunctionNACA

end # module
