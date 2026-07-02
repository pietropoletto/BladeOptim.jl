module BladeOptimFull

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
include("CST/CST.jl")

# ── XFOIL ──
include("XFOIL/XfoilRunner.jl")
include("XFOIL/Conversions.jl")
include("XFOIL/Viterna.jl")
include("XFOIL/AlphaOpt.jl")

# ── CCBlade ──
include("CCBlade/CCBladeAnalysis.jl")

# ── Plotting ──
include("Plotting/OptimizerPlots.jl")
include("Plotting/NACAPlots.jl")

# ── Optimization ──
include("OptimCST/ObjFunction.jl")
include("OptimCST/SequentialEval.jl")

# ── 3D Blade Visualization ──
include("BladeVisual/BladeVisual.jl")

# ── Using ──
using .Config
using .BladeGeom
using .CST
using .XfoilRunner
using .Conversions
using .Viterna
using .AlphaOpt
using .CCBladeAnalysis
using .OptimizerPlots
using .NACAPlots
using .ObjFunction
using .SequentialEval
using .BladeVisual

# ── Exports ──
export Config, BladeGeom, CST
export XfoilRunner, Conversions, Viterna, AlphaOpt
export CCBladeAnalysis
export OptimizerPlots, NACAPlots
export ObjFunction
export SequentialEval
export BladeVisual

end # module
