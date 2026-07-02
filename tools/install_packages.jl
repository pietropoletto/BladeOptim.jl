# =============================================================================
#  install_packages.jl — Install all dependencies for BladeOptim.jl
#
#  Usage:
#    julia tools/install_packages.jl
# =============================================================================

using Pkg

packages = [
    # ── Core BEM & aerodynamics ──
    "CCBlade",
    "Xfoil",
    "FLOWMath",

    # ── Optimization ──
    "NLopt",
    "Surrogates",
    "Optim",

    # ── Data ──
    "DataFrames",
    "CSV",
    "XLSX",

    # ── Plotting ──
    "Plots",
    "PlotlyJS",
    "GLMakie",
    "Colors",

    # ── Stdlib (usually bundled, but just in case) ──
    "Printf",
    "Statistics",
    "Random",
    "LinearAlgebra",
    "DelimitedFiles",
    "Dates",
]

println("=" ^ 60)
println("  BladeOptim.jl — Installing $(length(packages)) packages")
println("=" ^ 60)

for pkg in packages
    print("  $pkg ... ")
    try
        Pkg.add(pkg)
        println("✓")
    catch e
        println("✗  ($e)")
    end
end

println()
println("=" ^ 60)
println("  Done. You can now run the examples.")
println("=" ^ 60)