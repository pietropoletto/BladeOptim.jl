# Set working directory if needed: cd("/path/to/WT_opt")
println("Current working directory: ", pwd())

# Clean up leftover .dat files from previous runs
foreach(rm, filter(f -> endswith(f, ".dat"), readdir(".", join=true)))
println("Leftover .dat files removed")

include("../src/BladeOptimNACA.jl")   # standalone: no penalty, no optimizer deps
using .BladeOptimNACA

const PROJECT_ROOT = dirname(@__DIR__)
using DataFrames
using CSV
using Plots
using PlotlyJS
using Surrogates
using Statistics
using Printf

Plots.plotlyjs()


############ BOUNDS & RANGES ###########################################

lower_bound = [BladeOptimNACA.Config.lb1, BladeOptimNACA.Config.lb2]
upper_bound = [BladeOptimNACA.Config.ub1, BladeOptimNACA.Config.ub2]

lb1 = BladeOptimNACA.Config.lb1
ub1 = BladeOptimNACA.Config.ub1
lb2 = BladeOptimNACA.Config.lb2
ub2 = BladeOptimNACA.Config.ub2

x1_range = (lb1:1:ub1)
x2_range = (lb2:1:ub2)

# Grid of evaluation points
x_true_range = [(x1, x2) for x1 in x1_range for x2 in x2_range]

x1s = [xi[1] for xi in x_true_range]
x2s = [xi[2] for xi in x_true_range]

objective_str = string(BladeOptimNACA.ObjFunctionNACA.OBJECTIVE)

#=
############ NACA MAP — load from CSV or compute ###########################################

println("")
println("-----------------------------------")
println("         Building NACA map")
println("-----------------------------------")

results_dir = joinpath(PROJECT_ROOT, "Results_naca")
csv_path    = joinpath(results_dir, "true_vals_$(objective_str)_$(lb1)$(ub1)$(lb2)$(ub2).csv")

if isfile(csv_path)
    println("Found existing CSV for objective ':$objective_str': $(basename(csv_path))")
    println("Loading cached values (no XFOIL/CCBlade calls)...")
    y_true_df = CSV.read(csv_path, DataFrame)
    y_true    = Matrix(y_true_df)
    println("  → $(size(y_true,1)) × $(size(y_true,2)) valori caricati")
else
    println("No CSV found for objective ':$objective_str' and bounds [$lb1,$ub1]x[$lb2,$ub2].")
    println("Computing true map ($(length(x_true_range)) points)...")
    y_true_vec = BladeOptimNACA.ObjFunctionNACA.obj_func_load.(x_true_range; label="[Mappa]")
    y_true     = reshape(y_true_vec, (length(x2_range), length(x1_range)))
    mkpath(results_dir)
    CSV.write(csv_path, DataFrame(y_true, :auto))
    println("Saved to: $csv_path")
end
=#

############ KRIGING MAP CREATION ###########################################

println("")
println("-----------------------------------")
println("      Building kriging map")
println("-----------------------------------")

x_samples = Surrogates.sample(BladeOptimNACA.Config.nDOE, lower_bound, upper_bound,
                               Surrogates.LatinHypercubeSample())

println("Kriging samples: ", x_samples)

y_samples = BladeOptimNACA.ObjFunctionNACA.obj_func_load.(x_samples; label="[Kriging DOE]")

naca_surrogate = Kriging(x_samples, y_samples, lower_bound, upper_bound, p = [2.0, 2.0])

naca_surrogate_initial = Kriging(deepcopy(x_samples), deepcopy(y_samples),
                                  lower_bound, upper_bound, p = [2.0, 2.0])

x1s_doe = [xi[1] for xi in x_samples]
x2s_doe = [xi[2] for xi in x_samples]



############ KRIGING ERROR MINIMIZATION ###########################################

println("")
println("-----------------------------------")
println("   Minimization error process")
println("-----------------------------------")

tol_percentage = BladeOptimNACA.Config.NACA_TOL_PERCENTAGE
tol_frac = tol_percentage / 100
maxtol   = (tol_frac * mean(y_samples))^2
println("  maxtol = $(round(maxtol, sigdigits=7))  ($(tol_frac*100)% of mean $(round(mean(y_samples), digits=8)))")

naca_surrogate, x_samples, y_samples, outlier_indices = BladeOptimNACA.ObjFunctionNACA.krig_error_minim(
    x_true_range, naca_surrogate, BladeOptimNACA.ObjFunctionNACA.obj_func_load,
    x_samples, y_samples, lower_bound, upper_bound, maxtol)

sample_x1 = [xi[1] for xi in x_samples]
sample_x2 = [xi[2] for xi in x_samples]

great_array        = hcat(sample_x1, sample_x2, y_samples)
great_array_sorted = sortslices(great_array, dims=1, by=x->x[1])

xvals = sort(unique(great_array_sorted[:,1]))
yvals = sort(unique(great_array_sorted[:,2]))


############ SAVE BEST AIRFOIL ###########################################

println("")
println("-----------------------------------")
println("      Saving best airfoil")
println("-----------------------------------")

x1_fine = range(lb1, ub1, length=500)
x2_fine = range(lb2, ub2, length=500)
y_kriging_fine = [naca_surrogate([x1, x2]) for x2 in x2_fine, x1 in x1_fine]
best_ci_krig    = argmax(y_kriging_fine)
best_x_krig     = (x1_fine[best_ci_krig[2]], x2_fine[best_ci_krig[1]])
best_value_krig = y_kriging_fine[best_ci_krig]

clean_indices   = [i for i in eachindex(y_samples) if i ∉ outlier_indices]
best_sample_idx = clean_indices[argmax(y_samples[clean_indices])]
best_sample_val = y_samples[best_sample_idx]

if best_sample_val >= best_value_krig
    println("  ℹ Best real sample ($(round(best_sample_val, digits=4))) ≥ Kriging optimum ($(round(best_value_krig, digits=4))) → using sample")
    best_x     = x_samples[best_sample_idx]
    best_value = best_sample_val
else
    best_x     = best_x_krig
    best_value = best_value_krig
end

println("Best NACA (kriging): camber=$(round(best_x[1], digits=3))%, position=$(round(best_x[2], digits=3))0%")
println("Best $(objective_str) = $(round(best_value, digits=6))")

n_points = 25
best_x_coords, best_y_coords = BladeOptimNACA.NACACord.nacaXXXX_coordinates(best_x, n_points)
best_coord = hcat(best_x_coords, best_y_coords)

run_dir = joinpath(PROJECT_ROOT, "Results_naca")
mkpath(run_dir)

best_dat_path = joinpath(run_dir, "best_naca_airfoil.dat")
open(best_dat_path, "w") do io
    println(io, "NACA$(round(best_x[1], digits=2))$(round(best_x[2], digits=2))12")
    for i in 1:size(best_coord, 1)
        @printf(io, "%.6f  %.6f\n", best_coord[i, 1], best_coord[i, 2])
    end
end
println("Coordinates saved to: $best_dat_path")

best_summary_path = joinpath(run_dir, "best_naca_summary.csv")
CSV.write(best_summary_path, DataFrame(
    camber     = [best_x[1]],
    camber_pos = [best_x[2]],
    objective  = [objective_str],
    value      = [best_value]
))
println("Summary saved to: $best_summary_path")


############ PLOTS ###########################################

obj_label = Dict("cp_max" => "Cp_max", "clcd" => "Cl/Cd", "cp_robust" => "Cp_robust")[objective_str]

BladeOptimNACA.NACAPlots.plot_naca_airfoil(
    naca_surrogate, naca_surrogate_initial,
    x1_range, x2_range,
    x1s_doe, x2s_doe,
    sample_x1, sample_x2,
    xvals, yvals,
    best_x, best_value,
    obj_label,
    run_dir,
    x_true_range,
    nothing)  