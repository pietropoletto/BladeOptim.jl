





using DelimitedFiles
using Plots
using Printf
using Optim
using LinearAlgebra

     dat_file = joinpath(@__DIR__, "coordinates.dat")   # insert .dat file path here (relative or absolute, cross-platform)




# CHANGED: separate number of coefficients for upper and lower surface
n_wu = 3   # number of CST coefficients for upper surface
n_wl = 4    # number of CST coefficients for lower surface

# ============================= FUNZIONI CST BASE =============================

function ClassShape(w::AbstractVector{<:Real}, x::AbstractVector{<:Real}, N1::Real, N2::Real, dz::Real)
    C = [xi^N1 * (1 - xi)^N2 for xi in x]
    n = length(w) - 1
    K = [binomial(n, j - 1) for j in 1:(n + 1)]
    S = zeros(Float64, length(x))
    for i in eachindex(x)
        xi = x[i]
        for j in 1:(n + 1)
            S[i] += w[j] * K[j] * xi^(j - 1) * (1 - xi)^(n - (j - 1))
        end
    end
    y = [C[i] * S[i] + x[i] * dz for i in eachindex(x)]
    return y
end

function CST_airfoil(wu::AbstractVector{<:Real}, wl::AbstractVector{<:Real}, dz::Real, N::Integer)
    N1 = 0.5
    N2 = 1.0
    x = ones(Float64, N + 1)
    yall = zeros(Float64, N + 1)
    zeta = zeros(Float64, N + 1)
    for i in 1:(N + 1)
        zeta[i] = 2π / N * (i - 1)
        x[i] = 0.5 * (cos(zeta[i]) + 1)
    end
    zerind = findfirst(==(0.0), x)
    if zerind === nothing
        zerind = argmin(abs.(x))
    end
    xu = x[1:zerind-1]
    xl = x[zerind:end]
    yu = ClassShape(wu, xu, N1, N2, dz)
    yl = ClassShape(wl, xl, N1, N2, -dz)
    yall = vcat(yu, yl)
    return hcat(x, yall)
end

function read_airfoil_dat(filepath::String)
    lines = readlines(filepath)
    start_idx = 2
    coords = []
    for i in start_idx:length(lines)
        line = strip(lines[i])
        if isempty(line)
            continue
        end
        parts = split(line)
        if length(parts) >= 2
            x = parse(Float64, parts[1])
            y = parse(Float64, parts[2])
            push!(coords, [x, y])
        end
    end
    return reduce(hcat, coords)'
end

function split_airfoil(coord::AbstractMatrix{<:Real})
    x = coord[:, 1]
    y = coord[:, 2]
    le_idx = argmin(x)
    xu = x[1:le_idx]
    yu = y[1:le_idx]
    xl = x[le_idx:end]
    yl = y[le_idx:end]
    upper_sorted = sortperm(xu)
    lower_sorted = sortperm(xl)
    xu = xu[upper_sorted]
    yu = yu[upper_sorted]
    xl = xl[lower_sorted]
    yl = yl[lower_sorted]
    return xu, yu, xl, yl
end

function linear_interp(x::Vector, y::Vector, x_new::Vector)
    y_new = zeros(length(x_new))
    for i in eachindex(x_new)
        xi = x_new[i]
        if xi <= x[1]
            y_new[i] = y[1]
        elseif xi >= x[end]
            y_new[i] = y[end]
        else
            idx = findlast(x .<= xi)
            if idx === nothing || idx == length(x)
                y_new[i] = y[end]
            else
                x1, x2 = x[idx], x[idx+1]
                y1, y2 = y[idx], y[idx+1]
                y_new[i] = y1 + (y2 - y1) * (xi - x1) / (x2 - x1)
            end
        end
    end
    return y_new
end

# CHANGED: separate n_wu and n_wl instead of single n_coeff
function fit_cst_parameters(coord_target::AbstractMatrix{<:Real}; n_wu::Int=2, n_wl::Int=4, N::Int=220)
    xu_target, yu_target, xl_target, yl_target = split_airfoil(coord_target)
    wu_init = zeros(n_wu) .+ 0.1
    wl_init = zeros(n_wl) .+ 0.1
    dz_init = 0.0
    p0 = vcat(wu_init, wl_init, dz_init)
    
    function objective(p)
        wu = p[1:n_wu]
        wl = p[n_wu+1:n_wu+n_wl]
        dz = p[end]
        coord_cst = CST_airfoil(wu, wl, dz, N)
        xu_cst, yu_cst, xl_cst, yl_cst = split_airfoil(coord_cst)
        yu_interp = linear_interp(xu_cst, yu_cst, xu_target)
        yl_interp = linear_interp(xl_cst, yl_cst, xl_target)
        err_upper = norm(yu_interp .- yu_target) / sqrt(length(yu_target))
        err_lower = norm(yl_interp .- yl_target) / sqrt(length(yl_target))
        return err_upper + err_lower
    end
    
    result = Optim.optimize(objective, p0, LBFGS(), Optim.Options(iterations=1000, show_trace=false))
    p_opt = Optim.minimizer(result)
    wu_opt = p_opt[1:n_wu]
    wl_opt = p_opt[n_wu+1:n_wu+n_wl]
    dz_opt = p_opt[end]
    error = Optim.minimum(result)
    return wu_opt, wl_opt, dz_opt, error
end

function write_xfoil_dat(filepath::String, coord::AbstractMatrix{<:Real}; title::String="AIRFOIL")
    open(filepath, "w") do io
        println(io, title)
        for i in 1:size(coord, 1)
            @printf(io, "  %.6f  %.6f\n", coord[i, 1], coord[i, 2])
        end
    end
    println("File salvato: $filepath")
end

# CHANGED: separate n_wu and n_wl instead of single n_coeff
function fit_airfoil_to_matrix(dat_filepath::String; n_wu::Int=2, n_wl::Int=4, N::Int=220)
    coord_target = read_airfoil_dat(dat_filepath)
    wu, wl, dz, error = fit_cst_parameters(coord_target, n_wu=n_wu, n_wl=n_wl, N=N)
    
    # Calcola dt (maximum thickness) dall'airfoil originale
    xu, yu, xl, yl = split_airfoil(coord_target)
    L = min(length(yu), length(yl))
    dt = maximum(yu[1:L] .- yl[1:L])
    
    return wu, wl, dz, dt, error
end

# ===================== ESECUZIONE AUTOMATICA =====================

wu, wl, dz, dt, fit_error = fit_airfoil_to_matrix(dat_file, n_wu=n_wu, n_wl=n_wl, N=220)

wu = round.(wu, digits=3)
wl = round.(wl, digits=3)
dz = round(dz, digits=3)
dt = round(dt, digits=3)

# Output in formato Config.jl
output_text = """
wu = [$(join(wu, ", "))],
wl = [$(join(wl, ", "))],
dz = $dz,
dt = $dt,
"""

# Stampa a schermo
println("\n" * "="^60)
println("CST PARAMETERS (Config.jl format):")
println("  n_wu = $n_wu,  n_wl = $n_wl")
println("  fit error = $(round(fit_error, sigdigits=4))")
println("="^60)
println(output_text)
println("="^60)

# Salva su file (opzionale)
#output_path = joinpath(@__DIR__, "opt-values.txt")
#open(output_path, "w") do io
    #write(io, output_text)
#end

#println("\n✅ Saved to: $output_path")