module XfoilRunner

# CHANGED: replaced external XFOIL executable (XFOILExec module) with
#          Xfoil.jl built-in Julia package. No more file I/O, subprocess
#          spawning, or polar.out parsing. Returns a DataFrame directly.

using Xfoil
using DataFrames
using Plots

export run_xfoil

"""
    run_xfoil(x_coords, y_coords; Re, Mach, alpha_start, alpha_end, alpha_step, iter, ncrit)

Run XFOIL analysis using the built-in Xfoil.jl package.
Returns a DataFrame with columns: alpha, cl, cd, cdp, cm
(only converged points are returned).
"""
function run_xfoil(x_coords::AbstractVector{<:Real},
                   y_coords::AbstractVector{<:Real};
                   Re::Real          = 300000.0,
                   Mach::Real        = 0.0,
                   alpha_start::Real = -2.0,
                   alpha_end::Real   = 18,
                   alpha_step::Real  = 0.5,
                   iter::Int         = 100,
                   ncrit::Real       = 6.0,
                   reinitialize::Bool = true,
                   show_airfoil::Bool = true,
                   timeout::Real     = 15.0)   # ADDED: timeout in seconds (like old XFOIL_exec)

    # ADDED: live plot of the airfoil being analyzed.
    #        Set show_airfoil=false to disable.
    if show_airfoil
        p = Plots.plot(x_coords, y_coords,
                 label="", lw=2, color=:black,
                 aspect_ratio=:equal,
                 title="XFOIL — analyzing airfoil  (Re=$(round(Int, Re)))",
                 xlabel="x/c", ylabel="y/c",
                 size=(700, 250),               # CHANGED: larger plot to fit window
                 margin=5Plots.mm)
        display(p)
    end

    alpha_range = alpha_start:alpha_step:alpha_end

    # ADDED: timeout mechanism similar to the old XFOILExec (15s default).
    #        Xfoil.jl runs in-process so we can't kill it, but we can run it
    #        in a Task and wait with a timeout. If it hangs, we return nothing.
    result_channel = Channel{Any}(1)

    task = @async begin
        try
            cl, cd, cdp, cm, conv = Xfoil.alpha_sweep(
                x_coords, y_coords, alpha_range, Re;
                mach = Mach, iter = iter, ncrit = ncrit,
                reinit = reinitialize
            )
            put!(result_channel, (cl, cd, cdp, cm, conv))
        catch e
            put!(result_channel, e)
        end
    end

    # Wait up to `timeout` seconds
    result = timedwait(() -> isready(result_channel), timeout)

    if result == :timed_out
        println("⚠ XFOIL timeout after $timeout seconds — skipping airfoil")
        return nothing
    end

    data = take!(result_channel)

    if data isa Exception
        println("⚠ XFOIL error: $data — skipping airfoil")
        return nothing
    end

    cl, cd, cdp, cm, converged = data
    alpha_vec = collect(alpha_range)

    # Filter only converged points (same logic as before: discard failed alphas)
    mask = converged .== true
    if sum(mask) == 0
        println("⚠ XFOIL: no converged points in alpha range [$alpha_start, $alpha_end] — skipping airfoil")
        return nothing   # CHANGED: return nothing instead of error() to not crash the kriging loop
    end

    df = DataFrame(
        alpha = alpha_vec[mask],
        cl    = cl[mask],
        cd    = cd[mask],
        cdp   = cdp[mask],
        cm    = cm[mask]
    )

    return df
end

end # module