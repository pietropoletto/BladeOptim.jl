module Conversions

using Printf
using DataFrames
using CSV

# CHANGED: removed xfoil2csv function — no longer needed because Xfoil.jl
#          returns cl, cd arrays directly (no more polar.out file to parse).

export write_ccblade_airfoil, write_xfoil_dat

"Write a .dat file for XFOIL"
function write_xfoil_dat(path::AbstractString, coord::AbstractMatrix{<:Real}; 
                        title::AbstractString="CST_AIRFOIL")
    dirp = dirname(path)
    isdir(dirp) || mkpath(dirp)

    open(path, "w") do io
        println(io, title)
        @inbounds for i in 1:size(coord,1)
            @printf(io, "%.6f  %.6f\n", coord[i,1], coord[i,2])
        end
    end
    return path
end

"New CCBlade format requested"
function write_ccblade_airfoil(df::DataFrame, af_file::AbstractString, 
                               Re::Real, Mach::Real)
    dirp = dirname(af_file)
    isdir(dirp) || mkpath(dirp)

    open(af_file, "w") do io
        println(io, "Airfoil from XFOIL polar")
        @printf(io, "%.2f\n", Re)
        @printf(io, "%.4f\n", Mach)

        @inbounds for i in 1:nrow(df)
            @printf(io, "%.4f  %.6f  %.6f\n",
                    df.alpha[i], df.cl[i], df.cd[i])
        end
    end
    return af_file
end

end