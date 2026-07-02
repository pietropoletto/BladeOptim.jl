module Viterna

using DataFrames
using CSV
using CCBlade

export apply_viterna_extrapolation

function apply_viterna_extrapolation(df::DataFrame, r::AbstractVector,
                                    chord::AbstractVector, Rtip::Real,
                                    csv_out::Union{AbstractString, Nothing} = nothing)

    r_75 = 0.75 * Rtip
    idx_75 = argmin(abs.(r .- r_75))
    cr75 = chord[idx_75] / r_75

    alpha_0 = df.alpha .* (π/180)
    cl_0 = df.cl
    cd_0 = df.cd

    alpha_ext, cl_ext, cd_ext = CCBlade.viterna(alpha_0, cl_0, cd_0, cr75)

    alpha_ext_deg = alpha_ext .* (180/π)

    df_ext = DataFrame(alpha=alpha_ext_deg, cl=cl_ext, cd=cd_ext)

    if csv_out !== nothing
        CSV.write(csv_out, df_ext; writeheader=false)
    end

    return df_ext, cr75
end

end