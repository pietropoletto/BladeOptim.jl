module ObjFunction

using DataFrames
using ..Config: OBJECTIVE_CST, STALL_PENALTY_K, CM_CAP_TIP, CM_CAP_ROOT

export OBJECTIVE
export _interp_polar, soft_stall_penalty
export CM_CAP_TIP, CM_CAP_ROOT

const OBJECTIVE = OBJECTIVE_CST

"""
soft_stall_penalty(df; k) → Float64 in [0, 1]

Penalty for a non-monotonic CL curve in the pre-stall region.
"""
function soft_stall_penalty(df::DataFrames.DataFrame; k::Float64=STALL_PENALTY_K)
    k <= 0.0 && return 0.0

    alpha = df[:, :alpha]
    cl    = df[:, :cl]
    n     = length(alpha)
    n < 3 && return 0.0

    idx_clmax = argmax(cl)

    neg_grad_sum = 0.0
    n_neg = 0
    for i in 2:idx_clmax
        dalpha = alpha[i] - alpha[i-1]
        dalpha <= 0.0 && continue
        dcl_dalpha = (cl[i] - cl[i-1]) / dalpha
        if dcl_dalpha < 0.0
            neg_grad_sum += dcl_dalpha^2
            n_neg += 1
        end
    end

    n_neg == 0 && return 0.0

    penalty = 1.0 - exp(-k * neg_grad_sum)
    return clamp(penalty, 0.0, 1.0)
end

"""Interpolate a polar DataFrame onto a common alpha grid."""
function _interp_polar(df::DataFrames.DataFrame, alpha_grid::Vector{Float64})
    cl_out = zeros(length(alpha_grid))
    cd_out = zeros(length(alpha_grid))
    for (i, a) in enumerate(alpha_grid)
        if a <= df.alpha[1]
            cl_out[i] = df.cl[1];  cd_out[i] = df.cd[1]
        elseif a >= df.alpha[end]
            cl_out[i] = df.cl[end]; cd_out[i] = df.cd[end]
        else
            idx = findlast(df.alpha .<= a)
            if idx === nothing || idx == DataFrames.nrow(df)
                cl_out[i] = df.cl[end]; cd_out[i] = df.cd[end]
            else
                t = (a - df.alpha[idx]) / (df.alpha[idx+1] - df.alpha[idx])
                cl_out[i] = df.cl[idx] + t * (df.cl[idx+1] - df.cl[idx])
                cd_out[i] = df.cd[idx] + t * (df.cd[idx+1] - df.cd[idx])
            end
        end
    end
    return cl_out, cd_out
end

end # module
