
module AlphaOpt

using DataFrames

export find_optimal_alpha

function find_optimal_alpha(df::DataFrame)
    alpha = df[:, :alpha]
    cl    = df[:, :cl]
    cd    = df[:, :cd]

    ratio = (cl) ./ cd
    idx_max = argmax(ratio)
    
    # Find optimal alpha
    alpha_target = alpha[idx_max] - 0.0
    
    # Find the index closest to alpha_target
    idx = argmin(abs.(alpha .- alpha_target))

    cm_val = hasproperty(df, :cm) ? df.cm[idx] : NaN
    return alpha[idx], cl[idx], cd[idx], cm_val
end

end