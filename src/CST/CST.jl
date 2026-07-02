module CST
export airfoil_thickness, CST_airfoil, ClassShape

"Thickness scaling"
function airfoil_thickness(coord::AbstractMatrix{<:Real}, dt::Real)
    x = coord[:, 1]
    y = coord[:, 2]

    split_ind = findfirst(==(0.0), y)
    if split_ind === nothing
        split_ind = argmin(abs.(y))
    end

    upper_y = y[1:split_ind]
    lower_y = y[split_ind:end]

    L = min(length(upper_y), length(lower_y))
    airfoil_thk = maximum(upper_y[1:L] .- lower_y[1:L])

    scale_factor = dt / airfoil_thk
    scaled_y = y .* scale_factor

    return hcat(x, scaled_y)
end

"Generate coordinatet using CTS"
function CST_airfoil(wu::AbstractVector{<:Real},
                     wl::AbstractVector{<:Real},
                     dz::Real, N::Integer)

    N1 = 0.5
    N2 = 1.0

    x    = ones(Float64, N+1)
    yall = zeros(Float64, N+1)
    zeta = zeros(Float64, N+1)

    # cosine clustering toward leading/trailing edge
    for i in 1:(N + 1)
        zeta[i] = 2π / N * (i - 1)
        x[i]    = 0.5 * (cos(zeta[i]) + 1)
    end

    zerind = findfirst(==(0.0), x)
    if zerind === nothing
        zerind = argmin(abs.(x))
    end

    xu = x[1:zerind-1]
    xl = x[zerind:end]

    yu = ClassShape(wu, xu, N1, N2,  dz)
    yl = ClassShape(wl, xl, N1, N2, -dz)

    yall = vcat(yu, yl)
    return hcat(x, yall)
end

"Class-Shape function CST (Class + Shape)"
function ClassShape(w::AbstractVector{<:Real},
                    x::AbstractVector{<:Real},
                    N1::Real, N2::Real, dz::Real)

    C = [xi^N1 * (1 - xi)^N2 for xi in x]

    n = length(w) - 1
    K = [binomial(n, j-1) for j in 1:(n+1)]

    S = zeros(Float64, length(x))
    for i in eachindex(x)
        xi = x[i]
        for j in 1:(n+1)
            S[i] += w[j] * K[j] * xi^(j-1) * (1 - xi)^(n-(j-1))
        end
    end

    y = [C[i]*S[i] + x[i]*dz for i in eachindex(x)]
    return y
end

end