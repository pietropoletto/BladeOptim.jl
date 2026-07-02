module NACACord

using ..Config

function cosine_spacing(n)
    θ = range(0, π, length=n)
    return 0.5 * (1 .- cos.(θ))
end

function nacaXXXX_coordinates(x, n_points)
    m = x[1] / 100
    p = x[2] / 10
    t = Config.NACA_THICKNESS  # <- da Config_geom.jl

    z  = cosine_spacing(n_points)
    yt = 5 * t * (0.2969 * sqrt.(z) .- 0.1260 * z .- 0.3516 * z.^2 .+
                  0.2843 * z.^3 .- 0.1015 * z.^4)

    y_c    = similar(z)
    dyc_dx = similar(z)

    for i in eachindex(z)
        if z[i] < p
            y_c[i]    = (m / p^2) * (2p * z[i] - z[i]^2)
            dyc_dx[i] = (2m / p^2) * (p - z[i])
        else
            y_c[i]    = (m / (1 - p)^2) * ((1 - 2p) + 2p * z[i] - z[i]^2)
            dyc_dx[i] = (2m / (1 - p)^2) * (p - z[i])
        end
    end

    θ  = atan.(dyc_dx)
    xu = z .- yt .* sin.(θ)
    yu = y_c .+ yt .* cos.(θ)
    xl = z .+ yt .* sin.(θ)
    yl = y_c .- yt .* cos.(θ)

    x_coords = vcat(reverse(xu), xl[2:end])
    y_coords = vcat(reverse(yu), yl[2:end])

    return x_coords, y_coords
end

function naca_shape_func(x, n_points)
    x_cor, y_cor = nacaXXXX_coordinates(x, n_points)

    cor_pair = hcat(x_cor, y_cor)

    airname = "NACA$(round(Int,x[1]))$(round(Int,x[2]))$(round(Int, Config.NACA_THICKNESS*100))"
    fname   = joinpath(tempdir(), "$airname.dat")

    open(fname, "w") do f
        write(f, "NACA$(x[1])$(x[2])$(round(Int, Config.NACA_THICKNESS*100)) \n")
        for i in 1:size(cor_pair, 1)
            write(f, "$(cor_pair[i,1]) $(cor_pair[i,2])\n")
        end
    end

    return fname
end

end # module
