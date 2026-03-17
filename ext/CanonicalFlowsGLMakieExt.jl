module CanonicalFlowsGLMakieExt

using CanonicalFlows
import GLMakie

# ---------------------------------------------------------------------------
# Similarity-space helpers (used by plot(bl, vars))
#
# Map each variable symbol to its legend label and normalised profile on bl.η.
# All values are dimensionless fractions of the freestream quantity, using the
# normalisation already present in the similarity solution:
#   :streamwise_vel  → u/U_e  = f′(η)               (both flow types)
#   :temperature     → T/T_e  = θ(η)                 (compressible only)
#   :density         → ρ/ρ_e  = 1/θ(η)               (compressible only; perfect gas, const p)
#   :pressure        → p/p_e  = 1 everywhere          (dp/dy = 0 across the BL)
# ---------------------------------------------------------------------------

function _normalized_profile(bl::CanonicalFlow, var::Symbol)
    if var === :streamwise_vel
        return ("u / U_e", bl.fp)

    elseif var === :temperature
        bl isa CompressibleBoundaryLayer ||
            error(":temperature is only available for CompressibleBoundaryLayer")
        return ("T / T_e", bl.θ)

    elseif var === :density
        bl isa CompressibleBoundaryLayer ||
            error(":density is only available for CompressibleBoundaryLayer")
        return ("ρ / ρ_e", 1.0 ./ bl.θ)

    elseif var === :pressure
        return ("p / p_e", ones(eltype(bl.η), length(bl.η)))

    elseif var === :wallnormal_vel
        error(":wallnormal_vel depends on streamwise position x and cannot be plotted " *
              "in similarity coordinates. Use plot(bl, x, vars) instead.")

    else
        error("Unknown variable: $var. Supported: " *
              ":streamwise_vel, :temperature, :density, :pressure")
    end
end

# ---------------------------------------------------------------------------
# Physical-space helpers (used by plot(bl, x, vars, N))
# ---------------------------------------------------------------------------

# Invert the coordinate transform: given a similarity coordinate η_target,
# return the physical wall-normal position y/L at streamwise location x.
#
# Incompressible: η = y · scale  →  y = η / scale   (exact, linear)
# Compressible:   y · scale = Y(η) = ∫₀^η θ dη′    →  y = Y(η_target) / scale
#   Integral evaluated by the trapezoidal rule, identical to the forward
#   transform in _η_from_xy but stopped at η_target rather than inverted.

function _y_from_eta(bl::IncompressibleBoundaryLayer, x::Real, η_target::Real)
    scale = CanonicalFlows._η_scale(bl.geometry, Float64(bl.Re), x)
    return η_target / scale
end

function _y_from_eta(bl::CompressibleBoundaryLayer, x::Real, η_target::Real)
    η_grid = bl.η
    θ      = bl.θ
    Y      = 0.0
    for i in 2:length(η_grid)
        if η_grid[i] >= η_target
            # Partial interval: integrate only up to η_target
            t     = (η_target - η_grid[i-1]) / (η_grid[i] - η_grid[i-1])
            θ_end = θ[i-1] + t * (θ[i] - θ[i-1])
            Y    += 0.5 * (θ[i-1] + θ_end) * (η_target - η_grid[i-1])
            break
        end
        Y += 0.5 * (θ[i] + θ[i-1]) * (η_grid[i] - η_grid[i-1])
    end
    scale = CanonicalFlows._η_scale(bl.geometry, Float64(bl.Re), x)
    return Y / scale
end

# Return (legend label, freestream scalar) for normalising a dimensional profile
# produced by evaluate_profile.  Dividing the profile vector by the scalar gives
# a fraction of the freestream value.

function _freestream_normalizer(bl::IncompressibleBoundaryLayer, var::Symbol)
    U_e = bl.flow_prop.U_e
    if var === :streamwise_vel
        return ("u / U_e", U_e)
    elseif var === :wallnormal_vel
        return ("v / U_e", U_e)
    elseif var === :pressure
        return ("p / p_e", 1.0)   # evaluate_profile returns 0 (gauge) → ratio is 0
    elseif var === :temperature
        error(":temperature is not solved for IncompressibleBoundaryLayer")
    elseif var === :density
        error(":density is not solved for IncompressibleBoundaryLayer")
    else
        error("Unknown variable: $var")
    end
end

function _freestream_normalizer(bl::CompressibleBoundaryLayer, var::Symbol)
    T_e = bl.flow_prop.T∞
    μ_e = viscosity(bl.flow_prop.fluid, T_e)
    ρ_e = bl.Re * μ_e / bl.U_e
    p_e = ρ_e * bl.flow_prop.fluid.R * T_e

    if var === :streamwise_vel
        return ("u / U_e", bl.U_e)
    elseif var === :wallnormal_vel
        return ("v / U_e", bl.U_e)
    elseif var === :temperature
        return ("T / T_e", T_e)
    elseif var === :pressure
        return ("p / p_e", p_e)
    elseif var === :density
        return ("ρ / ρ_e", ρ_e)
    else
        error("Unknown variable: $var")
    end
end

# ---------------------------------------------------------------------------
# plot(bl, vars) → GLMakie.Figure
#
# Plots each variable as a fraction of its freestream value against η,
# from the wall (η = 0) to 1.2 × the boundary layer edge location.
# ---------------------------------------------------------------------------

function GLMakie.plot(bl::CanonicalFlow, vars::NTuple{N, Symbol}) where {N}
    η_edge    = boundary_layer_edge(bl)
    η_max_plt = 1.2 * η_edge

    η  = bl.η
    i1 = clamp(searchsortedlast(η, η_max_plt), 1, length(η))
    η_plt = η[1:i1]

    fig = GLMakie.Figure()
    ax  = GLMakie.Axis(fig[1, 1];
                       xlabel = "Fraction of freestream value",
                       ylabel = "η (similarity coordinate)")

    for var in vars
        label, data = _normalized_profile(bl, var)
        GLMakie.lines!(ax, data[1:i1], η_plt; label = label)
    end

    GLMakie.axislegend(ax)
    GLMakie.ylims!(ax, 0.0, η_max_plt)

    return fig
end

# ---------------------------------------------------------------------------
# plot(bl, x, vars, N) → GLMakie.Figure
#
# Plots each variable as a fraction of its freestream value against the
# physical wall-normal coordinate y/L at streamwise position x, from the
# wall (y = 0) to 1.2 × the physical boundary layer thickness at that x.
# N controls the number of wall-normal sample points (default 50).
# ---------------------------------------------------------------------------

function GLMakie.plot(bl::CanonicalFlow, x::Real,
                      vars::NTuple{M, Symbol}, N::Int = 50) where {M}
    η_edge  = boundary_layer_edge(bl)
    y_max   = _y_from_eta(bl, x, 1.2 * η_edge)
    y_range = LinRange(0.0, y_max, N)

    data_tuple = evaluate_profile(bl, x, y_range, vars)

    fig = GLMakie.Figure()
    ax  = GLMakie.Axis(fig[1, 1];
                       xlabel = "Fraction of freestream value",
                       ylabel = "y / L")

    for (i, var) in enumerate(vars)
        label, fs_val = _freestream_normalizer(bl, var)
        GLMakie.lines!(ax, data_tuple[i] ./ fs_val, collect(y_range); label = label)
    end

    GLMakie.axislegend(ax)
    GLMakie.ylims!(ax, 0.0, y_max)

    return fig
end

end
