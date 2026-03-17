
# Physical-space evaluation of boundary layer solutions.
#
# Coordinate conventions (all lengths non-dimensionalised by reference length L):
#   Incompressible flat plate:   η = y √(Re / x)
#   Incompressible wedge:        η = y √(Re x^(m−1))  where m = β/(2−β)
#   Incompressible cone:         η = y √(3 Re / (sin²α x³))   [Mangler transform]
#   Compressible (all bodies):   η defined implicitly via
#                                y * scale(x) = ∫₀^η θ(η′) dη′ = Y(η)

# ---------------------------------------------------------------------------
# Linear interpolation helper
# ---------------------------------------------------------------------------

function _interp(x_data::AbstractVector, y_data::AbstractVector, x::Real)
    n = length(x_data)
    x < x_data[1]   && return y_data[1]
    x > x_data[end] && return y_data[end]
    k = searchsortedlast(x_data, x)
    k == n && return y_data[n]
    t = (x - x_data[k]) / (x_data[k+1] - x_data[k])
    return y_data[k] + t * (y_data[k+1] - y_data[k])
end

# ---------------------------------------------------------------------------
# η scale factor: physical (x) → scalar multiplier so that η = y * _η_scale(...)
# Dispatch on Geometry subtype eliminates isa branching.
# ---------------------------------------------------------------------------

_η_scale(::FlatPlate, Re::Float64, x::Real) = sqrt(Re / x)

function _η_scale(g::Wedge, Re::Float64, x::Real)
    β = _β_parameter(g)
    m = β / (2 - β)          # power-law exponent: U_e ∝ x^m
    return sqrt(Re * x^(m - 1))
end

function _η_scale(g::Cone, Re::Float64, x::Real)
    x_tilde = sin(g.half_angle)^2 * x^3 / 3   # Mangler mapped coordinate
    return sqrt(Re / x_tilde)
end

# ---------------------------------------------------------------------------
# Similarity coordinate: physical (x, y) → η
# ---------------------------------------------------------------------------

function _η_from_xy(bl::IncompressibleBoundaryLayer, x::Real, y::Real)
    return y * _η_scale(bl.geometry, bl.Re, x)
end

function _η_from_xy(bl::CompressibleBoundaryLayer, x::Real, y::Real)
    η_grid = bl.η
    θ      = bl.θ

    # Cumulative integral Y(η) = ∫₀^η θ dη′ (trapezoidal rule)
    n = length(η_grid)
    Y = zeros(n)
    for i in 2:n
        Y[i] = Y[i-1] + 0.5 * (θ[i] + θ[i-1]) * (η_grid[i] - η_grid[i-1])
    end

    target = y * _η_scale(bl.geometry, bl.Re, x)
    return _interp(Y, η_grid, target)
end

# ---------------------------------------------------------------------------
# Physical variable computation
# ---------------------------------------------------------------------------

function _compute_variable(bl::IncompressibleBoundaryLayer, x::Real,
                            y_range, var::Symbol)
    U_e = bl.flow_prop.U_e
    Re  = bl.flow_prop.Re
    out = Vector{Float64}(undef, length(y_range))
    for (k, y) in enumerate(y_range)
        η  = _η_from_xy(bl, x, y)
        fp = _interp(bl.η, bl.fp, η)
        f  = _interp(bl.η, bl.f,  η)

        if var === :streamwise_vel
            out[k] = U_e * fp

        elseif var === :wallnormal_vel
            # v = (U_e / 2√(Re·x)) * (η f′ − f)
            out[k] = U_e / (2 * sqrt(Re * x)) * (η * fp - f)

        elseif var === :pressure
            out[k] = 0.0   # gauge pressure (dp/dy = 0 across BL)

        elseif var === :temperature
            out[k] = 0.0   # not solved in incompressible; return zero

        elseif var === :density
            out[k] = 0.0   # constant in incompressible; return zero

        else
            error("Unknown variable: $var. Supported: " *
                  ":streamwise_vel, :wallnormal_vel, :pressure, :temperature, :density")
        end
    end
    return out
end

function _compute_variable(bl::CompressibleBoundaryLayer, x::Real,
                            y_range, var::Symbol)
    T_e = bl.flow_prop.T∞
    μ_e = viscosity(bl.flow_prop.fluid, T_e)
    ρ_e = bl.Re * μ_e / bl.U_e
    p_e = ρ_e * bl.flow_prop.fluid.R * T_e

    out = Vector{Float64}(undef, length(y_range))
    for (k, y) in enumerate(y_range)
        η  = _η_from_xy(bl, x, y)
        fp = _interp(bl.η, bl.fp, η)
        f  = _interp(bl.η, bl.f,  η)
        θ  = _interp(bl.η, bl.θ,  η)

        if var === :streamwise_vel
            out[k] = bl.U_e * fp

        elseif var === :wallnormal_vel
            ν_e = μ_e / ρ_e
            out[k] = sqrt(ν_e * bl.U_e / (2x)) * θ * (η * fp - f)

        elseif var === :pressure
            out[k] = p_e

        elseif var === :temperature
            out[k] = T_e * θ

        elseif var === :density
            out[k] = ρ_e / θ   # perfect gas, constant pressure

        else
            error("Unknown variable: $var. Supported: " *
                  ":streamwise_vel, :wallnormal_vel, :pressure, :temperature, :density")
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    boundary_layer_edge(bl; tol = 1e-3) -> η_edge

Return the similarity coordinate η at the boundary layer edge, defined as the
first location where the dimensionless streamwise velocity f′(η) = u/U_e
satisfies `1 − f′(η) ≤ tol`.

The criterion is a velocity-deficit threshold: the flow is considered to have
reached the freestream once the remaining deficit below U_e is smaller than
`tol` (default 0.1%, i.e. tol = 1e-3). This is more stringent than the
classical 99% boundary layer thickness (tol = 0.01) and avoids the ambiguity
of derivative-based tests while remaining insensitive to solver noise.

The exact crossing is found by linear interpolation between the last
sub-threshold grid point and the first super-threshold grid point. If the
solution domain ends before the criterion is met, the last grid point is
returned with a warning.

The same `fp` profile (u/U_e) is used for both incompressible and compressible
layers; the Dorodnitsyn transformation preserves this normalisation.
"""
function boundary_layer_edge(bl::CanonicalFlow; tol::Real = 1e-3)
    η  = bl.η
    fp = bl.fp
    threshold = 1.0 - tol

    for i in eachindex(fp)
        if fp[i] >= threshold
            i == 1 && return η[1]
            # linear interpolation to the exact crossing η
            Δfp = fp[i] - fp[i-1]
            Δη  = η[i]  - η[i-1]
            return η[i-1] + (threshold - fp[i-1]) / Δfp * Δη
        end
    end

    @warn "boundary_layer_edge: f′ never reached 1 − tol = $threshold within the " *
          "solution domain (η_max = $(η[end])). Returning η_max."
    return η[end]
end

"""
    evaluate_profile(bl, x, y_range, variables) -> Tuple of Vectors

Evaluate the boundary layer solution at streamwise position `x` and a range of
wall-normal positions `y_range` for each requested physical variable.

# Arguments
- `bl`: a solved `CanonicalFlow` (e.g. the result of `BoundaryLayer(...)`)
- `x`: dimensionless streamwise coordinate (x / L)
- `y_range`: iterable of dimensionless wall-normal positions (y / L)
- `variables`: tuple of symbols, any of:
  - `:streamwise_vel`  — streamwise velocity u (m/s)
  - `:wallnormal_vel`  — wall-normal velocity v (m/s)
  - `:pressure`        — static pressure (Pa); 0 (gauge) for incompressible
  - `:temperature`     — static temperature (K); compressible only
  - `:density`         — density (kg/m³); compressible only

# Returns
A `Tuple` of `Vector{Float64}`, one per requested variable, each of length
`length(y_range)`.
"""
function evaluate_profile(bl::CanonicalFlow, x::Real, y_range,
                           variables::NTuple{N, Symbol}) where {N}
    return Tuple(_compute_variable(bl, x, y_range, var) for var in variables)
end
