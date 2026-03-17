
# Constructor for boundary layer problems.
# Multiple dispatch on FlowProperties, Geometry, and BoundaryCondition subtypes
# selects the correct physics and solver.

# ---------------------------------------------------------------------------
# Internal dispatch helpers
# ---------------------------------------------------------------------------

function _solve_incompressible(flow::IncompressibleFlow, geometry::Geometry,
                                bc::BoundaryCondition, Re; kw...)
    η, f, fp, fpp = solve_matrix(flow, geometry, bc; kw...)
    return IncompressibleBoundaryLayer(flow, geometry, bc, Re, η, f, fp, fpp)
end

function _solve_compressible(flow::CompressibleFlow, geometry::Geometry,
                              bc::BoundaryCondition, Re; kw...)
    fluid = flow.fluid
    U_e   = flow.M∞ * sqrt(fluid.γ * fluid.R * flow.T∞)
    η, f, fp, fpp, θ, θp = solve_matrix(flow, geometry, bc; kw...)
    return CompressibleBoundaryLayer(flow, geometry, bc, Re, η, f, fp, fpp, θ, θp, U_e)
end

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    BoundaryLayer(flow, geometry, Re, bcs; kwargs...) -> CanonicalFlow

Solve for a self-similar laminar boundary layer. The specific problem is
determined by multiple dispatch on the types of `flow`, `geometry`, and the
wall boundary condition.

# Arguments
- `flow`: `IncompressibleFlow(U_e, fluid)` or `CompressibleFlow(M∞, T∞, fluid)`
- `geometry`: `FlatPlate()`, `Wedge(half_angle)`, or `Cone(half_angle)`
- `Re`: Reynolds number (U_e L / ν_e for incompressible; ρ_e U_e L / μ_e for compressible)
- `bcs`: tuple of `BoundaryCondition`; the first element is the wall BC

# Example
```julia
flow = CompressibleFlow(5.0, 300.0, Air())
bl   = BoundaryLayer(flow, FlatPlate(), 1e6, (AdiabaticBC(),))
```
"""
BoundaryLayer(flow::IncompressibleFlow, geometry::Geometry, Re::Real,
              bcs; kw...) =
    _solve_incompressible(flow, geometry, first(bcs), Re; kw...)

BoundaryLayer(flow::CompressibleFlow, geometry::Geometry, Re::Real,
              bcs; kw...) =
    _solve_compressible(flow, geometry, first(bcs), Re; kw...)
