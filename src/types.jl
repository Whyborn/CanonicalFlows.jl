
"""
    CanonicalFlow

Abstract supertype for all solved canonical flow objects returned by solver constructors
such as [`BoundaryLayer`](@ref).
"""
abstract type CanonicalFlow end

"""
    FlowProperties

Abstract supertype for flow condition specifications. Concrete subtypes:
[`IncompressibleFlow`](@ref), [`CompressibleFlow`](@ref).
"""
abstract type FlowProperties end

"""
    Geometry

Abstract supertype for body geometry. Concrete subtypes:
[`FlatPlate`](@ref), [`Wedge`](@ref), [`Cone`](@ref).
"""
abstract type Geometry end

"""
    BoundaryCondition

Abstract supertype for wall thermal boundary conditions. Concrete subtypes:
[`AdiabaticBC`](@ref), [`IsothermalBC`](@ref).
"""
abstract type BoundaryCondition end

"""
    FluidProperties

Stores fluid thermodynamic and transport properties for a perfect gas.
Fields are concrete `Float64` since these are fixed physical constants.
The Prandtl number is not stored directly; use [`prandtl_number`](@ref) to derive
it from the specific heat, reference viscosity, and thermal conductivity.

# Fields
- `γ`: ratio of specific heats (cp/cv)
- `R`: specific gas constant (J/(kg·K))
- `cp`: specific heat at constant pressure (J/(kg·K))
- `μ_ref`: reference dynamic viscosity (Pa·s) at temperature `T_ref`
- `T_ref`: reference temperature for Sutherland's law (K)
- `S`: Sutherland temperature constant (K)
- `k`: reference thermal conductivity (W/(m·K)) at temperature `T_ref`
"""
struct FluidProperties
    γ::Float64
    R::Float64
    cp::Float64
    μ_ref::Float64
    T_ref::Float64
    S::Float64
    k::Float64
end

"""
    IncompressibleFlow{T} <: FlowProperties

Flow properties for an incompressible boundary layer.

# Fields
- `U_e`: freestream velocity (m/s); used for dimensional output (default 1)
- `fluid`: `FluidProperties` containing the transport and thermodynamic constants
"""
struct IncompressibleFlow{T} <: FlowProperties
    U_e::T
    fluid::FluidProperties
end
IncompressibleFlow(fluid::FluidProperties) = IncompressibleFlow(one(Float64), fluid)

"""
    CompressibleFlow{T} <: FlowProperties

Flow properties for a compressible boundary layer.

# Fields
- `M∞`: freestream Mach number
- `T∞`: freestream (edge) temperature (K)
- `fluid`: `FluidProperties` containing the transport and thermodynamic constants
"""
struct CompressibleFlow{T} <: FlowProperties
    M∞::T
    T∞::T
    fluid::FluidProperties
end

"""
    FlatPlate <: Geometry

Zero-pressure-gradient flat plate geometry. Produces the classical Blasius boundary
layer for incompressible flow (β = 0).
"""
struct FlatPlate <: Geometry end

"""
    Wedge{T} <: Geometry

Wedge geometry for Falkner-Skan boundary layers.

# Fields
- `half_angle`: wedge half-angle in radians (angle between wall and freestream)

The Hartree pressure-gradient parameter is computed as β = 2 * half_angle / π,
so a half-angle of 0 gives the flat-plate (Blasius) case and π/2 gives stagnation flow.
"""
struct Wedge{T} <: Geometry
    half_angle::T
end

"""
    Cone{T} <: Geometry

Axisymmetric cone geometry. Solved via the Mangler transformation, which maps
the cone boundary layer to a flat-plate problem with a modified similarity variable.

# Fields
- `half_angle`: cone half-angle in radians
"""
struct Cone{T} <: Geometry
    half_angle::T
end

"""
    _β_parameter(geometry) -> T

Return the Hartree pressure-gradient parameter β for the given geometry.
For a `Wedge`, β = 2 * half_angle / π.  For `FlatPlate` and `Cone` (Mangler
transform reduces the cone to a flat-plate ODE), β = 0.
"""
_β_parameter(::FlatPlate) = 0.0
_β_parameter(g::Wedge)    = 2 * g.half_angle / π
_β_parameter(::Cone)      = 0.0   # solved as flat plate; Mangler applied in evaluate

"""
    IncompressibleBoundaryLayer{T} <: CanonicalFlow

Stores the self-similar solution to an incompressible laminar boundary layer.

# Fields
- `flow_prop`: `IncompressibleFlow{T}` with U_e and embedded `FluidProperties`
- `geometry`: `FlatPlate`, `Wedge{T}`, or `Cone{T}`
- `bc`: boundary condition at the wall
- `Re`: Reynolds number U_e L / ν_e used to solve this boundary layer
- `η`: similarity coordinate grid
- `f`, `fp`, `fpp`: Falkner-Skan function f and its first and second derivatives
"""
struct IncompressibleBoundaryLayer{T} <: CanonicalFlow
    flow_prop::IncompressibleFlow{T}
    geometry::Geometry
    bc::BoundaryCondition
    Re::T
    η::Vector{T}
    f::Vector{T}
    fp::Vector{T}
    fpp::Vector{T}
end

"""
    CompressibleBoundaryLayer{T} <: CanonicalFlow

Stores the self-similar solution to a compressible laminar boundary layer using the
Dorodnitsyn-Howarth transformation.

# Fields
- `flow_prop`: `CompressibleFlow{T}` with M∞, T∞, and embedded `FluidProperties`
- `geometry`: `FlatPlate`, `Wedge{T}`, or `Cone{T}`
- `bc`: boundary condition at the wall
- `Re`: Reynolds number ρ_e U_e L / μ_e used to solve this boundary layer
- `η`: Dorodnitsyn similarity coordinate grid
- `f`, `fp`, `fpp`: transformed stream function and derivatives
- `θ`, `θp`: temperature ratio T/T_e and its derivative w.r.t. η
- `U_e`: freestream velocity, derived from `flow_prop` and `fluid`
"""
struct CompressibleBoundaryLayer{T} <: CanonicalFlow
    flow_prop::CompressibleFlow{T}
    geometry::Geometry
    bc::BoundaryCondition
    Re::T
    η::Vector{T}
    f::Vector{T}
    fp::Vector{T}
    fpp::Vector{T}
    θ::Vector{T}
    θp::Vector{T}
    U_e::T
end
