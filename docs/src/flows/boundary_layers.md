# Boundary Layers

CanonicalFlows solves self-similar laminar boundary layers via the Falkner-Skan (incompressible)
and Levy-Lees/Dorodnitsyn-Howarth (compressible) formulations. Multiple dispatch on the
[`FlowProperties`](@ref), [`Geometry`](@ref), and [`BoundaryCondition`](@ref) types
selects the correct system of equations automatically.

## Setting up a flow

The first step is to define the flow conditions. Choose between [`IncompressibleFlow`](@ref)
and [`CompressibleFlow`](@ref) depending on whether compressibility effects are important.
Both types embed a [`FluidProperties`](@ref) object, which provides all transport and
thermodynamic constants — including those that make up the Prandtl number.

### Incompressible flow

[`IncompressibleFlow`](@ref) takes a freestream velocity `U_e` (m/s) and a fluid. The
Reynolds number is not part of the flow definition; it is supplied to the solver as a
separate argument so that the same flow condition can be evaluated at any length scale.

```julia
fluid = Air()

# Freestream velocity 100 m/s
flow = IncompressibleFlow(100.0, fluid)

# Default U_e = 1 (non-dimensional output)
flow = IncompressibleFlow(fluid)
```

### Compressible flow

[`CompressibleFlow`](@ref) takes the freestream Mach number, freestream (edge) temperature
in Kelvin, and a fluid. The freestream velocity is derived internally from M∞ and the
speed of sound.

```julia
fluid = Air()

# Mach 2, T∞ = 300 K
flow = CompressibleFlow(2.0, 300.0, fluid)

# Mach 5 hypersonic case
flow = CompressibleFlow(5.0, 300.0, fluid)
```

## Choosing a geometry

Three body shapes are available. All produce self-similar solutions; the cone uses the
Mangler transformation to reduce to a modified flat-plate ODE.

| Constructor | Description |
|:------------|:------------|
| `FlatPlate()` | Zero-pressure-gradient flat plate (Blasius/Dorodnitsyn). |
| `Wedge(half_angle)` | Symmetric wedge; `half_angle` in radians gives the Hartree parameter β = 2 half_angle / π. |
| `Cone(half_angle)` | Axisymmetric cone; `half_angle` in radians. Mangler transformation applied during evaluation. |

```julia
plate  = FlatPlate()
wedge  = Wedge(π/4)   # β = 0.5, favourable pressure gradient
cone   = Cone(π/6)    # 30° half-angle cone
```

## Choosing a wall boundary condition

Wall thermal boundary conditions control the energy equation (compressible) or simply
mark the wall as insulating for heat-transfer diagnostics (incompressible).

```julia
# Zero heat flux at the wall
bc = AdiabaticBC()

# Fixed wall temperature (K)
bc = IsothermalBC(500.0)    # cooled wall relative to adiabatic temperature
bc = IsothermalBC(1500.0)   # heated wall
```

## Solving the boundary layer

Pass the flow, geometry, Reynolds number, and boundary condition to [`BoundaryLayer`](@ref).
The Reynolds number `Re` (U_e L / ν_e for incompressible; ρ_e U_e L / μ_e for compressible)
is a solve-time parameter that sets the physical length scale without altering the flow
or fluid description. The boundary conditions are supplied as a tuple to allow future
support for multiple BCs.

```julia
using CanonicalFlows

fluid = Air()
Re    = 1e6

# --- Incompressible examples ---

# Blasius flat plate
bl_blasius = BoundaryLayer(IncompressibleFlow(1.0, fluid), FlatPlate(), Re, (AdiabaticBC(),))

# Falkner-Skan wedge (β = 0.5)
bl_wedge = BoundaryLayer(IncompressibleFlow(1.0, fluid), Wedge(π/4), Re, (AdiabaticBC(),))

# Axisymmetric cone
bl_cone = BoundaryLayer(IncompressibleFlow(1.0, fluid), Cone(π/4), Re, (AdiabaticBC(),))

# --- Compressible examples ---

# Mach 5 adiabatic flat plate
bl_adiabatic = BoundaryLayer(CompressibleFlow(5.0, 300.0, fluid), FlatPlate(), Re, (AdiabaticBC(),))

# Mach 2 isothermal wall
bl_isothermal = BoundaryLayer(CompressibleFlow(2.0, 300.0, fluid), FlatPlate(), Re, (IsothermalBC(400.0),))

# Mach 3 wedge
bl_comp_wedge = BoundaryLayer(CompressibleFlow(3.0, 300.0, fluid), Wedge(π/6), Re, (AdiabaticBC(),))
```

The solver returns an [`IncompressibleBoundaryLayer`](@ref) or [`CompressibleBoundaryLayer`](@ref)
object depending on the flow type. Both store the similarity-coordinate grid `η` and the
solved fields directly as vectors:

```julia
bl = BoundaryLayer(IncompressibleFlow(1.0, Air()), FlatPlate(), 1e6, (AdiabaticBC(),))

bl.η    # similarity coordinate grid
bl.f    # stream function f(η)
bl.fp   # f′(η) = u / U_e
bl.fpp  # f″(η) — proportional to wall shear when evaluated at η = 0
bl.Re   # Reynolds number used for this solution
```

For compressible flows, the temperature ratio and its derivative are also available:

```julia
bl = BoundaryLayer(CompressibleFlow(5.0, 300.0, Air()), FlatPlate(), 1e6, (AdiabaticBC(),))

bl.θ    # temperature ratio T / T_e (η)
bl.θp   # dθ/dη

# Adiabatic wall temperature
T_aw = bl.flow_prop.T∞ * bl.θ[1]
```

## Evaluating profiles in physical space

Use [`evaluate_profile`](@ref) to obtain physical-space profiles at a given streamwise
location `x` (non-dimensionalised by the reference length L) over a range of
wall-normal positions `y`.

```julia
bl = BoundaryLayer(CompressibleFlow(5.0, 300.0, Air()), FlatPlate(), 1e6, (AdiabaticBC(),))

x       = 0.5                        # x/L = 0.5
y_range = LinRange(0.0, 0.01, 200)   # y/L grid

# Request multiple variables in one call
u, v, T, ρ = evaluate_profile(bl, x, y_range,
                               (:streamwise_vel, :wallnormal_vel, :temperature, :density))
```

Available variable symbols:

| Symbol | Description | Incompressible | Compressible |
|:-------|:------------|:--------------:|:------------:|
| `:streamwise_vel` | Streamwise velocity u (m/s) | yes | yes |
| `:wallnormal_vel` | Wall-normal velocity v (m/s) | yes | yes |
| `:pressure` | Static pressure (Pa); 0 (gauge) for incompressible | — | yes |
| `:temperature` | Static temperature (K) | — | yes |
| `:density` | Density (kg/m³) | — | yes |

### Finding the boundary layer edge

[`boundary_layer_edge`](@ref) returns the similarity coordinate η at which the
velocity profile reaches the freestream within the specified tolerance.

```julia
bl = BoundaryLayer(IncompressibleFlow(1.0, Air()), FlatPlate(), 1e6, (AdiabaticBC(),))

η_edge = boundary_layer_edge(bl)            # default tol = 1e-3 (0.1% deficit)
η_edge = boundary_layer_edge(bl; tol=0.01)  # classical 99% thickness definition
```

## API Reference

### Flow properties

```@docs
FlowProperties
IncompressibleFlow
CompressibleFlow
```

### Geometry

```@docs
Geometry
FlatPlate
Wedge
Cone
```

### Boundary conditions

```@docs
BoundaryCondition
AdiabaticBC
IsothermalBC
```

### Solution types

```@docs
CanonicalFlow
IncompressibleBoundaryLayer
CompressibleBoundaryLayer
```

### Solver

```@docs
BoundaryLayer
```

### Evaluation

```@docs
evaluate_profile
boundary_layer_edge
```
