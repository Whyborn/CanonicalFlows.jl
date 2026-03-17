# CanonicalFlows.jl

**CanonicalFlows.jl** is a Julia package for computing self-similar solutions to laminar
boundary layers. It solves the coupled ODE systems that govern incompressible and
compressible boundary layers over flat plates, wedges, and cones, and provides
utilities for evaluating the resulting flow profiles in physical space.

## Overview

A self-similar boundary layer solution collapses the two-dimensional velocity (and
temperature) field into a single ODE in the similarity coordinate η. CanonicalFlows
solves these ODEs using a second-order finite-difference Newton iteration and provides
a clean interface for constructing and querying the solutions.

The package is built around three concepts:

| Concept | Description |
|:--------|:------------|
| **Flow properties** | Dimensionless parameters that characterise the flow: Reynolds number, Mach number, Prandtl number, etc. |
| **Geometry** | The body shape — flat plate, wedge, or axisymmetric cone. |
| **Boundary condition** | The thermal condition at the wall — adiabatic or isothermal. |

A boundary layer is solved by combining these three ingredients with a fluid model:

```julia
using CanonicalFlows

# Mach 5 flat plate, adiabatic wall, Re = 10^6
flow = CompressibleFlow(5.0, 300.0, Air())
bl   = BoundaryLayer(flow, FlatPlate(), 1e6, (AdiabaticBC(),))
```

The returned object stores the similarity-coordinate grid and all solved fields. It
can then be evaluated at any physical location with [`evaluate_profile`](@ref).

## Contents

```@contents
Pages = ["flows/boundary_layers.md", "fluid_properties.md"]
Depth = 2
```
