# Fluid Properties

The [`FluidProperties`](@ref) struct holds the thermodynamic and transport constants for a
perfect gas. Rather than storing the Prandtl number directly, it stores the three physical
coefficients that define it — specific heat `cp`, dynamic viscosity `μ_ref`, and thermal
conductivity `k` — so that Pr = cp · μ / k is always traceable to first principles.
Use [`prandtl_number`](@ref) to compute Pr from a `FluidProperties` instance.

These properties feed directly into the Chapman-Rubesin viscosity correction and the
energy equation used in the compressible boundary layer solver.

## Built-in fluids

### Air

[`Air()`](@ref) returns standard dry-air properties:

| Property | Symbol | Value |
|:---------|:-------|:------|
| Ratio of specific heats | γ | 1.4 |
| Specific gas constant | R | 287.05 J/(kg·K) |
| Specific heat at constant pressure | cp | 1005.0 J/(kg·K) |
| Reference dynamic viscosity | μ_ref | 1.716 × 10⁻⁵ Pa·s |
| Reference temperature | T_ref | 273.15 K |
| Sutherland constant | S | 110.4 K |
| Reference thermal conductivity | k | 0.02428 W/(m·K) |

```julia
fluid = Air()

# Derived Prandtl number
Pr = prandtl_number(fluid)   # ≈ 0.71
```

## Computing viscosity

[`viscosity`](@ref) applies Sutherland's law to compute dynamic viscosity at any
temperature:

```julia
fluid = Air()

μ_300 = viscosity(fluid, 300.0)   # Pa·s at 300 K
μ_600 = viscosity(fluid, 600.0)   # Pa·s at 600 K — roughly double
```

## Custom fluids

A [`FluidProperties`](@ref) object can be constructed directly to define a custom fluid.
The Prandtl number is implicitly set by the choice of `cp`, `μ_ref`, and `k`.

```julia
# Custom gas: Pr = cp * μ_ref / k = 1100 * 2.0e-5 / 0.031 ≈ 0.71
my_gas = FluidProperties(
    1.35,      # γ
    260.0,     # R  (J/(kg·K))
    1100.0,    # cp (J/(kg·K))
    2.0e-5,    # μ_ref (Pa·s)
    280.0,     # T_ref (K)
    120.0,     # S (K)
    0.031,     # k (W/(m·K))
)

Pr = prandtl_number(my_gas)
μ  = viscosity(my_gas, 500.0)
```

## API Reference

```@docs
FluidProperties
Air
prandtl_number
viscosity
```
