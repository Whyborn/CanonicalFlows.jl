
"""
    AdiabaticBC <: BoundaryCondition

Adiabatic (zero heat-flux) wall boundary condition: dT/dy|_w = 0.
"""
struct AdiabaticBC <: BoundaryCondition end

"""
    IsothermalBC{T} <: BoundaryCondition

Isothermal wall boundary condition with specified wall temperature.

# Fields
- `T_w`: wall temperature in the same units as the freestream temperature
"""
struct IsothermalBC{T} <: BoundaryCondition
    T_w::T
end
