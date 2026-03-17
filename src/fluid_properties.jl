
"""
    Air()

Returns `FluidProperties` for dry air at standard conditions using standard values.
"""
function Air()
    FluidProperties(
        1.4,        # γ
        287.05,     # R (J/(kg·K))
        1005.0,     # cp (J/(kg·K))
        1.716e-5,   # μ_ref at T_ref = 273.15 K (Pa·s)
        273.15,     # T_ref (K)
        110.4,      # S (K)
        0.02428,    # k at T_ref (W/(m·K))
    )
end

"""
    prandtl_number(fluid::FluidProperties) -> Float64

Compute the (reference) Prandtl number from fluid transport properties:
Pr = cp · μ_ref / k.
"""
function prandtl_number(fluid::FluidProperties)
    return fluid.cp * fluid.μ_ref / fluid.k
end

"""
    viscosity(fluid::FluidProperties, T)

Compute dynamic viscosity (Pa·s) at temperature `T` (K) using Sutherland's law.
"""
function viscosity(fluid::FluidProperties, T::Real)
    return fluid.μ_ref * (T / fluid.T_ref)^1.5 * (fluid.T_ref + fluid.S) / (T + fluid.S)
end

"""
    chapman_rubesin(fluid, T_e, θ)

Compute C(θ) = ρμ/(ρ_e μ_e) for a perfect gas with Sutherland viscosity, where θ = T/T_e.

For constant-pressure flow: ρ/ρ_e = T_e/T = 1/θ and μ/μ_e given by Sutherland's law,
giving C(θ) = θ^(1/2) * (T_e + S) / (T_e*θ + S).
"""
function chapman_rubesin(fluid::FluidProperties, T_e::Real, θ::Real)
    S = fluid.S
    return sqrt(θ) * (T_e + S) / (T_e * θ + S)
end

"""
    chapman_rubesin_dθ(fluid, T_e, θ)

Compute dC/dθ analytically for the Chapman-Rubesin parameter C(θ).
"""
function chapman_rubesin_dθ(fluid::FluidProperties, T_e::Real, θ::Real)
    S = fluid.S
    A = T_e + S
    denom = T_e * θ + S
    # C = A * θ^(1/2) / denom
    # dC/dθ = A * [0.5 * θ^(-1/2) * denom - θ^(1/2) * T_e] / denom^2
    #        = A * θ^(-1/2) * [0.5*(T_e*θ + S) - T_e*θ] / denom^2
    #        = A * θ^(-1/2) * 0.5*(S - T_e*θ) / denom^2
    return A * 0.5 * (S - T_e * θ) / (sqrt(θ) * denom^2)
end
