
# ODE right-hand-side functions for the similarity equations.
# All use the Dorodnitsyn-Howarth normalisation:
#   η = y * √(U_e / (ν_e * x))   (flat plate; modified for wedge/cone)
#
# Incompressible Falkner-Skan:
#   f''' + (β+1)/2 * f * f'' + β * (1 - f'^2) = 0
#   β = 0 : flat plate (Blasius)
#   β = 1 : stagnation point (Hiemenz)
#
# Compressible (Levy-Lees / Dorodnitsyn):
#   (C f'')' + (β+1)/2 * f * f'' + β * (C - f'^2) = 0     [momentum]
#   (C/Pr * θ')' + (β+1)/2 * f * θ' + (γ-1) M² C (f'')² = 0 [energy]
#   where C = ρμ/(ρ_e μ_e) = chapman_rubesin(fluid, T_e, θ)

"""
    falkner_skan_ode!(du, u, p, η)

In-place ODE RHS for the incompressible Falkner-Skan equation.
State: `u = [f, f', f'']`.  Parameter: `p = [β]`.
"""
function falkner_skan_ode!(du, u, p, η)
    β = p[1]
    f, fp, fpp = u[1], u[2], u[3]
    du[1] = fp
    du[2] = fpp
    du[3] = -(β + 1) / 2 * f * fpp - β * (1 - fp^2)
end

"""
    compressible_bl_ode!(du, u, p, η)

In-place ODE RHS for the compressible boundary layer similarity equations.
State: `u = [f, f', f'', θ, θ']` where θ = T / T_e.
Parameters: `p = (β, Pr, γ, M_e, fluid, T_e)`.
"""
function compressible_bl_ode!(du, u, p, η)
    β, Pr, γ, M_e, fluid, T_e = p
    f, fp, fpp, θ, θp = u[1], u[2], u[3], u[4], u[5]

    C    = chapman_rubesin(fluid, T_e, θ)
    dCdθ = chapman_rubesin_dθ(fluid, T_e, θ)
    dCdη = dCdθ * θp      # chain rule

    # Momentum: C f''' + C' f'' + (β+1)/2 * f * f'' + β*(C - f'^2) = 0
    fppp = (-(β + 1) / 2 * f * fpp - β * (C - fp^2) - dCdη * fpp) / C

    # Energy: C/Pr * θ'' + (C/Pr)' * θ' + (β+1)/2 * f * θ' + (γ-1)*M²*C*f''^2 = 0
    # (C/Pr)' = dCdη/Pr
    θpp = -(dCdη / C * θp + Pr * (β + 1) / 2 * f * θp / C + Pr * (γ - 1) * M_e^2 * fpp^2)

    du[1] = fp
    du[2] = fpp
    du[3] = fppp
    du[4] = θp
    du[5] = θpp
end
