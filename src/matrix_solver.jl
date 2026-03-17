
# Matrix solver using a Newton-finite-difference method.
# The nonlinear BVP is discretised on a uniform η grid using a first-order
# forward-difference scheme. At each Newton step the resulting linear system
# J Δu = -F is solved with LinearAlgebra.lu.

# ---------------------------------------------------------------------------
# Boundary condition dispatch helpers
# ---------------------------------------------------------------------------

# Wall thermal BC residual for the compressible residual vector
_θ_bc_residual(u, ::AdiabaticBC, ::Float64)  = u[5]               # θ'(0) = 0
_θ_bc_residual(u, bc::IsothermalBC, T_e) = u[4] - bc.T_w / T_e   # θ(0) = T_w/T_e

# Wall temperature ratio estimate for the initial guess
_θ_wall_guess_matrix(::AdiabaticBC, flow::CompressibleFlow) =
    1.0 + sqrt(prandtl_number(flow.fluid)) * (flow.fluid.γ - 1) / 2 * flow.M∞^2
_θ_wall_guess_matrix(bc::IsothermalBC, flow::CompressibleFlow) =
    bc.T_w / flow.T∞

# ---------------------------------------------------------------------------
# Residual assembly
# ---------------------------------------------------------------------------

"""
    _residual_incompressible!(F, u, η_grid, β)

Fill the residual vector `F` for the incompressible Falkner-Skan system.

Layout: `u = [f₁, f′₁, f″₁, f₂, f′₂, f″₂, …]` (local ordering, 1-based).
Grid has `N` points; `length(u) = 3N`.
"""
function _residual_incompressible!(F, u, η_grid, β)
    N  = length(η_grid)
    Δη = η_grid[2] - η_grid[1]

    idx(i, j) = 3(i - 1) + j

    F[1] = u[idx(1, 1)]           # f(0) = 0
    F[2] = u[idx(1, 2)]           # f′(0) = 0
    F[3] = u[idx(N, 2)] - 1.0    # f′(η_max) = 1

    for i in 1:(N - 1)
        row = 3 + 3(i - 1)
        fi  = u[idx(i, 1)];  fip = u[idx(i + 1, 1)]
        gi  = u[idx(i, 2)];  gip = u[idx(i + 1, 2)]
        hi  = u[idx(i, 3)];  hip = u[idx(i + 1, 3)]

        F[row + 1] = (fip - fi) / Δη - gi
        F[row + 2] = (gip - gi) / Δη - hi
        F[row + 3] = (hip - hi) / Δη +
                     (β + 1) / 2 * fi * hi +
                     β * (1 - gi^2)
    end
end

"""
    _residual_compressible!(F, u, η_grid, β, Pr, γ, M_e, T_e, fluid, bc)

Fill the residual vector for the compressible BL system.

Layout: `u = [f₁, f′₁, f″₁, θ₁, θ′₁, f₂, …]` (local ordering, 1-based).
Grid has `N` points; `length(u) = 5N`.
"""
function _residual_compressible!(F, u, η_grid, β, Pr, γ, M_e, T_e, fluid, bc)
    N  = length(η_grid)
    Δη = η_grid[2] - η_grid[1]

    idx(i, j) = 5(i - 1) + j

    F[1] = u[idx(1, 1)]           # f(0) = 0
    F[2] = u[idx(1, 2)]           # f′(0) = 0
    F[3] = _θ_bc_residual(u, bc, T_e)
    F[4] = u[idx(N, 2)] - 1.0    # f′(η_max) = 1
    F[5] = u[idx(N, 4)] - 1.0    # θ(η_max) = 1

    for i in 1:(N - 1)
        row = 5 + 5(i - 1)
        fi   = u[idx(i, 1)];  fip  = u[idx(i + 1, 1)]
        gi   = u[idx(i, 2)];  gip  = u[idx(i + 1, 2)]
        hi   = u[idx(i, 3)];  hip  = u[idx(i + 1, 3)]
        θi   = u[idx(i, 4)];  θip  = u[idx(i + 1, 4)]
        θpi  = u[idx(i, 5)];  θpip = u[idx(i + 1, 5)]

        Ci    = chapman_rubesin(fluid, T_e, θi)
        dCdθi = chapman_rubesin_dθ(fluid, T_e, θi)
        dCdηi = dCdθi * θpi

        F[row + 1] = (fip - fi) / Δη - gi
        F[row + 2] = (gip - gi) / Δη - hi
        F[row + 3] = (hip - hi) / Δη +
                     (β + 1) / 2 * fi * hi / Ci +
                     β * (1 - gi^2 / Ci) +
                     dCdηi * hi / Ci
        F[row + 4] = (θip - θi) / Δη - θpi
        F[row + 5] = (θpip - θpi) / Δη +
                     dCdηi / Ci * θpi +
                     Pr * (β + 1) / 2 * fi * θpi / Ci +
                     Pr * (γ - 1) * M_e^2 * hi^2
    end
end

# ---------------------------------------------------------------------------
# Numerical Jacobian (column-wise finite differences)
# ---------------------------------------------------------------------------

function _numerical_jacobian!(J, F!, u, F0, params)
    n  = length(u)
    ε  = 1e-7
    F1 = similar(F0)
    for j in 1:n
        u[j] += ε
        F!(F1, u, params...)
        @. J[:, j] = (F1 - F0) / ε
        u[j] -= ε
    end
end

# ---------------------------------------------------------------------------
# Newton solver
# ---------------------------------------------------------------------------

function _newton_solve!(F!, u, params; tol = 1e-10, max_iter = 100)
    n  = length(u)
    F0 = zeros(n)
    J  = zeros(n, n)

    for iter in 1:max_iter
        F!(F0, u, params...)
        res = maximum(abs, F0)
        res < tol && return iter

        _numerical_jacobian!(J, F!, u, F0, params)
        fact = lu(J)
        ldiv!(fact, F0)
        u .-= F0
    end
    @warn "Newton solver did not converge within $max_iter iterations."
end

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    solve_matrix(flow::IncompressibleFlow, geometry, bc; η_max, N_pts)

Solve the incompressible Falkner-Skan BVP with a Newton-FD matrix method.
Returns `(η, f, fp, fpp)`.
"""
function solve_matrix(flow::IncompressibleFlow, geometry::Geometry,
                      bc::BoundaryCondition;
                      η_max = 15.0, N_pts = 200)
    β      = _β_parameter(geometry)
    η_grid = range(0.0, η_max, N_pts) |> collect

    u = zeros(3N_pts)
    η_c = 3.0
    for i in 1:N_pts
        ηi = η_grid[i]
        if ηi < η_c
            u[3(i-1)+2] = ηi / η_c
            u[3(i-1)+1] = ηi^2 / (2η_c)
            u[3(i-1)+3] = 1.0 / η_c
        else
            u[3(i-1)+2] = 1.0
            u[3(i-1)+1] = ηi - η_c / 2
            u[3(i-1)+3] = 0.0
        end
    end

    _newton_solve!(_residual_incompressible!, u, (η_grid, β))

    f   = [u[3(i-1)+1] for i in 1:N_pts]
    fp  = [u[3(i-1)+2] for i in 1:N_pts]
    fpp = [u[3(i-1)+3] for i in 1:N_pts]
    return η_grid, f, fp, fpp
end

"""
    solve_matrix(flow::CompressibleFlow, geometry, bc; η_max, N_pts)

Solve the compressible BL BVP with a Newton-FD matrix method.
Returns `(η, f, fp, fpp, θ, θp)`.
"""
function solve_matrix(flow::CompressibleFlow, geometry::Geometry,
                      bc::BoundaryCondition;
                      η_max = 20.0, N_pts = 300)
    fluid  = flow.fluid
    β      = _β_parameter(geometry)
    Pr     = prandtl_number(fluid)
    γ      = fluid.γ
    M_e    = flow.M∞
    T_e    = flow.T∞
    η_grid = range(0.0, η_max, N_pts) |> collect

    θ_w_est = _θ_wall_guess_matrix(bc, flow)

    u = zeros(5N_pts)
    η_c = 3.0
    for i in 1:N_pts
        ηi = η_grid[i]
        if ηi < η_c
            u[5(i-1)+2] = ηi / η_c
            u[5(i-1)+1] = ηi^2 / (2η_c)
            u[5(i-1)+3] = 1.0 / η_c
        else
            u[5(i-1)+2] = 1.0
            u[5(i-1)+1] = ηi - η_c / 2
            u[5(i-1)+3] = 0.0
        end
        decay = exp(-ηi / 3.0)
        u[5(i-1)+4] = 1.0 + (θ_w_est - 1.0) * decay
        u[5(i-1)+5] = -(θ_w_est - 1.0) / 3.0 * decay
    end

    _newton_solve!(_residual_compressible!, u, (η_grid, β, Pr, γ, M_e, T_e, fluid, bc))

    f   = [u[5(i-1)+1] for i in 1:N_pts]
    fp  = [u[5(i-1)+2] for i in 1:N_pts]
    fpp = [u[5(i-1)+3] for i in 1:N_pts]
    θ   = [u[5(i-1)+4] for i in 1:N_pts]
    θp  = [u[5(i-1)+5] for i in 1:N_pts]
    return η_grid, f, fp, fpp, θ, θp
end
