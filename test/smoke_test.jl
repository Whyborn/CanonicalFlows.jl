using CanonicalFlows

fluid = Air()

println("=== Incompressible flat plate ===")
flow_inc = IncompressibleFlow(1.0, fluid)
bl_inc = BoundaryLayer(flow_inc, FlatPlate(), 1e6, (AdiabaticBC(),))
println("f''(0) ≈ ", bl_inc.fpp[1], "  (Blasius: 0.3321)")
println("f'(0)  = ", bl_inc.fp[1],  "  (expected 0)")
println("f'(∞)  ≈ ", bl_inc.fp[end], "  (expected 1)")

println("\n=== Incompressible wedge, half_angle=π/4 → β=0.5 ===")
bl_wedge = BoundaryLayer(flow_inc, Wedge(π/4), 1e6, (AdiabaticBC(),))
println("f''(0) ≈ ", bl_wedge.fpp[1])

println("\n=== Incompressible cone α=45° ===")
bl_cone = BoundaryLayer(flow_inc, Cone(π/4), 1e6, (AdiabaticBC(),))
println("f''(0) ≈ ", bl_cone.fpp[1], "  (same as flat plate, Mangler applies in evaluate_profile)")

println("\n=== Compressible flat plate, adiabatic ===")
flow_comp = CompressibleFlow(5.0, 300.0, fluid)
bl_comp = BoundaryLayer(flow_comp, FlatPlate(), 1e6, (AdiabaticBC(),))
println("f'(0)  = ", bl_comp.fp[1],  "  (expected 0)")
println("f'(∞)  ≈ ", bl_comp.fp[end], "  (expected 1)")
println("θ(∞)   ≈ ", bl_comp.θ[end], "  (expected 1)")
T_aw = bl_comp.flow_prop.T∞ * bl_comp.θ[1]
println("T_aw   ≈ ", round(T_aw, digits=1), " K  (adiabatic wall temperature)")

println("\n=== evaluate_profile (compressible) ===")
p, T, u, v = evaluate_profile(bl_comp, 0.2, LinRange(0.0, 0.02, 6),
                               (:pressure, :temperature, :streamwise_vel, :wallnormal_vel))
println("p (Pa) = ", round.(p, sigdigits=4))
println("T (K)  = ", round.(T, digits=1))
println("u (m/s)= ", round.(u, digits=2))
println("v (m/s)= ", round.(v, sigdigits=4))

println("\n=== prandtl_number derived from fluid ===")
println("Pr = ", round(prandtl_number(fluid), digits=4), "  (expected ≈ 0.71)")
