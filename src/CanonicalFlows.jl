module CanonicalFlows

using LinearAlgebra: lu, ldiv!

include("types.jl")
include("fluid_properties.jl")
include("boundary_conditions.jl")
include("similarity_odes.jl")
include("matrix_solver.jl")
include("problems.jl")
include("evaluate.jl")

export CanonicalFlow
export FlowProperties, IncompressibleFlow, CompressibleFlow
export Geometry, FlatPlate, Wedge, Cone
export FluidProperties, Air, viscosity, prandtl_number
export BoundaryCondition, AdiabaticBC, IsothermalBC
export IncompressibleBoundaryLayer, CompressibleBoundaryLayer
export BoundaryLayer
export evaluate_profile, boundary_layer_edge

end
