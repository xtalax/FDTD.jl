module FDTD
    using LinearAlgebra, StaticArrays, Strided, Rotations, SparseArrays
    using DiscreteAxis, DiffEqOperators, GeometryTypes
    using ProgressMeter
    const EHTuple{T,N} = Tuple{Tuple{Array{T,N},Array{T,N},Array{T,N}},Tuple{Array{T,N},Array{T,N},Array{T,N}}}
    const VecArray{T,N} = Tuple{Array{T,N},Array{T,N},Array{T,N}}

    include("pde_utils.jl")
    include("yee.jl")
    include("fdtd_utils.jl")

    include("update_equations.jl")
    include("cpml.jl")
    include("propagation.jl")

    export FDTD_propagate, radar_propagate, setup_spacetime2D, setup_spacetime3D
    export EMField, VecField3, Medium, Coefficients, PML, pontying, power_density, VecArray
    export PMLAux
    export ε₀, μ₀, η₀, c₀ #Export Constants
end
