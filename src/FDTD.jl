module FDTD
    using LinearAlgebra
    using Strided
    using DiscreteAxis

    const EHTuple{T,N} = Tuple{Tuple{Array{T,N},Array{T,N},Array{T,N}},Tuple{Array{T,N},Array{T,N},Array{T,N}}}
    const VecArray{T,N} = Tuple{Array{T,N},Array{T,N},Array{T,N}}

    include("pde_utils.jl")
    include("yee.jl")
    include("fdtd_utils.jl")

    include("update_equations.jl")
    include("cpml.jl")
    include("propagation.jl")

    export FDTD_propagate, radar_propagate, setup_spacetime
    export EMField, VecField3, Medium, Coefficients, PML, pontying, power_density, vec_array
    export PMLAux
    export ε₀, μ₀, η₀, c₀
end
