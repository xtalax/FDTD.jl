function setup_spacetime2D(fmax::Number, multiplier::Number, npadding::Int = 0; xlims::NTuple{2,T}, ylims::NTuple{2,T}) where T
    pts_per_wavelength = 16
    factor = c₀*1.2/Sc
    Δ = c₀ /(pts_per_wavelength*fmax)
    Δt = Δ/factor
    x = LinearAxis(xlims[1]-Δ*(npadding+1),xlims[2]+Δ*(npadding+1),Δ) # compensating axis to the new units
    y = LinearAxis(ylims[1]-Δ*(npadding+1),ylims[2]+Δ*(npadding+1),Δ)

    space = Space2D(x,y)

    Nt = round(Int, multiplier*(maximum(domain.(space)))/(c₀*Δt))

    tmin = 0.0 ; tmax = tmin + Δt*Nt
    time = LinearAxis(tmin, tmax, Δt)
    for axis in space
        @assert axis.N > 2*npadding
    end

    return (space, time)
end

function setup_spacetime3D(fmax::Number, multiplier::Number, npadding::Int = 0; xlims::NTuple{2,T}, ylims::NTuple{2,T}, zlims::NTuple{2,T}) where T
    pts_per_wavelength = 16
    factor = c₀*1.2/Sc
    Δ = c₀ /(pts_per_wavelength*fmax)
    Δt = Δ/factor
    x = LinearAxis(xlims[1]-Δ*(npadding+1),xlims[2]+Δ*(npadding+1),Δ) # compensating axis to the new units
    y = LinearAxis(ylims[1]-Δ*(npadding+1),ylims[2]+Δ*(npadding+1),Δ)
    space = Space2D(x,y,z)

    Nt = round(Int, multiplier*(maximum(domain.(space)))/(c₀*Δt))

    tmin = 0.0 ; tmax = tmin + Δt*Nt
    time = LinearAxis(tmin, tmax, Δt)
    for axis in space
        @assert axis.N > 2*npadding
    end

    return (space, time)
end


curl(F, δ, dims) = ⊗(δ[dims[2]], F[dims[1]], dims[2]) .-
                        ⊗(δ[dims[1]], F[dims[2]],dims[1])

@inline function cover3d(A::AbstractVector, dim::Integer)
    extension = ones(Int64, 3)
    if dim == 1
        perm = (1, 2, 3)
    elseif dim == 2
        perm = (2, 1, 3)
    elseif dim == 3
        perm = (2, 3, 1)
    else
        throw("dim must be less than 3")
    end
    return permutedims(repeat(A, 1, 1, 1), perm)
end

@inline function cover2d(A::AbstractVector, dim::Integer)
    extension = ones(Int64, 2)
    if dim == 1
        perm = (1, 2)
    elseif dim == 2
        perm = (2, 1)
    else
        throw("dim must be less than 3")
    end
    return permutedims(repeat(A, 1, 1), perm)
end

perpsize(A::AbstractArray, i::Int) = size(A)[setdiff(1:ndims(A), i)]
