"""
Step time for a TEM grid (2D)
"""
function step_time(F::EHTuple{T,2}, Φ::PML, C::Coefficients, TEsource, TMsource, sourceindex) where T # F is the EM field, C is a wrapper for the Coefficients used here
    E, H = F

    H = update_H(E, H, C)
    E = update_E(E, H, C)

    ΦEaux = updatePML(Φ.E, H, C, 1, Φ.N)
    ΦHaux = updatePML(Φ.H, E, C, 2, Φ.N)


    H = applyPML(H, ΦHaux, C, 1, Φ.N)
    E = applyPML(E, ΦEaux, C, 2, Φ.N)

    Φ = PML(ΦEaux,ΦHaux,Φ.N)

    E[3][sourceindex...] = E[3][sourceindex...] .+ TEsource #hard source acts like antenna

    H[3][sourceindex...] = H[3][sourceindex...] .+ TMsource #hard source acts like antenna

    return ((E, H), Φ)
end

"""
Step time for a full 3D grid
"""
function step_time(F::EHTuple{T}, Φ::PML, C::Coefficients, Esourcenow, Hsourcenow, sourceindex) where T # F is the EM field, C is a wrapper for the Coefficients used here
    E, H = F

    H = update_H(E, H, C)
    E = update_E(E, H, C)

    ΦEaux = updatePML(Φ.E, H, C, 1, Φ.N)
    ΦHaux = updatePML(Φ.H, E, C, 2, Φ.N)

    H = applyPML(H, ΦHaux, C, 1, Φ.N)
    E = applyPML(E, ΦEaux, C, 2, Φ.N)

    Φ = PML(ΦEaux,ΦHaux,Φ.N)
    E[sourceindex[1]][sourceindex[2]...] = E[sourceindex[1]][sourceindex[2]...] .+ Esourcenow #hard source acts like antenna
    Hnew[sourceindex[1]+ (sourceindex[1]==3 ? -2 : 1)][sourceindex[2]...] = F[2][sourceindex[1]+ (sourceindex[1]==3 ? -2 : 1)][sourceindex[2]...] + Hsourcenow
    return ((E, H), Φ)
end


function FDTD_propagate(space::ProductAxis, time::LinearAxis, f₀::T, nPML::Int = 0;
                        source,
                        medium::Medium,
                        sourceindex = origin(space),
                        detectorindex =  (:, [:, :, :]),
                        reduction = identity
                        ) where T


      Ex = zeros(T, size(space)...)
      Ey = zeros(T, size(space)...)
      Ez = zeros(T, size(space)...)

      Hx = zeros(T, size(space)...)
      Hy = zeros(T, size(space)...)
      Hz = zeros(T, size(space)...)

      F = (E,H) = ((Ex,Ey,Ez), (Hx,Hy,Hz))

   Φ = PML(space, nPML)
   C = Coefficients(space, time, medium, f₀, nPML)

    if detectorindex isa Tuple{Colon, Array{Colon,1}}
        out = Array{EMField{T}}(undef,length(time))
        #out = Array{PML{T}}(undef,length(time))
        EM_out = true
    else
        Sout = size(F[detectorindex[1]][detectorindex[2]...])
        Nout = length(Sout)
        out = zeros(Sout..., length(time))
        colons = fill(:, Nout)
        EM_out = false
    end
    @showprogress "Simulating with FDTD..." for i in time.i
        Fnew, Φnew = step_time(F, Φ, C, source[i], source[i], sourceindex)
        if EM_out
             out[i] =EMField(Fnew)
        else
             out[colons...,i] = Fnew[detectorindex[1]][detectorindex[2]][detectorindex[3]...] # Don't mess with this, it works but I can't remember how
        end
        F = Fnew;
        Φ = Φnew;
    end
    return out
end


function detect(F::EHTuple{T}, detectorindex::CartesianIndices) where T
    out = zero(T)
    for I in detectorindex
        out = out + F[1][I]^2
        out = out + (F[2][I])^2
        out = out + F[3][I]^2
    end
    return out/(2*η₀ * prod(size(detectorindex)))
end

"""
A draft function for generating the radar cross section of a particular medium
"""
function radar_propagate(space::ProductAxis, time::LinearAxis, f₀::T, nPML::Int = 0;
                        source,
                        medium::Medium,
                        sourceindex = origin(space),
                        detectorindex::CartesianIndices = CartesianIndices(sourceindex[2])
                        ) where T
   Ex = zeros(T, size(space)...)
   Ey = zeros(T, size(space)...)
   Ez = zeros(T, size(space)...)

   Hx = zeros(T, size(space)...)
   Hy = zeros(T, size(space)...)
   Hz = zeros(T, size(space)...)

   F = (E,H) = ((Ex,Ey,Ez), (Hx,Hy,Hz))
   Φ = PML(space, nPML)
   C = Coefficients(space, time, medium, f₀, nPML)
   out = Vector{T}(undef, time.N)
   for i in 1:time.N
        Fnew, Φnew = step_time(F, Φ, C, source[i], source[i], sourceindex) #You actually need Fnew here to get the compiler to not tidy the loop away
        @inbounds out[i] = sum([sum(view(Fnew[1][k], detectorindex).^2)./(2*η₀ * prod(size(detectorindex))) for k in 1:3])
        F = Fnew ; Φ = Φnew
    end
    return out
end
