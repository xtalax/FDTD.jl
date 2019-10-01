#Plane wave expansion

function plane_wave_expansion(F::EHTuple{T,N}, d::Integer) where T
    F_k = (zeros(Complex{T}, size(F)), zeros(Complex{T}, size(F)))
    for (i,field) in F
        F_k[i] .= fft(field[d])./length(F)
    end
    return F_k
end
function inv_plane_wave_expansion(F_k::EHTuple{Complex{T},N}) where T
    F_k = (zeros(T, size(F)), zeros(T, size(F)))
    for (i,field) in F_k
        F[i][j] .= real.(ifft(field).*length(F_k))
    return F
end

function square(s::AbstractVector{L}, I::CartesianIndex) where {L<:LinearAxis}
    return @inbounds [s[I[i]]^2 for i in 1:length[I]]]
end

function compact_range_transform(F::EHTuple{T1,N}, space::ProductAxis, detector_range::Number, d::Int) where {T,T2,N}
    F_k = plane_wave_expansion(F, 3) #Calculate PWE
    k_inds = CartesianIndices(F_k)
    k_space = [1/x for x in space]
    for I in k_inds
        k_sq = square(k_space, I)
        if sum(k_sq) ≥ sum(k_sq[setdiff(1:N, d)])
            F_k[I] = zero(T)
        end
        @inbounds F_k[I] = F_k[I]*(-exp(im*k_space[d][I[d]]*detector_range)) #Propagate the wave back to the reciever and reflect
    end
    Eout, Hout = inv_plane_wave_expansion(F_k)[origin(space)...] #Read the value of the wave at the detector
    return covernd([Ezout[3], Hzout[3]] N, d)
end

function loves_ff_transform_x(F::EHTuple{T,N}, Fold::EHTuple, d::Int, Δ::Number) where {T,N}
    E, H = ((Ex,Ez,Ez), (Hx,Hy,Hz)) = F
    Hz_old = Fold[2][3]
    Hy_old = Fold[2][2]
    B = boundary_indices(CartesianIndices(Ex), d)[1][1]
    Ez∞ = zero(T)
    Ey∞ = zero(T)
    for I in B
        Ez∞ += Ez[I] + η₀*(Hy[I]+Hy_old[I])/2
        Ey∞ += Ey[I] + η₀*(Hz[I]+Hz_old[I])/2
    end

    out = Array[Ez∞, Ey∞] .* Δ^2

    return out
end
