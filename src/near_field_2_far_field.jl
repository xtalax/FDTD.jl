#Plane wave expansion

function plane_wave_expansion(F::EHTuple{T,N}) where T
    F_k = EHTuple{Complex{T}}(size(F))
    for (i,field) in F
        for (j,component) in field
            F_k[i][j] .= fft(component)./length(F)
        end
    end
    return F_k
end
function inv_plane_wave_expansion(F_k::EHTuple{Complex{T},N}) where T
    F = EHTuple{T}(size(F_k))
    for (i,field) in F_k
        for (j,component) in field
            F[i][j] .= abs.(ifft(component).*length(F_k))
        end
    end
    return F
end


function compact_range_transform(F::EHTuple{T,N}, d::Int) where {T,N}
    broadcastable_sum(F) = sum(F, dims=setdiff(1:N, d)) # anonymous function to broadcast over

    mid = size(F, d) ÷ 2
    inds = Array{Union{Colon, AbstractRange}}(fill(Colon(), N))
    inds[d] = mid+1:end

    F_k = plane_wave_expansion(F) #Calculate PWE
    upper_half_plane = selectdim(CartesianIndices(F_k), d, inds...)
    for I in upper_half_plane
        @inbounds F_k[I] = EHElement(zero(Complex{T}))
    end
    k_space = Space(1/x for x in space)

    function backpropagate(F, d)
        for I in selectdim(CartesianIndices(F), d, inds...)
            @inbounds F_k[I] = deep_times(F_k[I], (-exp(im*k[d][I[d]]*r))) #Propagate the wave back to the reciever and reflect
        end
    end

    Fnew = inv_plane_wave_expansion(F_k)
    Fnew = broadcast(broadcastable_sum, Fnew)

    return $$ What $$
end
