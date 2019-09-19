function δ₂(F::AbstractVector{N}, axis::DAxis) where N <: Number
    δ = spsimilar(F)
    for i in axis.i[2:end-1]
        δ[i] = 2*F[i]-F[i-1]-F[i-2]/axis
    end
    return δ
end

function δ⁻(size::Int, step::Number) # reverse difference matrix generator
    δ = spzeros(size,size)
    for i in 1:size
        δ[i,i] = 1
        if i > 1
            δ[i,i-1] = -1
        end
    end
    return (1/step).*δ
end

function δ⁺(F::AbstractVector{N}, axis::DAxis) where N <: Number # forward difference matrix generator
    δ = similar(F)
    for i in axis.i[1:end-1]
        δ[i] = (F[i+1]-F[i])/axis.Δ[i]
    end
    δ[end] = -F[end]
    return δ
end

function δ⁻(F::AbstractVector{N}, axis::DAxis) where N <: Number
    δ = similar(F)
    for i in axis.i[2:end]
        δ[i] = (F[i]-F[i-1])/axis.Δ[i]
    end
    return δ
end



function δ⁺(x::LinearAxis)
    δ = spzeros(x.N,x.N)
    for i in 1:x.N
        δ[i,i] = -1
        if i < x.N
            δ[i,i+1] = 1
        end
    end
    return δ
end

function δ⁺(x::Int)
    δ = spzeros(x,x)
    for i in 1:x
        δ[i,i] = -1
        if i < x
            δ[i,i+1] = 1
        end
    end
    return δ
end


function δ⁻(x::LinearAxis)
    δ = spzeros(x.N,x.N)
    for i in 1:x.N
        δ[i,i] = 1
        if i > 1
            δ[i,i-1] = -1
        end
    end
    return δ
end


function δ⁻(x::Int)
    δ = spzeros(x,x)
    for i in 1:x
        δ[i,i] = 1
        if i > 1
            δ[i,i-1] = -1
        end
    end
    return δ
end


function δδ(size::Int, step::Number) # computes the 2nd order differentiator matrix
     δ = spzeros(size,size)
     for i in 1:size
         δ[i,i] = -2.0

         if i >= 2
             δ[i-1,i] = δ[i,i-1] = 1.0
         end
     end
    # a = [203.0/45.0,-87.0/5.0,117.0/4.0, -254.0/9.0,33.0/2.0, -27.0/5.0, 137.0/180.0]
     #l = length(a)
     #δ[1,1:l] = a
     # δ[end, end-(l-1):end] = a[end:-1:1]
     return δ./step^2
  end

  @noinline function _slice_rmul!(out::AbstractArray,A::AbstractArray{B,M}, u::AbstractArray{T,N}, dim::Int, pre, post) where {T,B,N,M}
      for J in post
          for I in pre
              out[I,:,J] = A*u[I, :, J]
          end
      end
      return out
  end

  function ⊗(A::AbstractArray{B,2}, u::AbstractArray{T,N}, dim::Int) where {T, B, N,M}
      @assert N != 1
      @assert size(A,2) == size(u,dim) "The second dim of `A` must match the sliced dim of `u`, got $(size(A,2)) $(size(u,dim))"
      out = similar(u)
      _slice_rmul!(out, A, u, dim, CartesianIndices(axes(u)[1:dim-1]), CartesianIndices(axes(u)[(dim+1):end]))
      return out
  end


# this function splits the matrix B in to vectors along the specified direction, and then multiplies each vector by the matrix A
 #=function ⊗(A::AbstractMatrix,B::AbstractArray{N,2}, direction::Int) where N<:Number
     s = size(B)

    if direction == 1
        Bt = transpose(B)
        Ct = zeros(N,s[end:-1:1])
        for i in 1:s[1]
            @views Ct[:,i] = A*Bt[:,i]
        end
        C = transpose(Ct)
    elseif direction == 2
        C = zeros(N,s)
        for i in 1:s[2]
            @views C[:,i] = A*B[:,i]
        end
    else
        throw("Direction given is too large")
    end
    return C
end
function ⊗(δ::AbstractMatrix, F::AbstractArray{T,3}, dim::Int) where T
    s = size(F)
    out = similar(F)
    if dim == 1
        for j in 1:s[3]
            for i in 1:s[2]
                @views out[:,i,j] = δ*F[:,i,j]
            end
        end
    elseif dim == 2
        for j in 1:s[3]
            for i in 1:s[1]
                @views out[i,:,j] = δ*F[i,:,j]
            end
        end
    elseif dim == 3
        for j in 1:s[2]
            for i in 1:s[1]
                @views out[i,j,:] = δ*F[i,j,:]
            end
        end
    end
    return out
end

# this is a test that does exactly the same thing as the function above
 function SliceProduct(A::Array{N,2},B::Array{N,2}, direction::Int) where N<:Number
     s = size(B)
     C = zeros(N,s)
     if direction == 1
         for i in 1:s[1]
            @views C[:,i] = A*B[:,i]
         end
    elseif direction == 2
        for i in 1:s[2]
            @views C[i,:] = A*B[i,:]
        end
    else
        throw("direction too large, choose 1 for through columns and 2 for through rows")
    end
    return C
end
#=
function ⊗(A::CuArray,B::CuArray, direction::Int64)
    s = size(B)
    C = cu(zeros(s))
    if direction == 1
        for i in 1:s[2]
            @inbounds C[:,i] = A*B[:,i]
        end
   elseif direction == 2
       for i in 1:s[1]
           @inbounds C[i,:] = A*B[i,:]
       end
   end
   return C
end
=#
=#
∇(u, δxx, δyy) = ⊗(δxx,u,1).+⊗(δyy,u,2)
