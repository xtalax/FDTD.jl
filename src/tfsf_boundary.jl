#TFSF - at the moment it is fixed to impinging from the -x direction

struct TFSF{T}
    Ez::Array{T,2}
    Hy::Array{T,2}
    d::Int
    N::Int
end

function TFSF(space::LinearAxis, time::LinearAxis, source::AbstractVector{T}, N::Int) where T
      Ez = zeros(T, space.N)
      Hy = zeros(T, space.N)

      F = (Ez,Hy)
      Esource, Hsource = source
      C = Coefficients(space, time)

      temp = sqrt(C.E[6][1]*C.H[5][1])
      abcCoefLeft = (temp - 1.0) / (temp + 1.0)
      temp = sqrt(C.E[6][end]*C.H[5][end-2])
      abcCoefRight = (temp - 1.0) / (temp + 1.0)

      abc = (abcCoefLeft, abcCoefRight)

      Ezout = Matrix{T}(undef, space.N, time.N)
      Hyout = Matrix{T}(undef, space.N, time.N)

  for i in time.i
      Fnew= step_time(F, abc, C, Esource[i], Hsource[i], sourceindex)
      Ezout[:,i] .= Fnew[1]
      Hyout[:,i] .= Fnew[2]
      F = Fnew;
  end
  return TFSF{T}(Ezout, Hyout, N)
end
function boundary_indices(A::CartesianIndices{N}, d::Int) where N
    S = size(A)
    I = CaretesianIndices(A)
    out = []
    for dim in 1:N
        Ilowview = selectdim(I, dim, d)
        Ihighview = selectdim(out, dim, S[dim]-d+1)
        for (index, otherdim) in enumerate(setdiff(dimset, dim))
            Ilowview = selectdim(Ilowview, index, d:(S[otherdim]-d+1))
            Ihighview = selectdim(Ihighview, index, d:(S[otherdim]-d+1))
        end
        push!(out, (Ilowview, Ihighview))
    end
    return out
end

@pure function unit_indices(N::Int) #create unit CartesianIndex for each dimension
    @assert N>-1
    out = Vector{CartesianIndex{N}}(undef, N)
    null = zeros(Int64, N)
    for i in 1:N
        unit_i = copy(null)
        unit_i[i] = 1
        out[i] = CartesianIndex(Tuple(unit_i))
    end
    Tuple(out)
end

function apply_TFSF_boundary!(F::EHTuple{T,3}, C::Coefficients{T,3} inc::TFSF, ti::Int) where{T}
    ((Ex, Ey,Ez), (Hx, Hy, Hz)) = deepcopy(F)
    boundary = boundary_indices(CartesianIndices(A), 1)

    e = unit_indices(3)

    for I in boundary[1][1]
        Ez[I] = Ey[I] - C.E[6][I] * inc.Hy[I[1], ti]
        Hy[I] = Hy[I-e[1]] - C.H[5][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[2][1]
        Hx[I] = Hx[I-e[i]] + C.H[4][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[3][1]
        Ex[I] = Ex[I] + C.E[4][I] * inc.Hy[I[1], ti]
    end

    for I in boundary[1][2]
        Ez[I] = Ey[I] + C.E[6][I] * inc.Hy[I[1], ti]
        Hy[I] = Hy[I] + C.H[5][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[2][2]
        Hx[I] = Hx[I] - C.H[4][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[3][2]
        Ex[I] = Ex[I] - C.E[4][I] * inc.Hy[I[1], ti]
    end

    return ((Ex, Ey,Ez), (Hx, Hy, Hz))
end

function apply_TFSF_boundary!(F::EHTuple{T,2}, C::Coefficients{T,2} inc::TFSF, ti::Int) where{T}
    ((Ex, Ey,Ez), (Hx, Hy, Hz)) = deepcopy(F)
    boundary = boundary_indices(CartesianIndices(A), 1)

    e = unit_indices(2)

    for I in boundary[1][1]
        Ez[I] = Ey[I] - C.E[6][I] * inc.Hy[I[1], ti]
        Hy[I] = Hy[I-e[1]] - C.H[5][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[2][1]
        Hx[I] = Hx[I-e[i]] + C.H[4][I] * inc.Ez[I[1], ti]
    end

    for I in boundary[1][2]
        Ez[I] = Ey[I] + C.E[6][I] * inc.Hy[I[1], ti]
        Hy[I] = Hy[I] + C.H[5][I] * inc.Ez[I[1], ti]
    end
    for I in boundary[2][2]
        Hx[I] = Hx[I] - C.H[4][I] * inc.Ez[I[1], ti]
    end

    return ((Ex, Ey,Ez), (Hx, Hy, Hz))
end
