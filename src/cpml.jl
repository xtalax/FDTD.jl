################################################################################
# PML update in 2D - see Coefficients in yee.jl if you are confused
################################################################################
###
# pml[z][x] is the effect on x field from pml at z boundary!
# the nodes are grouped differently at the lower boundary

# ΦE[y][x] is backwards from the notation in books - ΦExy
function updatePML(Φ::PMLAux{T}, F::VecArray{T, 2}, C::Coefficients, comp::Int, N::Int) where T
    #allocate arrays
    maxxx =  zeros(T, N, C.s.y.N); maxyy = zeros(T, C.s.x.N, N);
    minxx =  zeros(T, N, C.s.y.N); minyy = zeros(T, C.s.x.N, N);

    x,y,z = (1,2,3)
    maxyx = cover2d(C.bu[comp], y).*Φ.max[y][x]
                    .+ cover2d(C.cu[comp], y) .*selectdim(⊗(C.δ[comp], selectdim(F[z], y, size(F[z],y)-N:size(F[z],y)), y), y, 2:N+1)
    minyx = cover2d(C.bl[comp], y).*Φ.min[y][x]
                    .+ cover2d(C.cl[comp], y) .*selectdim(⊗(C.δ[comp], selectdim(F[z], y, 1:N+1), y), y, 1:N)

    x,y,z = (2,3,1)

    maxxy = cover2d(C.bu[comp], z).*Φ.max[z][x]
                    .+ cover2d(C.cu[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, (size(F[y],z)-N):size(F[y],z)), z)), z, 2:N+1)
    minxy = cover2d(C.bl[comp], z).*Φ.min[z][x]
                    .+ cover2d(C.cl[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, 1:N+1), z)), z, 1:N)
    x,y,z = (3,1,2)
    maxxz = cover2d(C.bu[comp], y).*Φ.max[y][x]
                    .+ cover2d(C.cu[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, (size(F[z],y)-N):size(F[z],y)), y)), y, 2:N+1)
    minxz = cover2d(C.bl[comp], y).*Φ.min[y][x]
                    .+ cover2d(C.cl[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, 1:N+1), y)), y, 1:N)

    maxyz = cover2d(C.bu[comp], z).*Φ.max[z][x]
                    .+ cover2d(C.cu[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, (size(F[y],z)-N):size(F[y],z)), z)), z, 2:N+1)
    minyz = cover2d(C.bl[comp], z).*Φ.min[z][x]
                    .+ cover2d(C.cl[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, 1:N+1), z)), z, 1:N)

    return PMLAux([VecField3(maxxx, maxxy, maxxz), VecField3(maxyx, maxyy, maxyz)],
                      [VecField3(minxx, minxy, minxz), VecField3(minyx, minyy, minyz)])
end

function applyPML(F::VecArray{T,2}, Φ::PMLAux{T}, C::Coefficients, comp::Int, N::Int) where T
    Fx = deepcopy(F[1]); Fy = deepcopy(F[2]); Fz = deepcopy(F[3])

    x,y,z = (1,2,3)
    selectdim(Fx, y, 1:N) .= selectdim(F[x], y, 1:N)
                                    .+ C.s.x.Δ .* selectdim(C[comp][x+3], y, 1:N).*Φ.min[y][x]
    selectdim(Fx, y, (size(F[x],y)-N+1):size(F[x],y)) .= selectdim(F[x], y, (size(F[x],y)-N+1):size(F[x],y))
                                                            .+  C.s.x.Δ .* selectdim(C[comp][x+3], y, (size(F[x],y)-N+1):size(F[x],y)).*Φ.max[y][x]

   x,y,z = (2,3,1)

   selectdim(Fy, z, 1:N) .= selectdim(F[x], z, 1:N)
                               .-  C.s.x.Δ .* selectdim(C[comp][x+3], z, 1:N).*Φ.min[z][x]
   selectdim(Fy, z, (size(F[x],z)-N+1):size(F[x],z)) .= selectdim(F[x], z, (size(F[x],z)-N+1):size(F[x],z))
                                                           .-  C.s.x.Δ .*selectdim(C[comp][x+3], z, (size(F[x],z)-N+1):size(F[x],z)).*Φ.max[z][x]
   x,y,z = (3,1,2)
   selectdim(Fz, y, 1:N) .= selectdim(F[x], y, 1:N)
                                   .+ C.s.x.Δ .* selectdim(C[comp][x+3], y, 1:N).*Φ.min[y][x]
   selectdim(Fz, y, (size(F[x],y)-N+1):size(F[x],y)) .= selectdim(F[x], y, (size(F[x],y)-N+1):size(F[x],y))
                                                           .+  C.s.x.Δ .* selectdim(C[comp][x+3], y, (size(F[x],y)-N+1):size(F[x],y)).*Φ.max[y][x]

   selectdim(Fz, z, 1:N) .= selectdim(F[x], z, 1:N)
                               .-  C.s.x.Δ .* selectdim(C[comp][x+3], z, 1:N).*Φ.min[z][x]
   selectdim(Fz, z, (size(F[x],z)-N+1):size(F[x],z)) .= selectdim(F[x], z, (size(F[x],z)-N+1):size(F[x],z))
                                                           .-  C.s.x.Δ .*selectdim(C[comp][x+3], z, (size(F[x],z)-N+1):size(F[x],z)).*Φ.max[z][x]
    return (Fx, Fy, Fz)
end



function updatePML(Φ::PMLAux{T}, F::VecArray{T,3}, C::Coefficients, comp::Int, N::Int) where T
   #allocate arrays
   maxxx =  zeros(T, N, C.s.y.N, C.s.z.N); maxyy = zeros(T, C.s.x.N, N, C.s.z.N); maxzz = zeros(T, C.s.x.N, C.s.y.N, N);
   minxx =  zeros(T, N, C.s.y.N, C.s.z.N); minyy = zeros(T, C.s.x.N, N, C.s.z.N); minzz = zeros(T, C.s.x.N, C.s.y.N, N);


   x,y,z = (1,2,3)
   maxyx = cover3d(C.bu[comp], y).*Φ.max[y][x]
                   .+ cover3d(C.cu[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, (size(F[z])[y]-N):size(F[z])[y]), y)), y, 2:N+1)
   minyx = cover3d(C.bl[comp], y).*Φ.min[y][x]
                   .+ cover3d(C.cl[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, 1:N+1), y)), y, 1:N)

   maxzx = cover3d(C.bu[comp], z).*Φ.max[z][x]
                   .+ cover3d(C.cu[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, (size(F[y])[z]-N):size(F[y])[z]), z)), z, 2:N+1)
   minzx = cover3d(C.bl[comp], z).*Φ.min[z][x]
                   .+ cover3d(C.cl[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, 1:N+1), z)), z, 1:N)

   x,y,z = (2,3,1)
   maxzy = cover3d(C.bu[comp], y).*Φ.max[y][x]
                   .+ cover3d(C.cu[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, (size(F[z])[y]-N):size(F[z])[y]), y)), y, 2:N+1)
   minzy = cover3d(C.bl[comp], y).*Φ.min[y][x]
                   .+ cover3d(C.cl[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, 1:N+1), y)), y, 1:N)

   maxxy = cover3d(C.bu[comp], z).*Φ.max[z][x]
                   .+ cover3d(C.cu[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, (size(F[y])[z]-N):size(F[y])[z]), z)), z, 2:N+1)
   minxy = cover3d(C.bl[comp], z).*Φ.min[z][x]
                   .+ cover3d(C.cl[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, 1:N+1), z)), z, 1:N)
   x,y,z = (3,1,2)
   maxxz = cover3d(C.bu[comp], y).*Φ.max[y][x]
                   .+ cover3d(C.cu[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, (size(F[z])[y]-N):size(F[z])[y]), y)), y, 2:N+1)
   minxz = cover3d(C.bl[comp], y).*Φ.min[y][x]
                   .+ cover3d(C.cl[comp], y) .*selectdim(( ⊗(C.δ[comp], selectdim(F[z], y, 1:N+1), y)), y, 1:N)

   maxyz = cover3d(C.bu[comp], z).*Φ.max[z][x]
                   .+ cover3d(C.cu[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, (size(F[y])[z]-N):size(F[y])[z]), z)), z, 2:N+1)
   minyz = cover3d(C.bl[comp], z).*Φ.min[z][x]
                   .+ cover3d(C.cl[comp], z) .*selectdim(( ⊗(C.δ[comp], selectdim(F[y], z, 1:N+1), z)), z, 1:N)

   return PMLAux{T}([VecField3(maxxx, maxxy, maxxz), VecField3(maxyx, maxyy, maxyz), VecField3(maxzx, maxzy, maxzz)],
                     [VecField3(minxx, minxy, minxz), VecField3(minyx, minyy, minyz), VecField3(minzx, minzy, minzz)])
end
################################################################################
# PML update in 3d - see Coefficients in yee.jl if you are confused
################################################################################
###
# pml[z][x] is the effect on x field from pml at z boundary!
# the nodes are grouped differently at the lower boundary

# ΦE[y][x] is backwards from the notation in books - ΦExy

function applyPML(F::VecArray{T,3}, Φ::PMLAux{T}, C::Coefficients, comp::Int, N::Int) where T
   Fx = deepcopy(F[1]); Fy = deepcopy(F[2]); Fz = deepcopy(F[3])

   x,y,z = (1,2,3)
   selectdim(Fx, y, 1:N) .= selectdim(F[x], y, 1:N)
                                   .+ C.s.x.Δ .* selectdim(C[comp][x+3], y, 1:N).*Φ.min[y][x]
   selectdim(Fx, y, (size(F[x],y)-N+1):size(F[x],y)) .= selectdim(F[x], y, (size(F[x],y)-N+1):size(F[x],y))
                                                           .+  C.s.x.Δ .* selectdim(C[comp][x+3], y, (size(F[x],y)-N+1):size(F[x],y)).*Φ.max[y][x]

   selectdim(Fx, z, 1:N) .= selectdim(F[x], z, 1:N)
                               .-  C.s.x.Δ .* selectdim(C[comp][x+3], z, 1:N).*Φ.min[z][x]
   selectdim(Fx, z, (size(F[x],z)-N+1):size(F[x],z)) .= selectdim(F[x], z, (size(F[x],z)-N+1):size(F[x],z))
                                                           .-  C.s.x.Δ .*selectdim(C[comp][x+3], z, (size(F[x],z)-N+1):size(F[x],z)).*Φ.max[z][x]

  x,y,z = (2,3,1)
  selectdim(Fy, y, 1:N) .= selectdim(F[x], y, 1:N)
                                  .+ C.s.x.Δ .* selectdim(C[comp][x+3], y, 1:N).*Φ.min[y][x]
  selectdim(Fy, y, (size(F[x],y)-N+1):size(F[x],y)) .= selectdim(F[x], y, (size(F[x],y)-N+1):size(F[x],y))
                                                          .+  C.s.x.Δ .* selectdim(C[comp][x+3], y, (size(F[x],y)-N+1):size(F[x],y)).*Φ.max[y][x]

  selectdim(Fy, z, 1:N) .= selectdim(F[x], z, 1:N)
                              .-  C.s.x.Δ .* selectdim(C[comp][x+3], z, 1:N).*Φ.min[z][x]
  selectdim(Fy, z, (size(F[x],z)-N+1):size(F[x],z)) .= selectdim(F[x], z, (size(F[x],z)-N+1):size(F[x],z))
                                                          .-  C.s.x.Δ .*selectdim(C[comp][x+3], z, (size(F[x],z)-N+1):size(F[x],z)).*Φ.max[z][x]
  x,y,z = (3,1,2)
  selectdim(Fz, y, 1:N) .= selectdim(F[x], y, 1:N)
                                  .+ C.s.x.Δ .* selectdim(C[comp][x+3], y, 1:N).*Φ.min[y][x]
  selectdim(Fz, y, (size(F[x],y)-N+1):size(F[x],y)) .= selectdim(F[x], y, (size(F[x],y)-N+1):size(F[x],y))
                                                          .+  C.s.x.Δ .* selectdim(C[comp][x+3], y, (size(F[x],y)-N+1):size(F[x],y)).*Φ.max[y][x]

  selectdim(Fz, z, 1:N) .= selectdim(F[x], z, 1:N)
                              .-  C.s.x.Δ .* selectdim(C[comp][x+3], z, 1:N).*Φ.min[z][x]
  selectdim(Fz, z, (size(F[x],z)-N+1):size(F[x],z)) .= selectdim(F[x], z, (size(F[x],z)-N+1):size(F[x],z))
                                                          .-  C.s.x.Δ .*selectdim(C[comp][x+3], z, (size(F[x],z)-N+1):size(F[x],z)).*Φ.max[z][x]
   return (Fx, Fy, Fz)
end
