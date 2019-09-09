################################################################################
# Update equations in 2D
################################################################################

function update_H(E::VecArray{T,2}, H::VecArray{T,2}, C::Coefficients{T,2}) where T
    s = size(H[1])
    Hx = zeros(T, size(H[1]))
    Hy = zeros(T, size(H[2]))
    Hz = zeros(T, size(H[3]))

    Threads.@threads for j in 1:(s[2]-1)
        for i in 1:(s[1]-1)
            @inbounds Hx[i,j] = C.H[1][i,j]*H[1][i,j] + C.H[4][i,j]*(- C.Km[2][j]*(E[3][i, j+1] - E[3][i, j]))
            @inbounds Hy[i,j] = C.H[2][i,j]*H[2][i,j] + C.H[5][i,j]*(C.Km[1][i]*(E[3][i+1, j] - E[3][i, j]))
            @inbounds Hz[i,j] = C.H[3][i,j]*H[3][i,j] + C.H[6][i,j]*(C.Km[2][j]*(E[1][i, j+1] - E[1][i, j]) - C.Km[1][i]*(E[2][i+1, j] - E[2][i, j]))
        end
    end
    for h in (Hx,Hy,Hz) #apply dirichlet boundary conditions
        for i in 1:2
            selectdim(h, i, size(h,i)) .= zeros(T, perpsize(h,i)...)
            selectdim(h, i, 1) .= zeros(T, perpsize(h,i)...)
        end
    end
    for (i,h) in enumerate((Hx,Hy,Hz)) #truncating the grid for the overhanging components
        for j in setdiff(1:2, i)
            selectdim(h, j, size(h,j)-1) .= zeros(T, perpsize(h,j)...)
        end
    end
    return (Hx, Hy, Hz)
end

function update_E(E::VecArray{T,2}, H::VecArray{T,2}, C::Coefficients{T,2}) where T
    s = size(E[1])
    Ex = zeros(T, size(H[1]))
    Ey = zeros(T, size(H[2]))
    Ez = zeros(T, size(H[3]))
    Threads.@threads for j in 2:s[2]
        for i in 2:s[1]
            @inbounds Ex[i,j] = C.E[1][i,j]*E[1][i,j] + C.E[4][i,j]*(C.Ke[2][j]*(H[3][i,j] - H[3][i,j-1]))
            @inbounds Ey[i,j] = C.E[2][i,j]*E[2][i,j] + C.E[5][i,j]*(-C.Ke[1][i]*(H[3][i, j] - H[3][i-1, j]))
            @inbounds Ez[i,j] = C.E[3][i,j]*E[3][i,j] + C.E[6][i,j]*(C.Ke[1][i]*(H[2][i,j] - H[2][i-1,j]) - C.Ke[2][j]*(H[1][i, j] - H[1][i, j-1]))
        end
    end
    for e in (Ex,Ey,Ez) #apply dirichlet boundary conditions
        for i in 1:2
            selectdim(e, i, 1) .= zeros(T, perpsize(e,i)...)
            selectdim(e, i, size(e,i)) .= zeros(T, perpsize(e,i)...)
        end
    end
    for (i,e) in enumerate((Ex,Ey))
        selectdim(e, i, size(e,i)-1) .= zeros(T, perpsize(e,i)...)
    end
    return (Ex, Ey, Ez)
end

################################################################################
# Update equations in 3D
################################################################################

function update_H(E::VecArray{T,3}, H::VecArray{T,3}, C::Coefficients{T,3}) where T
    s = size(H[1])
    Hx = zeros(T, size(H[1]))
    Hy = zeros(T, size(H[2]))
    Hz = zeros(T, size(H[3]))

    Threads.@threads for k in 1:(s[3]-1)
        for j in 1:(s[2]-1), i in 1:(s[1]-1)
        @inbounds Hx[i,j,k] = C.H[1][i,j,k]*H[1][i,j,k] + C.H[4][i,j,k]*(C.Km[3][k]*(E[2][i, j, k+1] - E[2][i, j, k]) - C.Km[2][j]*(E[3][i, j+1, k] - E[3][i, j, k]))
        @inbounds Hy[i,j,k] = C.H[2][i,j,k]*H[2][i,j,k] + C.H[5][i,j,k]*(C.Km[1][i]*(E[3][i+1, j, k] - E[3][i, j, k]) - C.Km[3][k]*(E[1][i, j, k+1] - E[1][i, j, k]))
        @inbounds Hz[i,j,k] = C.H[3][i,j,k]*H[3][i,j,k] + C.H[6][i,j,k]*(C.Km[2][j]*(E[1][i, j+1, k] - E[1][i, j, k]) - C.Km[1][i]*(E[2][i+1, j, k] - E[2][i, j, k]))
    end
    for h in (Hx,Hy,Hz) #apply dirichlet0 boundary conditions
        for i in 1:3
            selectdim(h, i, size(h,i)) .= zeros(T, perpsize(h,i)...)
            selectdim(h, i, 1) .= zeros(T, perpsize(h,i)...)
        end
    end
    for (i,h) in enumerate((Hx,Hy,Hz)) #truncating the grid for the overhanging components - implemented this way for simplicity
        for j in setdiff(1:3, i)
            selectdim(h, j, size(h,j)-1) .= zeros(T, perpsize(h,j)...)
        end
    end
    return (Hx, Hy, Hz)
end

function update_E(E::VecArray{T,3}, H::VecArray{T,3}, C::Coefficients{T,3}) where T
    s = size(E[1])
    Ex = zeros(T, size(H[1]))
    Ey = zeros(T, size(H[2]))
    Ez = zeros(T, size(H[3]))
    Threads.@threads for k in 2:s[3],
        for j in 2:s[2], i in 2:s[1]
        @inbounds Ex[i,j,k] = C.E[1][i,j,k]*E[1][i,j,k] + C.E[4][i,j,k]*(C.Ke[2][j]*(H[3][i,j,k] - H[3][i,j-1,k]) - C.Ke[3][k]*(H[2][i, j, k] - H[2][i, j, k-1]))
        @inbounds Ey[i,j,k] = C.E[2][i,j,k]*E[2][i,j,k] + C.E[5][i,j,k]*(C.Ke[3][k]*(H[1][i,j,k] - H[1][i,j,k-1]) - C.Ke[1][i]*(H[3][i, j, k] - H[3][i-1, j, k]))
        @inbounds Ez[i,j,k] = C.E[3][i,j,k]*E[3][i,j,k] + C.E[6][i,j,k]*(C.Ke[1][i]*(H[2][i,j,k] - H[2][i-1,j,k]) - C.Ke[2][j]*(H[1][i, j, k] - H[1][i, j-1, k]))
    end
    for e in (Ex,Ey,Ez) #apply dirichlet0 boundary conditions
        for i in 1:3
            selectdim(e, i, 1) .= zeros(T, perpsize(e,i)...)
            selectdim(e, i, size(e,i)) .= zeros(T, perpsize(e,i)...)
        end
    end
    for (i,e) in enumerate((Ex,Ey,Ez))
        selectdim(e, i, size(e,i)-1) .= zeros(T, perpsize(e,i)...)
    end
    return (Ex, Ey, Ez)
end
