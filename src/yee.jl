abstract type Field{T, N} <: FieldVector{N, Array{T, N}} end
abstract type CompositeField{T, J, N} <: FieldVector{J, Field{T,N}} end

const ε₀ = 8.85418782*10^-12
const μ₀ = 4*pi*10^-7
const c₀ = 1/√(ε₀*μ₀)
const η₀ = √(μ₀/ε₀)
const Sc = 1/√2


struct VecField3{T<:Number, N} <: FieldVector{3, Array{T, N}}
    x::Array{T,N}
    y::Array{T,N}
    z::Array{T,N}
    VecField3(x::Array{T,N}, y::Array{T,N}, z::Array{T,N}) where {T,N} = new{T,N}(x,y,z)
    VecField3(s::Space2D{T}) where {T<:Number} =  new{T,2}(zeros(T, s.x.N, s.y.N),
                                                         zeros(T, s.x.N, s.y.N),
                                                         zeros(T, s.x.N, s.y.N))

    function VecField3(s::Space2D{T}, nLayer::Int, dim::Int) where {T<:Number}
        if dim == 1
            return new{T,2}(zeros(T, nLayer, s.y.N),
                            zeros(T, nLayer, s.y.N),
                            zeros(T, nLayer, s.y.N))
        elseif dim == 2
            return new{T,2}(zeros(T, s.x.N, nLayer),
                          zeros(T, s.x.N, nLayer),
                          zeros(T, s.x.N, nLayer))
        else
            throw("Only Integer dimension of 1-3 is supported! Check the dim argument.\n")
        end
    end
end

Base.size(field::Field{T,N}) where {T,N} = (size(field[1]), N)

colons(n::Int) = fill(:, n)

function flatten(F::Field{T,N}) where {T,N}
    flatfield = zeros(T, size(F)...)
    for (i,component) in enumerate(Field)
        flatfield[i, colons(N)...] .= component
    end
    return flatfield
end
Base.eltype(A::VecField3) = eltype(A[1])

struct PMLAux{T<:Number,N}
    max::Vector{VecField3{T}}
    min::Vector{VecField3{T}}
    function PMLAux(s::Space2D{T}, nLayer) where {T<:Number}
        max = [VecField3(s, nLayer, 1), VecField3(s, nLayer, 2)]
        min = [VecField3(s, nLayer, 1), VecField3(s, nLayer, 2)]
        new{T,2}(max,min)
    end
    PMLAux(max::Vector{VecField3{T,N}}, min::Vector{VecField3{T,N}}) where {T,N} = new{T,N}(max,min)
end

struct PML{T<:Number,M}
    E::PMLAux{T,M}
    H::PMLAux{T,M}
    N::Int64
    function PML(s::Space2D{T}, nLayer) where {T<:Number}
        E, H = (PMLAux(s, nLayer+1), PMLAux(s,nLayer+1))
        new{T,2}(E,H, nLayer+1)
    end
    PML(E::PMLAux{T,M}, H::PMLAux{T,M},N::Int64) where {T,M} = new{T,M}(E,H,N)
end

struct EMField{T<:Number,N} <:FieldVector{2, VecField3{T}}
    E::VecField3{T,N}
    H::VecField3{T,N}
    EMField(E::VecField3{T,N}, H::VecField3{T,N}) where {T,N} = new{T,N}(E,H)
    EMField(space::Space2D{T}) where {T<:Number}= new{T,2}(VecField3(space), VecField3(space))
    EMField(F::Tuple{Tuple{Array{T,N},Array{T,N},Array{T,N}},Tuple{Array{T,N},Array{T,N},Array{T,N}}}) where {T,N} = new{T,N}(VecField3(F[1]...), VecField3(F[2]...))
end

    function split_array(A::AbstractArray{T}) where {T <: AbstractVector}
        out = [zeros(eltype(first(A)), size(A)...) for i in 1:length(first(A))]
        for I in CartesianIndices(A)
            for (target, n) in enumerate(out)
                target[I] = A[I][n]
            end
        end
        return Tuple(out)
    end

    function vec_array(A::NTuple{N, T}) where {T <: AbstractArray, N}
        @assert all([size(A[i]) == size(A[1]) for i in 2:N])
        @assert all([eltype(A[i]) == eltype(A[1]) for i in 2:N])
        out = Array{SVector{N, T}, ndims(A[1])}(undef, size(A[1])...)
        for I in CartesianIndices{A}
            out[I] = SVector([A[i][I] for i in 1:N])
        end
        return out
    end


    power_density(F::EMField) = @. (F.E.x^2 + F.E.y^2 + F.E.z^2)/(2*η₀)
    function pontying(F::EMField{T}) where T
        out = Array{SVector{3,T}}(undef, size(F)...)
        for I in CartesianIndices(F.E)
            out[I] = SVector([F.E[2][I]*F.H[3][I] - F.E[3][I]*F.H[2][I],
                 F.E[3][I]*F.H[1][I] - F.E[1][I]*F.H[3][I],
                 F.E[1][I]*F.H[2][I] - F.E[2][I]*F.H[1][I]])

        end
        return out
    end
function deep_times(F::VecArray, k::Number)
    Fx = F[1].*k
    Fy = F[2].*k
    Fz = F[3].*k
    return (Fx,Fy,Fz)
end

function Base.getproperty(F::EHTuple, s::Symbol)
    if symbol == :E
        return F[1]
    elseif symbol == :H
        return F[2]
    elseif symbol == :D
        return deep_times(F[1], ε₀)
    elseif symbol == :W
        return power_density(F)
    elseif symbol == :S
        return pontying(F)
    else
        throw(ArgumentError("Symbol $s is not a property of F"))
    end
end

function Base.getproperty(F::VecArray, s::Symbol)
    if symbol == :x
        return F[1]
    elseif symbol == :y
        return F[2]
    elseif symbol == :z
        return F[3]
    else
        throw(ArgumentError("Symbol $s is not a property of F"))
    end
end

Base.getindex(F::EHTuple{T,2}, i::Int, j::Int) where T = (SVector(F.E.x[i,j], F.E.y[i,j], F.E.z[i,j]), SVector(F.H.x[i,j], F.H.y[i,j], F.H.z[i,j]))

Base.getindex(F::EHTuple{T,3}, i::Int, j::Int, k::Int) where T = (SVector(F.E.x[i,j,k], F.E.y[i,j,k], F.E.z[i,j,k]), SVector(F.H.x[i,j,k], F.H.y[i,j,k], F.H.z[i,j,k]))
Base.getindex(F::EHTuple, I::CartesianIndex) = (SVector(F.E.x[I], F.E.y[I], F.E.z[I]), SVector(F.H.x[I], F.H.y[I], F.H.z[I]))
Base.getindex(F::EHTuple, I::CartesianIndex, i::Int) = SVector(F.E.x[I], F.E.y[I], F.E.z[I])
power_density(x::AbstractVector) = sum(x.^2)/2*η₀

Base.CartesianIndices(F::EMField) = CartesianIndices(F.H.x)
Base.size(F::EMField) = size(F.E.x)
με_default(r::AbstractVector) = one(Float64)
σ_default(r::AbstractVector) = zero(Float64)

struct Medium{T<:Number, N} <: Field{T,4}
    ε::Array{T, N}
    μ::Array{T, N}
    σ::Array{T, N}
    σm::Array{T, N}
    function Medium(s::Space2D{T}; θ=0.0, rotcentre = [0.0,0.0], ε = με_default, μ = με_default, σ = σ_default, σm = σ_default) where {T<:Number}
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)] #generate rotation matrix
        X, Y= (s.x.pts .- rotcentre[1], s.y.pts .- rotcentre[2])
        r = [SVector{2,T}(R*[x,y]) for x in X, y in Y]
        eps = ε₀.*ε.(r)
        mu = μ₀.*μ.(r)
        sig = σ.(r)
        sigm = σm.(r)
        new{eltype(eps),2}(eps, mu, sig, sigm)
    end
    function Medium(s::Space3D{T}; azimuth = 0.0, elevation = 0.0, rotcentre = [0.0,0.0,0.0], ε = με_default, μ = με_default, σ = σ_default, σm = σ_default) where {T<:Number}
        R = RotZXY(azimuth, elevation, 0.0) #generate rotation matrix
        X, Y, Z = (s.x.pts .- rotcentre[1], s.y.pts .- rotcentre[2], s.z.pts .- rotcentre[3])
        r′ = [SVector{3,T}(R*[x,y,z]) for x in X, y in Y, z in Z]
        eps = ε₀.*ε.(r′)
        mu = μ₀.*μ.(r′)
        sig = σ.(r′)
        sigm = σm.(r′)
        new{eltype(eps),3}(eps, mu, sig, sigm)
    end
end

struct Coefficients{T,N} <: FieldVector{11, UnionAll}
    E::NTuple{6, Array{T, N}}
    H::NTuple{6, Array{T, N}}
    s::Space2D{T}
    t::LinearAxis{T}
    δ::NTuple{2, Matrix{T}}
    bu::NTuple{2,Vector{T}}
    bl::NTuple{2,Vector{T}}
    cu::NTuple{2,Vector{T}}
    cl::NTuple{2,Vector{T}}
    Km::NTuple{2,Vector{T}}
    Ke::NTuple{2,Vector{T}}
end

    # This function assumes that Δx Δy and Δz are the same - Not too keen on the idea of implementing a non uniform grid right now

function Coefficients(s::Space2D, t, m::Medium{T,2}, f₀, d::Int) where T
    δtmp = (δ⁺(d+2), δ⁻(d+2))
    δ = map(D -> D./s.x.Δ, δtmp)
    # Calculate Coefficients for the layer
    layer = (0:d)

    kmax = 15
    M = 4
    Ma = 1
    #https://www.degruyter.com/downloadpdf/j/jee.2017.68.issue-1/jee-2017-0006/jee-2017-0006.pdf !!!!!!
    σemax = (M+1)*0.8/(η₀*s.x.Δ)
    σmmax = (M+1)*0.8*η₀/s.x.Δ
    aemax = 2π*f₀*ε₀/10
    ammax = 2π*f₀*μ₀/10
    for (i, z) in enumerate(PML)
        kmu[i] = 1+(kmax-1)*((z+0.5)/d)^M
        kml[i] = 1+(kmax-1)*(z/d)^M
        keu[i] = 1+(kmax-1)*(z/d)^M
        kel[i] = 1+(kmax-1)*((z+0.5)/d)^M
        σmu[i] = σmmax*((z+0.5)/d)^M
        σml[i] = σmmax*(z/d)^M
        σeu[i] = σemax*(z/d)^M
        σel[i] = σemax*((z+0.5)/d)^M
        αmu[i] = ammax*(((d - z+0.5)/d)^Ma)
        αml[i] = ammax*(((d - z)/d)^Ma)
        αeu[i] = aemax*(((d - z)/d)^Ma)
        αel[i] = aemax*(((d - z+0.5)/d)^Ma)
    end

    # Calculate auxillary Coefficients for the layer
    beu = exp.(-(σeu./keu .+ αeu).*t.Δ/ε₀)
    ceu = (beu.-1.0).*σeu ./ (σeu.*keu + keu.^2 .*αeu)
    bel = exp.(-(σel./kel .+ αel).*t.Δ/ε₀)
    cel = (bel.-1.0).*σel ./ (σel.*kel + kel.^2 .*αel)


    bmu = exp.(-(σmu./kmu .+ αmu).*t.Δ/μ₀)
    cmu = (bmu.-1.0).*σmu ./ (σmu.*kmu + kmu.^2 .*αmu)
    bml = exp.(-(σml./kml .+ αml).*t.Δ/μ₀)


    cml = (bml.-1.0).*σml ./ (σml.*kml + kml.^2 .*αml)


    tmp = 1.0 .+ m.σm.*t.Δ ./(2.0.*m.μ)
    #Δt ⪕ √(μ₀ε₀)
    Chxh = (2.0 .- tmp)./tmp #source chapter 9 page 8 https://www.eecs.wsu.edu/~schneidj/ufdtd/, with modifications to support inhomogeneous media
    Chyh = Chzh = Chxh # Start here if you want to implement anisotropy

    Chxe = t.Δ ./ (tmp .*m.μ .* s.x.Δ )
    Chye = Chze = Chxe

    tmp = 1.0 .+ m.σ ./(2.0.*m.ε)
    Cexe = (2.0.-tmp)./tmp
    Ceye = Ceze = Cexe # these are equal in the anisotropic uniform grid case, this kind of assignment also only allocates 1 matrix of memory, meaning the constants are passed as pointers to the same place in memory

    Cexh = t.Δ ./ (tmp .* m.ε.*s.x.Δ)
    Ceyh = Cezh = Cexh

    # put all the coeffs in to a nice neat wrapper
    H = (Chxh, Chyh, Chzh, Chxe, Chye, Chze)
    E = (Cexe, Ceye, Ceze, Cexh, Ceyh, Cezh)

    Ke = Vector[ones(s.x.N),ones(s.y.N)]
    Km = Vector[ones(s.x.N),ones(s.y.N)]
    for i in 1:2
        Km[i][1:(d+1)] .= 1 ./reverse(kml)
        Km[i][(end-d):end] .= 1 ./kmu
        Ke[i][1:d+1] .= 1 ./reverse(kel)
        Ke[i][(end-d):end] .= 1 ./keu
    end
    rbel = reverse(bel)
    rcel = reverse(cel)
    rbml = reverse(bml)
    rcml = reverse(cml)
    Coefficients{T,2}(E, H,s,t, δ, (beu,bmu), (rbel,rbml), (ceu,cmu), (rcel,  rcml), Tuple(Km), Tuple(Ke))
end

    function Coefficients(s, t, m::Medium{T,3}, f₀, d::Int) where {T}
        δtmp =(δ⁺(d+2), δ⁻(d+2))
        δ = map(D -> D./s.x.Δ, δtmp)
        # Calculate Coefficients for the PML
        PML = (0:d)


        kmax = 15
        M = 4
        Ma = 1
        #https://www.degruyter.com/downloadpdf/j/jee.2017.68.issue-1/jee-2017-0006/jee-2017-0006.pdf !!!!!!
        σemax = (M+1)*0.8/(η₀*s.x.Δ)
        σmmax = (M+1)*0.8*η₀/s.x.Δ
        aemax = 2π*f₀*ε₀/10
        ammax = 2π*f₀*μ₀/10


        for (i, z) in enumerate(PML)
            kmu[i] = 1+(kmax-1)*((z+0.5)/d)^M
            kml[i] = 1+(kmax-1)*(z/d)^M
            keu[i] = 1+(kmax-1)*(z/d)^M
            kel[i] = 1+(kmax-1)*((z+0.5)/d)^M
            σmu[i] = σmmax*((z+0.5)/d)^M
            σml[i] = σmmax*(z/d)^M
            σeu[i] = σemax*(z/d)^M
            σel[i] = σemax*((z+0.5)/d)^M
            αmu[i] = ammax*(((d - z+0.5)/d)^Ma)
            αml[i] = ammax*(((d - z)/d)^Ma)
            αeu[i] = aemax*(((d - z)/d)^Ma)
            αel[i] = aemax*(((d - z+0.5)/d)^Ma)
        end


        # Calculate auxillary Coefficients for the PML
        beu = exp.(-(σeu./keu .+ αeu).*t.Δ/ε₀)
        ceu = (beu.-1.0).*σeu ./ (σeu.*keu + keu.^2 .*αeu)
        bel = exp.(-(σel./kel .+ αel).*t.Δ/ε₀)
        cel = (bel.-1.0).*σel ./ (σel.*kel + kel.^2 .*αel)


        bmu = exp.(-(σmu./kmu .+ αmu).*t.Δ/μ₀)
        cmu = (bmu.-1.0).*σmu ./ (σmu.*kmu + kmu.^2 .*αmu)
        bml = exp.(-(σml./kml .+ αml).*t.Δ/μ₀)


        cml = (bml.-1.0).*σml ./ (σml.*kml + kml.^2 .*αml)


        tmp = 1.0 .+ m.σm.*t.Δ ./(2.0.*m.μ)
        #Δt ⪕ √(μ₀ε₀)
        Chxh = (2.0 .- tmp)./tmp #source chapter 9 page 8 https://www.eecs.wsu.edu/~schneidj/ufdtd/, with modifications to support inhomogeneous media
        Chyh = Chzh = Chxh # Start here if you want to implement anisotropy

        Chxe = t.Δ ./ (tmp .*m.μ .* s.x.Δ )
        Chye = Chze = Chxe

        tmp = 1.0 .+ m.σ ./(2.0.*m.ε)
        Cexe = (2.0.-tmp)./tmp
        Ceye = Ceze = Cexe # these are equal in the anisotropic uniform grid case, this kind of assignment also only allocates 1 matrix of memory, meaning the constants are passed as pointers to the same place in memory

        Cexh = t.Δ ./ (tmp .* m.ε.*s.x.Δ)
        Ceyh = Cezh = Cexh

        # put all the coeffs in to a nice neat wrapper
        H = (Chxh, Chyh, Chzh, Chxe, Chye, Chze)
        E = (Cexe, Ceye, Ceze, Cexh, Ceyh, Cezh)

        Ke = Vector[ones(s.x.N),ones(s.y.N),ones(s.z.N)]
        Km = Vector[ones(s.x.N),ones(s.y.N),ones(s.z.N)]
        for i in 1:3
            Km[i][1:(d+1)] .= 1 ./reverse(kml)
            Km[i][(end-d):end] .= 1 ./kmu
            Ke[i][1:d+1] .= 1 ./reverse(kel)
            Ke[i][(end-d):end] .= 1 ./keu
        end
        rbel = reverse(bel)
        rcel = reverse(cel)
        rbml = reverse(bml)
        rcml = reverse(cml)
        Coefficients{T}(E, H,s,t, δ, (beu,bmu), (rbel,rbml), (ceu,cmu), (rcel,  rcml), Tuple(Km), Tuple(Ke))
    end
