pwd() ∉ LOAD_PATH && push!(LOAD_PATH, pwd())
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src")
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\TMFDTD" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\TMFDTD")

"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
using FDTD_2D
using DiscreteAxis
using DiffEqOperators

function epsilon(r::AbstractVector{T}) where T
    d = [0.305, 0.4]
    if all(0.0.-d./2 .<= r .<= d./2)
        return convert(T, 2.3)
    else
        return one(T)
    end
end
dB(x::Number) = 20*log10(abs(x))

function dipole_test2D(θ=0.0)
    fmax = 0.5*10^9
    f₀ = fmax/2
    nPML = 14
    space, time = setup_spacetime(fmax, 0.5, nPML;
                                    xlims = (-1.0,1.0),
                                    ylims = (-1.0,1.0))
    println("$(space.x.N), $(space.y.N), $time.N")

    println("Δt is $(time.Δ)")

    ###############################################################################
    # MEDIUM DEFINITION
    ###############################################################################

    @time medium = Medium(space; θ = θ, ε = epsilon) #freespace medium
    Esource = sin.(2π*f₀.*time)
    Hsource = zeros(time.N)

    sourceloc = [x.N ÷ 2 for x in space]
    sourceindex = (2, sourceloc)
    detectorindex = (:,[:,:])

    @time field = FDTD_propagate(space, time, f₀, nPML; #Try with a PML
                                    source = (Esource, Hsource),
                                    medium = medium,
                                    sourceindex = sourceindex,
                                    detectorindex = detectorindex
                                    )



end
    using DiffEqOperators
function diff_eq_wave()

function main()
s = x, y, z= (-5:0.2:5, -5:0.2:5,  -5:0.2:5)
dx = dy = dz = x[2] - x[1]

Dxx = CenteredDifference{1}(2, 4, dx, length(x))
Dyy = CenteredDifference{2}(2, 4, dy, length(y))

A = Dxx+Dyy+Dzz
Q = compose(Dirichlet0BC(Float64, length.(s))...)

dt = dx/(sqrt(3)*3e8)
t = 0.0:dt:10/3e8

f(u,p,t) = (3e8)^2 .*(A*Q*u)

uolder = deepcopy(u0)
uold = deepcopy(u0)
u = deepcopy(u0)

function steptime(u,uold,uolder)
    return ((dt^2 .*f(u,0,0) .+ 5u .- 4uold .+ uolder)./2, u, uold)
end
u,uold,uolder = steptime(u0,uold,uolder)
prog = Progress(length(t))
@gif for ti in eachindex(t) #4th order time stepper
        u,uold,uolder = steptime(u,uold,uolder)
        heatmap(u)
    update!(prog, ti)
end

end
main()

function Medium_Test2D(θ=0.0)
    fmax = 0.5*10^9
    f₀ = fmax/2
    nPML = 14
    space, time = setup_spacetime(fmax, 1.5, nPML;
                                    xlims = (-1.0,1.0),
                                    ylims = (-1.0,1.0))
    println("$(space.x.N), $(space.y.N), $time.N")

    println("Δt is $(time.Δ)")

    ###############################################################################
    # MEDIUM DEFINITION
    ###############################################################################

    @time medium = Medium(space; θ = θ, ε = epsilon) #freespace medium
    Esource = sinc.(4π.*fmax.*(time.-maximum(time.pts)/4))

    Hsource = zeros(time.N)
    interior = interior_range(space, nPML)
    sourceloc = [nPML+2, interior[2]]
    sourceindex = (2, sourceloc)
    detectorindex = (:,[:,:,:])

    @time field = FDTD_propagate(space, time, f₀, nPML; #Try with a PML
                                    source = (Esource, Hsource),
                                    medium = medium,
                                    sourceindex = sourceindex,
                                    detectorindex = detectorindex
                                    )
    #s_size = size(space)


    return field

end
field = Medium_Test2D(0.1)

using Plots
plotly()
heatmap(dB.(field[200][1][3]))
