pwd() ∉ LOAD_PATH && push!(LOAD_PATH, pwd())
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src")
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\TMFDTD" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\TMFDTD")

"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
using FDTD_2D
using DiscreteAxis

function epsilon(r::AbstractVector{T}) where T
    d = [0.305, 0.4]
    if all(0.0.-d./2 .<= r .<= d./2)
        return convert(T, 2.3)
    else
        return one(T)
    end
end
dB(x::Number) = 20*log10(abs(x))

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
