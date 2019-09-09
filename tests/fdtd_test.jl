


    using FDTD
    using DiscreteAxis

    function epsilon(r::AbstractVector{T}) where T
        d = [0.305, 0.305, 0.4]
        if all(0.0.-d./2 .<= r .<= d./2)
            return convert(T, 2.3)
        else
            return one(T)
        end
    end
    dB(x::Number) = 20*log10(abs(x))


        function PML_Test(detectorindex)
            f = 0.5*10^9
            f₀ = f/2
            nPML = 20
            space, time = setup_spacetime(f, 1.5, nPML;
                                            xlims = (-1.0,1.0),
                                            ylims = (-1.0,1.0),
                                            zlims = (-1.0,1.0))

            println("$(space.x.N), $(space.y.N), $(space.z.N), $time.N")

            println("Δt is $(time.Δ)")

            ###############################################################################
            # MEDIUM DEFINITION
            ###############################################################################

            @time medium = Medium(space) #freespace medium
            sourceloc = [div(space.x.N, 2), div(space.y.N, 2), div(space.z.N, 2)]
            Esource = [sin(2π*f*T)*exp(-(c₀*T)^2) for T in time]
            Hsource = 0.0

            # at the moment sources/detectors homogoeneously add/read to/from the field at the specified location, support for antennae like structures may be added in later versions with support for a J field
            sourceindex = (2, sourceloc) # source index follows the form (E component to stimulate, [location of source])

            @time field_pml = FDTD_propagate(space, time, f₀, nPML; #Try with a PML
                                            source = (Esource, Hsource),
                                            medium = medium,
                                            sourceindex = sourceindex,
                                            detectorindex = detectorindex
                                            )
            #s_size = size(space)
            bigspace, longtime = setup_spacetime(f, 1.5, 0;
                                            xlims = (-1.5,1.5),
                                            ylims = (-1.5,1.5),
                                            zlims = (-1.5,1.5))

            @time mediumbig = Medium(bigspace) #freespace medium



            @time field_extended = @time field = FDTD_propagate(
                                                                bigspace, time, f₀, nPML; #Try with a PML
                                                                source = (Esource, Hsource),
                                                                medium = mediumbig,
                                                                sourceindex = sourceindex,
                                                                detectorindex = detectorindex
                                                                )

        return field_pml, field_extended

    end

    function Medium_Test(θ=0.0)
        f = 0.5*10^9
        f₀ = f/2
        nPML = 14
        space, time = setup_spacetime(f, 1.5, nPML;
                                        xlims = (-1.0,1.0),
                                        ylims = (-1.0,1.0),
                                        zlims = (-1.0,1.0))

        println("$(space.x.N), $(space.y.N), $(space.z.N), $time.N")

        println("Δt is $(time.Δ)")

        ###############################################################################
        # MEDIUM DEFINITION
        ###############################################################################

        @time medium = Medium(space; azimuth = θ, ε = epsilon) #freespace medium
        Esource = [sin(2π*f*T)*exp(-(c₀*T)^2) for T in time]
        Hsource = zeros(time.N)
        interior = interior_range(space, nPML)
        sourceloc = [nPML+2, interior[2], interior[3]]
        sourceindex = (2, sourceloc)
        detectorindex = (1,2,[interior[1], div(space.x.N,2), interior[3]])

        @time field = FDTD_propagate(space, time, f₀, nPML; #Try with a PML
                                        source = (Esource, Hsource),
                                        medium = medium,
                                        sourceindex = sourceindex,
                                        detectorindex = detectorindex
                                        )
        #s_size = size(space)


    return field

end




    function RCS_Generate_test()
        f = 0.5*10^9
        f₀ = f/2
        nPML = 20
        space, time = setup_spacetime(f, 1.5, nPML;
                                        xlims = (-1.0,1.0),
                                        ylims = (-1.0,1.0),
                                        zlims = (-1.0,1.0))

        ###############################################################################
        # MEDIUM DEFINITION
        ###############################################################################

        @time medium = Medium(space; ε = epsilon) #freespace medium
        interior = interior_range(space, nPML)
        sourceloc = [interior[1], interior[2], nPML+2]
        tmax = maximum(time.pts)
        Esource = [sinc(4π*f*(T-tmax/4)) for T in time]
        Hsource = 0.0

        # at the moment sources/detectors homogoeneously add/read to/from the field at the specified location, support for antennae like structures may be added in later versions with support for a J field
        sourceindex = (2, sourceloc) # source index follows the form (E component to stimulate, [location of source])

        @time field_pml = RCSGenerate(space, time, f₀, nPML; #Try with a PML
                                        source = (Esource, Hsource),
                                        medium = medium,
                                        sourceindex = sourceindex,
                                        detectorindex = detectorindex
                                        )

    return field_pml, field_extended

end
