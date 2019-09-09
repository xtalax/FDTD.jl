# FDTD.jl
An implementation of the Finite Difference Time Domain (FDTD) method in 2D and 3D for Electromagnetic Simulation in julia. 

Adapted from [uFDTD (John Schneider)](https://www.eecs.wsu.edu/~schneidj/ufdtd/) for julia, with a CFS-PML.

# Usage:
1. Set up the space and time.
2. Define the medium.
3. Set up the source and detector.
4. Simulate.

This package uses my other package [DiscreteAxis](https://github.com/xtalax/DiscreteAxis.jl) for the space and time at the moment.`LinearAxis` is used for time and is a thin wrapper around a `StepRangeLen`, `Space2D/3D` is essentially a named tuple of `LinearAxis`.


## Setting up the space
There is no need to do this manually, the required time step and spatial step for an accurate and stable simulation is a function of the maximum frequency that will need to be simulated, therefore it is reccomended to use the `setup_spacetime` helper function:
```
fmax = 0.5*10^9 #Maximum frequency present in your source

xmax = ymax = zmax = 0.5
xmin = ymin = zmin = -0.5

nPML = 10 # Number of points for the perfectly matched layer, 0 causes the grid to be truncated with...
          ## perfect electric conductor boundaries.
time_multiplier = 2.0 # How long should the simulation run, in units of the propagation time at the speed of light from...
                      ## one end of the longest axis to the other
space, time = setup_spacetime(fmax, time_multiplier, nPML; 
                                xlims = (xmin, xmax),
                                ylims = (ymin, ymax),
                                zlims = (zmin, zmax))
                      
```
Choose your number of dimensions for the space with the number of `lims` that you pass - just xlims and ylims would result in a 2D grid.

## Initializing your medium
Declare functions for your relative permittivity, relative permeability, electric conductivity and magnetic conductivity:
```
function rel_permittivity(r::AbstractVector{T}) where T = 1 + exp(-sum(r.^2)) # A gaussian shaped relative permittivity
```
It is a requirement that these functions accept a vector with length equal to number of dimensions.
Pass them to the `Medium` constructor with the keyword arguments `ε` (\varepsilon), `μ`(\mu), `σ`(\sigma) and `σm` respectively:

```
Medium(space; ε = rel_permittivity) # defaults to ones everywhere for ε and μ, zeros everywhere for σ and σm
```
There are additional keyword arguments in the case that you would like to rotate your usual medium functions. For 3D there is `azimuth` and `elevation` (specified in radians), in 2D `θ`. 
The centre of rotation is specified by a vector of length equaling the number of dimensions, `rotcentre`
```
Medium(s::Space3D{T}; azimuth = 0.0, elevation = 0.0, rotcentre = [0.0,0.0,0.0])
```

## Setting up the source and detector
You need to define the source E field and the source H field for the whole length of time:

```
Esource = [sin(2π*f*T)*exp(-(c₀*T)^2) for T in time] #Set up the source for the Efield
Hsource = zeros(time.N) #Set the H field to zeros (can be whatever you wish, as long as it is time.N long
```
Then you need to define where it will be applied:
```
sourceindex = (2, [div(x.N,2), div(y.N,2), div(z.N,2)])
```
The above will apply the source on the `y` (2nd) component of the `E` and `H` fields, at the location defined in `sourceindex[2]`.
It is also possible to use ranges and colons to apply the source to a region of the domain. to this end thee is a helper function `interior_range` that will return the part of the grid which is outside of the perfectly matched layer.

```
interior = interior_range(space, nPML) # A function returning the interior of the space
sourceloc = [nPML+2, interior[2], interior[3]] #Applies the source on a whole plane in the interior
sourceindex = (2, sourceloc) #builds the full source index
```
The detector index is defined similarly. The first and second components are which field and which component respectively. Colons may be used to select both fields/all components. The third is again an index, which may contain integers, colons and ranges. It is also possible to use a `CartesianIndex`.
```
detectorindex = (1,2,[interior[1], div(space.x.N,2), interior[3]])
```
## Running the simulation
The simulation is run like this:
```
field = FDTD_propagate(space, time, f₀, nPML;
                                    source = (Esource, Hsource),
                                    medium = medium,
                                    sourceindex = sourceindex,
                                    detectorindex = detectorindex
                                    )
```                                    
`f₀` is the centre frequency of your source, which is used to optimize the PML.

ProgressMeter is in use to give you an idea of the stopping time for the simulation, this can potentially be a very long time for 3D grids at high frequencies or time multipliers.

# TO DO
- Allow anisotropic media
- Allow non-uniform grids
- Add a `J` field and optional `Jsource` for simulating antenna like structures
- Refine the grid based on the medium to ensure enough points per wavelength even inside dielectric media
- Implement a Transmitted Field Scattered Field (TFSF) boundary
- Implement a near field to far field transformation

# References:

[1] John B. Schneider, Understanding the Finite-Difference Time-Domain Method, www.eecs.wsu.edu/~schneidj/ufdtd, 2010. 

[2] S. D. Gedney, “Scaled CFS-PML: It is More Accurate, More Efficient, and Simple to Implement.  Why Aren’t YOU Using It?,” IEEE Trans. Antennas Propagat., Vol. AP-14, pp. 302-307, 1966.  
