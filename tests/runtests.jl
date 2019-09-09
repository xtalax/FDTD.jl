pwd() ∉ LOAD_PATH && push!(LOAD_PATH, pwd())
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src")
"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\3DFDTD" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\src\\3DFDTD")

"C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests" ∉ LOAD_PATH && push!(LOAD_PATH, "C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
#"/home/zander/Code/julia/RCS-Invert/src" ∉ LOAD_PATH && push!(LOAD_PATH, "/home/zander/Code/julia/RCS-Invert/src")
using Test
using Plots
using Yee

dB(x::Number) = 20*log10(abs(x))
include("C:\\Users\\Jones\\Desktop\\home\\Code\\julia\\RCS-Invert\\tests\\fdtd_test.jl")

detectorindex = (1, 2, [:, :, :]) # detectorindex follows the form (Field, which component to read, [location of detector])
typeof(detectorindex)
field, field_extended = PML_Test(detectorindex)

plotly()
heatmap(dB.(field[100][1][2][:,:,div(54,2)] .-
        field_extended[100][1][2][div(107,2)-div(54,2)+1:div(107,2)+div(54,2),div(107,2)-div(54,2)+1:div(107,2)+div(54,2),div(107,2)]))

heatmap(dB.(field_extended[100][1][2][:,div(107,2),:]))

axislen = size(FDTD_out[1][1][1], 1)
timelen = size(FDTD_out, 1)

lightcone = zeros(axislen, timelen-50)
for i in 1:timelen-50
    lightcone[:,i] .= FDTD_out[i][1][2][div(54,2),div(54,2),:]
end
heatmap(dB.(lightcone))
#plot((wierd))
#@test !any(isnan.(FDTD_out))
#@test !any(isinf.(FDTD_out))
#=
pmllightcone = zeros(size(FDTD_out[1].E.min[1][2],1), timelen)
for i in 1:timelen
    pmllightcone[:,i] = FDTD_out[i].E.max[1][3][:,div(54,2),div(54,2)]
end
heatmap(dB.(pmllightcone))
=#
