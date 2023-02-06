using Pkg
Pkg.activate("$(homedir())/.julia/dev/NeuronBifurcate")

using NeuronBifurcate
args = MLS_Param()
println(sn(args))
