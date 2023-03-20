using BSON, JSON

data_DS = BSON.load("data/MLDS-reconstruct-prcs.bson")
gon = data_DS[:gon]
θin_DS = data_DS[:θin]
Δθ_DS = data_DS[:Δθsim]
ΔV_DS = data_DS[:ΔV]
data_out = Dict("gon"=>gon, "θin"=>θin_DS, "Δθ"=>Δθ_DS, "ΔV"=>ΔV_DS)
json_string = JSON.json(data_out)

open("data/MLDS-reconstruct-prcs.json","w") do f 
  write(f, json_string) 
end