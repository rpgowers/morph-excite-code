using BSON, PyPlot

function roundint(x)
	if isinteger(x) == true
		return Int(x)
	else
		return x
	end
end

τ_s = 5.0
fig_width = 9 # width of figure
panel_height = 3 # height of a single panel
font_title = 16
font_axis = 12
font_inter = 14
font_legend = 10
linewidth = 2.0

C = get_cmap("viridis")
colors = [C((i-1)*0.225) for i=1:5]
alpha = 0.75