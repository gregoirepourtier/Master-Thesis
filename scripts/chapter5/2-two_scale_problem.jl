### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 15e08880-f25b-11ed-02b5-17edc8ccc7b1
begin
	using Pkg

	Pkg.activate("Project.toml")

	using Revise, SkeelBerzins, Plots, PlutoUI, LaTeXStrings, GLMakie, CairoMakie
end

# ╔═╡ 36745590-eae9-4694-b905-d6ff56960ee7
PlutoUI.TableOfContents()

# ╔═╡ a0ad10e1-5bd2-4e5e-a6d6-9f4dffca97e2
md"""
# Two-Scale Problem using SkeelBerzins.jl
"""

# ╔═╡ 8d4a0592-34b8-4f89-ac4d-76da642cc52e
md"""
## Problem Definition
"""

# ╔═╡ 0f6f3f7d-f7c9-4dba-9d6a-c678585ff180
md"""
Consider the two-scale problem $(x,r) \in [0,1]\times[0,1]$:

$$\partial_t u(x,t) = \partial_x ( f_u(u,\partial_x u) ) + g(u,v(r=1,t;x)) \quad \text{ with } g = -u + v(r=1,x;t) \quad \text{ and }$$
$$\hspace{6.3cm}f_u=\partial_x u$$
$$\partial_t v(r,t;x) = \partial_r ( f_v(v,\partial_r v) ), \quad \text{ with } f_v = \partial_r v$$
$$\hspace{1cm} \text{with BC } f_v|_{r=1} = -g(u(x,t),v(r=1,t;x))$$
"""

# ╔═╡ 762c45a4-a4b5-48bc-ba71-921e6d70e648
begin # Generating Meshes
	nx     = 50
	mesh_x = LinRange(0,1,nx)
	
	nr     = 50
	mesh_r = LinRange(0,1,nr)
end ;

# ╔═╡ eb83df24-80a8-4e49-ba1d-4aa9f2dbbb04
begin # Equation defined on the Macro-Scale
	function pdefun_macro(x,t,u,dudx)
		c = [1,1,1]
		f = [0.5,0.1,1] .* dudx
		s = [-u[1]+u[2], u[1]-u[2], u[2]-u[3]] # ignored only on marked meshpoints --> given by coupling_macro function
		
		return c,f,s
	end
	function icfun_macro(x)
		u0 = [0,0,0]
		
		return u0
	end
	function bdfun_macro(xl,ul,xr,ur,t)
		pl = [ul[1]-1.0, 0, ul[3]]
		ql = [0, 1, 0]
		pr = [0, ur[2], 0]
		qr = [1,0,1]
	
		return pl,ql,pr,qr
	end
end

# ╔═╡ 02802db0-9674-4cf5-9443-739b99da7f6c
begin # Equation defined on the Micro-Scale
	function pdefun_micro(x,t,u,dudx)
		c = 1
		f = dudx
		s = 0
		
		return c,f,s
	end
	function icfun_micro(x,ν)
		u0 = 0.0 # 2
	end
	function bdfun_micro(xl,ul,xr,ur,t)
		pl = 0 # ignored symmetry condition
		ql = 1 # ignored symmetry condition
		pr = 0 # ignored --> given by coupling_micro function
		qr = 1
	
		return pl,ql,pr,qr
	end
end

# ╔═╡ 54d2ffbf-6064-4248-b170-c41a590f63e8
md"""
## Two-Scale test
"""

# ╔═╡ 313e6c05-ebff-4b51-a3ec-f5b7bb0b02d7
function main_two_scale(Nx,x_mesh,Nr,r_mesh)
	
	tspan  = (0, 2)

	m_x = 0
	m_r = 2

	coupling_macro(x,t,u,v) = [-u[1] + v[end], u[1] - u[2], u[2]-u[3]]
	coupling_micro(x,t,u,v) = -u[1] + v[end]

	domain  = ((0.0, 0.3),(0.3,0.6),(0.6,1.0))
	d1,d2,d3 = domain
	markers_macro = hcat([d1[1] <= x <= d1[2] for x in x_mesh],
				   [d1[1] <= x <= d3[2] for x in x_mesh],
				   [d3[1] <= x <= d3[2] for x in x_mesh])

	params = SkeelBerzins.Params()
	params.solver   = :DiffEq
	params.sparsity = :symb
	# params.markers_macro = markers_macro

	markers = [x_mesh[i]<0.3 || x_mesh[i]>0.6 for i ∈ eachindex(x_mesh)]
	

	pb = pdepe(m_x,pdefun_macro,icfun_macro,bdfun_macro,x_mesh,tspan ; 										params=params, 
							mr=m_r, 
							rmesh=r_mesh, 
							pdefun_micro=pdefun_micro,
							icfun_micro=icfun_micro,
							bdfun_micro=bdfun_micro,
							coupling_macro=coupling_macro,
							coupling_micro=coupling_micro, 
							markers_micro=markers)
			
	problem   = DifferentialEquations.ODEProblem(pb)
	sol = DifferentialEquations.solve(problem,Rosenbrock23())

	return sol, pb, markers, pb.npde
end

# ╔═╡ 9b4989cb-6158-4208-93bf-a0fc54567eff
sol, pb, markers_test, npde_macro = main_two_scale(nx,mesh_x,nr,mesh_r) ; tsol = reshape(sol,pb) ;

# ╔═╡ 1bdb850d-8637-4ed2-be23-72c11ad56512
md"""
Timestep: $(@bind t_plot PlutoUI.Slider(1:length(sol.t),default=1,show_value=true))
"""

# ╔═╡ 4e7ed6ca-a097-4e7d-8b90-b01fd166ab99
res_two_scale_3species = [tsol[i+npde_macro].u[t_plot] for i in eachindex(mesh_x[markers_test])] ;

# ╔═╡ d78e4050-5b5b-450f-9895-f1bbf49f76f8
res_final_3species = reduce(hcat,res_two_scale_3species) ; 

# ╔═╡ ac1d3fe2-ccb5-45eb-b1c1-1b031f991d68
@bind resetLR PlutoUI.Button("Reset left/right camera")

# ╔═╡ b6d76a4c-ef31-4062-8b15-303837c4f948
resetLR ; md"""
Rotation left/right: $(@bind rot PlutoUI.Slider(1.0:0.001:2.5,default=1.85,show_value=true))
"""

# ╔═╡ f12d4444-6904-48b6-8620-7c3059542ded
@bind resetUD PlutoUI.Button("Reset up/down camera")

# ╔═╡ 6b9fd32e-990b-4185-8db4-f6be33a38de8
resetUD ; md"""
Rotation up/down: $(@bind elev PlutoUI.Slider(-0.5:0.0001:1.0,default=1/8,show_value=true))
"""

# ╔═╡ 98ed5a83-efa6-40ad-ad52-1b5c50746f3f
let
	idx = 15
	fig = Figure()
	
	CairoMakie.surface(fig[1,1],mesh_r, mesh_x[markers_test][1:idx], res_final_3species[1:end,1:idx], axis=(type=Axis3,xlabel="Micro Scale r", ylabel="Macro Scale x",zlabel="v(r,t,x)",zlabeloffset=50,title="t=$(sol.t[t_plot])",limits=(nothing,nothing,(0,1)),azimuth = rot * pi,elevation=elev*pi),colormap=:jet,colorrange=(0,1))
	
	CairoMakie.surface!(fig[1,1],mesh_r, mesh_x[markers_test][idx+1:end], res_final_3species[1:end,idx+1:end], colormap=:jet,colorrange=(0,1))
	fig[1,2] = Colorbar(fig; vertical=true, width=20,colormap=:jet)

	
	# save("v_specie.svg", fig, pt_per_unit = 2) 
	fig
end

# ╔═╡ 5353feef-11c1-4410-804e-da30cd78899f
let
	f = Figure()
	ax = Axis(f[1, 1], title="t=$(tsol[1].t[t_plot])", xlabel="x", ylabel="u(x,t)", xlabelsize=21.0f0, ylabelsize=21.0f0)

	CairoMakie.xlims!(-0.01,1.01)
	CairoMakie.ylims!(-0.01,1.01)
	
	CairoMakie.lines!(ax,mesh_x,tsol[1].u[t_plot],label=L"u_1")
	if npde_macro != 1
		for i = 2:npde_macro
			CairoMakie.lines!(ax,mesh_x,tsol[i].u[t_plot],label=L"u_{%$(i)}")
		end
	end
	axislegend(ax,position=:rt)
	
	# save("u_species.pdf", f, pt_per_unit = 2)
	f
end

# ╔═╡ b41fffc2-fd3b-479a-8da1-96f8d56621d2
size_macro = length(mesh_x[markers_test]) ;

# ╔═╡ 1d955324-beb1-4ec8-b52c-8b2c6045afc2
let
	cpt = 1
	Plots.plot(mesh_r,tsol[npde_macro+1].u[t_plot],label=L"v_{1}")
	for i=npde_macro+2:size_macro+npde_macro
		cpt += 1
		Plots.plot!(mesh_r,tsol[i].u[t_plot],label=L"v_{%$(cpt)}")
	end
	Plots.xlims!(-0.05, 1.05)
	Plots.ylims!(-0.05, 1.1)
	Plots.plot!(legend=Symbol(:outer, :topleft),size=(500,300))
end

# ╔═╡ Cell order:
# ╠═15e08880-f25b-11ed-02b5-17edc8ccc7b1
# ╟─36745590-eae9-4694-b905-d6ff56960ee7
# ╟─a0ad10e1-5bd2-4e5e-a6d6-9f4dffca97e2
# ╟─8d4a0592-34b8-4f89-ac4d-76da642cc52e
# ╟─0f6f3f7d-f7c9-4dba-9d6a-c678585ff180
# ╠═762c45a4-a4b5-48bc-ba71-921e6d70e648
# ╠═eb83df24-80a8-4e49-ba1d-4aa9f2dbbb04
# ╠═02802db0-9674-4cf5-9443-739b99da7f6c
# ╟─54d2ffbf-6064-4248-b170-c41a590f63e8
# ╠═313e6c05-ebff-4b51-a3ec-f5b7bb0b02d7
# ╠═9b4989cb-6158-4208-93bf-a0fc54567eff
# ╠═4e7ed6ca-a097-4e7d-8b90-b01fd166ab99
# ╠═d78e4050-5b5b-450f-9895-f1bbf49f76f8
# ╟─98ed5a83-efa6-40ad-ad52-1b5c50746f3f
# ╟─1bdb850d-8637-4ed2-be23-72c11ad56512
# ╟─ac1d3fe2-ccb5-45eb-b1c1-1b031f991d68
# ╟─b6d76a4c-ef31-4062-8b15-303837c4f948
# ╟─f12d4444-6904-48b6-8620-7c3059542ded
# ╟─6b9fd32e-990b-4185-8db4-f6be33a38de8
# ╟─5353feef-11c1-4410-804e-da30cd78899f
# ╟─1d955324-beb1-4ec8-b52c-8b2c6045afc2
# ╟─b41fffc2-fd3b-479a-8da1-96f8d56621d2
