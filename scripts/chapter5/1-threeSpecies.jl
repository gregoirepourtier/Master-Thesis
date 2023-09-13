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

# ╔═╡ 9f689c38-36c5-11ee-2b64-8fea2955f5ea
begin
	using Pkg

	Pkg.activate("Project.toml")

	using Revise, SkeelBerzins, PlutoUI, Plots, ExtendableGrids, SparseDiffTools, CairoMakie
end

# ╔═╡ b9244ad7-469f-47ca-b7d3-17f43ca7c136
md"""
# Three species in different subregions, 1D
"""

# ╔═╡ 4ed755d0-1a48-48a6-8c00-8c34d130d28d
# Problem Formulation
begin
	function pdefun_dom(x,t,u,dudx,domain)
		eps = [1, 1, 1]
	    k   = [1, 1, 1]
		d1,d2,d3 = domain
	
		c = SVector(1,1,1)
		f = [eps[1]*dudx[1], eps[2]*dudx[2], eps[3]*dudx[3]]
		
		if d1[1] <= x <= d1[2]
			s = SVector(-k[1]*u[1] + 1.0e-4*(3.0-x) , k[1]*u[1], 0)
		elseif d3[1] <= x <= d3[2]
			s = SVector(0, -k[3]*u[2], k[3]*u[2])
		else
			s = SVector(0, 0, 0)
		end
		
		c,f,s
	end
	
	function icfun(x)
		u0 = [0,0,0]
	end
	
	function bdfun(xl,ul,xr,ur,t)
		pl = SVector(0, 0, 0)
		ql = SVector(1, 1, 1)
		pr = SVector(0, 0, ur[3])
		qr = SVector(1, 1, 0)

		pl,ql,pr,qr
	end
end

# ╔═╡ 69716f1e-e28d-45f0-8949-5b5c13c55bdd
function main(N_x)
		
	L = 3
	T = 10

	x_mesh = collect(range(0,L,length=N_x))
	tspan  = (0, T)

	m = 0

	domain  = ((0.0, 1.0),(1.0,2.1),(1.9,3.0))
	d1,d2,d3 = domain
	markers = hcat([d1[1] <= x <= d1[2] for x in x_mesh],
				   [d1[1] <= x <= d3[2] for x in x_mesh],
				   [d3[1] <= x <= d3[2] for x in x_mesh])

	pdefun(x,t,u,dudx) = pdefun_dom(x,t,u,dudx,domain)

	params = SkeelBerzins.Params()
	params.solver = :DiffEq
	params.sparsity = :symb
	params.markers_macro = markers

	pb = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan ; params=params)
	problem   = DifferentialEquations.ODEProblem(pb)
	sol = DifferentialEquations.solve(problem,QNDF2())

	x_mesh,sol,pb
end

# ╔═╡ ed7a85f0-5c1d-4277-b2bb-472620f2d6d2
n = 100

# ╔═╡ d27dbad8-d806-49da-9f77-1f11b2c4ee06
xmesh, sol, pb = main(n) ; tres_sys = reshape(sol,pb) ;

# ╔═╡ 45aadeee-702b-4464-b986-122c5fcf110e
# Define the grid and subgrids
begin
	grid=simplexgrid(xmesh)
	cellmask!(grid, [0.0], [1.0], 1)
	cellmask!(grid, [1.0], [2.1], 2)
	cellmask!(grid, [1.9], [3.0], 3)
	
	subgrid1 = subgrid(grid, [1])
	subgrid2 = subgrid(grid, [1, 2, 3])
	subgrid3 = subgrid(grid, [3])
end

# ╔═╡ 12138fbe-c50b-450f-a489-e2370e262f6b
md"""
Timestep: $(@bind t_plot_system PlutoUI.Slider(1:length(sol.t),default=1,show_value=true))
"""

# ╔═╡ 36dc6952-71c7-43ff-b46c-86bcbecb9ccd
let
	mesh1 = subgrid1.components[Coordinates][:]
	mesh2 = subgrid2.components[Coordinates][:]
	mesh3 = subgrid3.components[Coordinates][:]

	U1 = view(tres_sys.u[t_plot_system][1,:], subgrid1)[1:end]
	U2 = view(tres_sys.u[t_plot_system][2,:], subgrid2)[1:end]
	U3 = view(tres_sys.u[t_plot_system][3,:], subgrid3)[1:end]
	
	f = Figure()
	ax = Axis(f[1, 1],title="t=$(tres_sys.t[t_plot_system])", xlabel="x", ylabel=L"u_{approx}", xlabelsize=21.0f0, ylabelsize=21.0f0)

	CairoMakie.xlims!(-0.01,3.01)
	CairoMakie.ylims!(-5e-5, 1e-3)

	CairoMakie.lines!(ax,mesh1,U1,label=L"u_1")
	CairoMakie.lines!(ax,mesh2,U2,label=L"u_2")
	CairoMakie.lines!(ax,mesh3,U3,label=L"u_3")

	axislegend(ax,position=:rt)

	# save("threeSpeciesSol.pdf", f, pt_per_unit = 2)
	f	
end

# ╔═╡ 3a66f32a-91fd-41ec-aaef-1b42f1e187f9
# SkeelBerzins.assemble!(ones(3*n),ones(3*n),pb,0)
# forwarddiff_color_jacobian!(pb.jac,
 #                   (du,u)->SkeelBerzins.assemble!(du,u,pb,0),
           #                 ones(n))

# ╔═╡ Cell order:
# ╠═9f689c38-36c5-11ee-2b64-8fea2955f5ea
# ╟─b9244ad7-469f-47ca-b7d3-17f43ca7c136
# ╠═45aadeee-702b-4464-b986-122c5fcf110e
# ╠═4ed755d0-1a48-48a6-8c00-8c34d130d28d
# ╟─69716f1e-e28d-45f0-8949-5b5c13c55bdd
# ╠═ed7a85f0-5c1d-4277-b2bb-472620f2d6d2
# ╠═d27dbad8-d806-49da-9f77-1f11b2c4ee06
# ╟─12138fbe-c50b-450f-a489-e2370e262f6b
# ╟─36dc6952-71c7-43ff-b46c-86bcbecb9ccd
# ╟─3a66f32a-91fd-41ec-aaef-1b42f1e187f9
