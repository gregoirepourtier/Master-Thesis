using SkeelBerzins
using BenchmarkTools
using CairoMakie

function process_times_matlab(path)
	f = open(path,"r")
	all_lines = readlines(f)

	close(f)

	filter!(x->!isempty(x),all_lines)
	
	matlab_times = Float64[]
	
	line_count = 0 
	for line in all_lines
		line = filter(x -> !isspace(x), line)
		line = split(line,';')
		if line_count == 0
			nothing
		else
			push!(matlab_times,parse(Float64, line[2]))
		end
		line_count += 1
	end

	matlab_times
end


function solve_PDE_bv3_yA(xmesh,Ch ; DA_tilde=1.0, omega_A=10, gamma_A=13, step=nothing, time_int=TRBDF2(), stopping_crit=1e-5)
	
    function pdefun(x,t,u,dudx)
        c = Ch/DA_tilde

		denominator = (1-u)*((1/omega_A)*u + (1-u))
		thermodynamic_factor = 1/denominator + gamma_A*(16*u^3 - 22*u^2 + (25/3)*u)
        f = thermodynamic_factor*dudx

        s = 0
		
        c,f,s
    end

    icfun(x) = 0.1e-5

    function bcfun(xl,ul,xr,ur,t)
        pl = 0
        ql = 1
        pr = -Ch/DA_tilde
        qr = 1

        pl,ql,pr,qr
    end

    # Parameter definitions
    m = 0
    params_sol = SkeelBerzins.Params()
    params_sol.solver = :DiffEq

	tspan = (0,1)

	# Define callback function to stop the time integration as soon as y_{A|AE} < 1 - stopping_crit
	condition(u,t,integrator) = u[length(xmesh)] - (1 - stopping_crit)
    affect!(integrator)       = terminate!(integrator)
    cb                        = ContinuousCallback(condition,affect!)

    # Solve PDE using DifferentialEquations.jl
    pb = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan ; params=params_sol)
    problem = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(problem, time_int, callback=cb, save_everystep=false)
	
    sol
end



function plot_scaling(powmax)
	
    N = (2).^collect(5:powmax)

    t_pdepe_diffEq_yA = zeros(length(N))
	times_matlab_yA   = zeros(length(N))

	params 		  = Params()
	params.solver = :DiffEq
	
	for i ∈ 1:(powmax-4)
		n     = N[i]
        xmesh = collect(range(0,1,length=n))

		Ch = 200
		
		
		t_pdepe_diffEq_yA[i] = @belapsed solve_PDE_bv3_yA($xmesh,$Ch)
	end

	path_yA = "matlab_code/time_scaling_yA.txt"
	
	times_matlab_yA = process_times_matlab(path_yA)


    f = CairoMakie.Figure(resolution = (600, 500))
	ax = Axis(f[1, 1],
			  xlabel = "size n",
    		  ylabel = "time/s",
			  xscale = log10,
			  yscale = log10,
			  xminorticksvisible = true,
			  yminorticksvisible = true,
			  xminorgridvisible = false,
        	  xminorticks = IntervalsBetween(10),
			  yminorgridvisible = false,
        	  yminorticks = IntervalsBetween(10))

    CairoMakie.lines!(ax,N,t_pdepe_diffEq_yA,label="SkeelBerzins.jl")
    CairoMakie.scatter!(ax,N,t_pdepe_diffEq_yA)

    CairoMakie.lines!(ax, N[1:length(times_matlab_yA)], times_matlab_yA, label="pdepe")
    CairoMakie.scatter!(ax,N[1:length(times_matlab_yA)],times_matlab_yA)

    axislegend(ax,position=:lt)

    f
end


test_pow_max = 14;

f_plot = plot_scaling(test_pow_max)

# save("plot_benchmark.pdf", f_plot, pt_per_unit = 2)
