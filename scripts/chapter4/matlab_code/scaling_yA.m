function scaling_yA(powmax)

    N = (2).^(5:powmax) ;

    t_pdepe = zeros(1,length(N));

    params.DA_tilde = 1.0;
    params.omega_A  = 10;
    params.gamma_A  = 13;
    params.Ch       = 200.0;
    
    for i=1:(powmax-4)
		n     = N(i) ;
        xmesh = linspace(0,1,n) ;
        tspan = linspace(0,1,100);

        f = @() pde_bv3_yA(xmesh,tspan,params);
        t_pdepe(i) = timeit(f);
    end

    loglog(N,t_pdepe,"o-")
    legend('pdepe','Location','northwest')

    fileID = fopen('time_scaling_yA.txt','w');
    fprintf(fileID,'mesh size \t\t time \n');
    fprintf(fileID,'%6.0f ; %15.10f \n',[N;t_pdepe]);
    fclose(fileID);
end
