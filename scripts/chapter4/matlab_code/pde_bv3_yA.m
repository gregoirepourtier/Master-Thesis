function u = pde_bv3_yA(xmesh,t,params)
    m   = 0;
   
    options = odeset('Events',@pdevents);

    eqn = @(x,t,u,dudx) pdefun(x,t,u,dudx,params);
    bc = @(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,params);
    
    sol = pdepe(m,eqn,@icfun,bc,xmesh,t,options);
    u   = sol(:,:,1);
end

%% Nested Functions
function [c,f,s] = pdefun(x,t,u,DuDx,params)
    c = params.Ch/params.DA_tilde;

    deno = (1-u)*((1/params.omega_A)*u + (1-u));
    tmp = 1/deno + params.gamma_A*(16*u^3 - 22*u^2 + (25/3)*u);
    f = tmp*DuDx;

    s = 0;
end
% --------------------------------------------------------------
function u0 = icfun(x)
    u0 = 0.1e-5;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t,params)
    pl = 0;
    ql = 1;
    pr = -params.Ch/params.DA_tilde;
    qr = 1;
end
%----------------------------------------------
function [value, isterminal, direction] = pdevents(m,t,xmesh,umesh)
value = umesh-(1 - 1e-5);
isterminal = ones(size(umesh));
direction  = zeros(size(umesh));
end
%----------------------------------------------
