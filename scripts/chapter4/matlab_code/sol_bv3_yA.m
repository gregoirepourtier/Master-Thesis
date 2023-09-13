%% Solution for molar fraction yA at interface

close all; clear;  clc;


%% Mesh and C-rates
x = linspace(0,1,100);
t = linspace(0,1,500);
m = 0;

Ch_list         = [0.04,1.0,5.0,10.0,25.0,50.0,100.0,200.0]; 
params.DA_tilde = 1.0;
params.omega_A  = 10;
params.gamma_A  = 13;

options = odeset('Events',@pdevents);

solutions = {};
times     = {};

%% Problem definition
tic
for k = 1:length(Ch_list)
    params.Ch       = Ch_list(k);

    eqn = @(x,t,u,dudx) pdesys(x,t,u,dudx,params);
    bc = @(xl,ul,xr,ur,t) pdeSysBc(xl,ul,xr,ur,t,params);
    [sol,tsol,sole,te,ie] = pdepe(m,eqn,@pdeSysIc,bc,x,t,options);
    solutions = [solutions, sol];
    times     = [times    , tsol];
end
toc

%% Plotting
figure()
plot(times{1},solutions{1}(:,end),'DisplayName','C_h = 0.04')
hold on
for i = 2:length(solutions)
    tmp = solutions{i}(:,end);
    plot(times{i},tmp,'DisplayName',strcat('C_h=',num2str(Ch_list(i))))
end
hold off

xlim([-0.01 1.01])
ylim([-0.01 1.01])
legend('show','Location','southeast')



%% Nested functions
function [c,f,s] = pdesys(x,t,u,DuDx,params)
    c = params.Ch/params.DA_tilde;

    deno = (1-u)*((1/params.omega_A)*u + (1-u));
    tmp = 1/deno + params.gamma_A*(16*u^3 - 22*u^2 + (25/3)*u);
    f = tmp*DuDx;

    s = 0;
end
% --------------------------------------------------------------
function u0 = pdeSysIc(x)
    u0 = 0.1e-5;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdeSysBc(xl,ul,xr,ur,t,params)
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
