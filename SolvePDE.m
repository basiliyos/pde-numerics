function [ t, u ] = SolvePDE( u0, PDE_name, solver, abstol, reltol, x, ... 
    tspan, deg_nonlinearity, N, h, amp0, visc )
% Solve choice of PDE given initial data

% Truncate initial data if necessary
if deg_nonlinearity ~= 1
    u0 = truncate(u0,N,deg_nonlinearity);
end

% Set absolute and relative tolerances if necessary
if abstol == 0 && reltol == 0
   options = [];
else
   options = odeset('AbsTol',1^-abstol,'RelTol',1^-reltol); 
end

% Run appropriate choice of ODE solver 
switch solver
    case '45'
        [t, u] = ode45(@RHS,tspan,u0,options,N,visc*h*amp0, ... 
            PDE_name, deg_nonlinearity);
    case '23'
        [t, u] = ode23(@RHS,tspan,u0,options,N,visc*h*amp0, ... 
            PDE_name, deg_nonlinearity);
    case '23tb'
        [t, u] = ode23tb(@RHS,tspan,u0,options,N,visc*h*amp0, ... 
            PDE_name, deg_nonlinearity);
    case '15s'
        [t, u] = ode15s(@RHS,tspan,u0,options,N,visc*h*amp0, ... 
            PDE_name, deg_nonlinearity);
end

end

