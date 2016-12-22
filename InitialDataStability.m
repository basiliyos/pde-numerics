function [ error, i ] = InitialDataStability(k, dt, tmax, visc, ...
    PDE_name, solver, deg_nonlinearity, abstol, reltol, ...
    firstu0switch, firstamp0, firstwidth, firstpow, ... 
    secondu0switch, secondamp0, secondwidth, secondpow, i)
% Compare resulting solutions of PDE with different initial data by
% plotting L2 norm of difference as a function of t


% Setup and solve PDE for first choice of initial data
[N, h, ~, ~, x, tspan, firstu0, ~, maximum] = ... 
    Setup(k, dt, tmax, firstu0switch, firstamp0, firstwidth, firstpow);

[ t, firstu ] = SolvePDE( firstu0, PDE_name, solver, abstol, reltol, x, ... 
    tspan, deg_nonlinearity, N, h, firstamp0, visc );


% Setup and solve PDE for second choice of initial data
[N, h, ~, ~, x, tspan, secondu0, ~, maximum] = ... 
    Setup(k, dt, tmax, secondu0switch, secondamp0, secondwidth, secondpow);

[ t, secondu ] = SolvePDE( secondu0, PDE_name, solver, abstol, reltol, x, ... 
    tspan, deg_nonlinearity, N, h, secondamp0, visc );

% Compare two initial functions, find Lp norm difference
p = 2;
[error, i] = Lpnorm(secondu - firstu, x, tspan, p, i);


end

