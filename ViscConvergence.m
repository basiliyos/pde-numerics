function [ i, stufftoplot ] = ... 
    ViscConvergence( viscvec, k, dt, tmax, u0switch, relativecutoff, method, ... 
    amp0, width, pow, PDE_name, solver, abstol, reltol, deg_nonlinearity, i )
% Check convergence of both types of wavefront approximations as a 
% function of viscosity parameter

% Set up initial data
[N, h, xmin, xmax, x, tspan, u0, ~, maximum] = ...
    Setup(k, dt, tmax, u0switch, amp0, width, pow);

% Initialize cell to hold legend titles for plotting
legendvec = cell(length(viscvec),1);

% Initialize vector of wavefront positions
stufftoplot = zeros(length(viscvec),length(tspan));

for l = 1:length(viscvec)
    
    % Get viscosity from viscosity vector
    visc = viscvec(l);
    
    % Solve PDE with function call
    [ t, u ] = SolvePDE( u0, PDE_name, solver, abstol, reltol, x, ...
        tspan, deg_nonlinearity, N, h, amp0, visc );

    % Calculate wavenumber
    phi = angle(u);
    wavenumber = zeros(size(phi));
    for m = 1:length(t)
        wavenumber(m,:) = deriv(phi(m,:)',N,1)';
    end
    
    % Call function to find wavefront position
    [ wavefront, ~, i ] = ...
        WavefrontPosition( u, N, t, x, xmin, xmax, relativecutoff, method, ..., 
        wavenumber, i, amp0, width, 'no');
    
    % Fill plot vector with wavefront position
    stufftoplot(l,:) = wavefront(:);
        
    % Fill legend vector
    legendvec{l} = strcat(num2str(visc,'Viscosity = %d'),' *h*amp0');
    
end 

% Plot results
colorvec = ['b','r','g','m','c','k'];
figure(i), hold on
for l = 1:length(viscvec)
    plot(tspan,stufftoplot(l,:),colorvec(mod(l-1,6)+1),'DisplayName',legendvec{l}) 
end
legend('show'),
title('Viscosity Parameter Wavefront Position Study'),
xlabel('t'),ylabel('Wavefront Position')

% Increment i
i = i+1;

end

