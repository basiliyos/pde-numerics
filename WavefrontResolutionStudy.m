function [ i, stufftoplot ] =  ... 
    WavefrontResolutionStudy(kvec, dt, tmax, u0switch, amp0, width, pow, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, ... 
    relativecutoff, method, i )
% Perform resolution study on wavefront motion

% Get length of tspan vector
tlength = tmax / dt + 1;

% Initialize cell to hold legend titles for plotting
legendvec = cell(length(kvec),1);

% Introduce solution vector
stufftoplot = zeros(length(kvec),tlength);

for l = 1:length(kvec)
    k = kvec(l);
    
    % Set up initial data
    [N, h, xmin, xmax, x, tspan, u0, ~, maximum] = ... 
    Setup(k, dt, tmax, u0switch, amp0, width, pow);
   
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
        WavefrontPosition( u, N, t, x, xmin, xmax, relativecutoff, method, wavenumber, i, ...
        amp0, width, 'no');
    
    % Fill plot vector with wavefront position
    stufftoplot(l,:) = wavefront(:);
    
    % Fill legend vector
    legendvec{l} = num2str(k,'2^{%d} modes');

end

% Plot results
colorvec = ['b','r','g','m','k'];
figure(i), hold on
for l = 1:length(kvec)
    plot(tspan(:),stufftoplot(l,:),colorvec(l),'DisplayName',legendvec{l})
end
legend('show','Location','northwest'),
title('Wavefront Position Resolution Study'),
xlabel('t'),ylabel('Wavefront Position')
hold off
% Increment i
i = i+1;

end

