function [ i, stufftoplot, relativeerrorvec ] = ... 
    WavefrontPowerlaw( powvec, dt, tmax, u0switch, amp0, width, k, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, ... 
    relativecutoff, method, i, minindexfrac, maxindexfrac  )
% Study wavefront propagation speed as a function of power in powerlaw
% initial data, find and return error in best linear fit

% Get length of tspan vector
tlength = tmax / dt + 1;

% Convert index ratios to integer indices
minindex = floor(tlength * minindexfrac);
maxindex = floor(tlength * maxindexfrac);

% Initialize plot array
stufftoplot = zeros(length(powvec),tlength);

for l = 1:length(powvec)
    
    pow = powvec(l);
    
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
    
    % Fill plot
    stufftoplot(l,:) = wavefront(:);
    
    % Fill legend vector
    legendvec{l} = num2str(pow, 'Initial data power: %.4f');

end

% Plot various wavefront positions
colorvec = ['b','r','g','m','c','y','k'];
figure(i), hold on
for l = 1:length(powvec)
    plot(tspan(minindex:maxindex),stufftoplot(l,minindex:maxindex), ... 
        colorvec(l),'DisplayName',legendvec{l})
end
legend('show','Location','northwest'),
title('Powerlaw Initial Data Wavefront Behavior'),
xlabel('t'),ylabel('Wavefront Position')
hold off

% Increment i
i = i+1;

% Initiate arrays for linear fit coefficients and function, error vector
polyfitcoefficients = zeros(length(powvec),2);
amilinear = zeros(length(powvec),maxindex-minindex+1);
relativeerrorvec = zeros(1,length(powvec));

% Obtain and evaluate linear fit, compute error
for j = 1:length(powvec)
    polyfitcoefficients(j,:) = polyfit(tspan(minindex:maxindex),stufftoplot(j,minindex:maxindex),1);
    amilinear(j,:) = polyval(polyfitcoefficients(j,:)',tspan(minindex:maxindex));
    relativeerrorvec(j) = norm(stufftoplot(j,minindex:maxindex) - amilinear(j,:)) ... 
        /norm(stufftoplot(j,minindex:maxindex));
end

end

