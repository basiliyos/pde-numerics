function [ i, plotwithme ] = AmpWidth( k, xmin, xmax, tmax, dt, pow, visc, ... 
    PDE_name, deg_nonlinearity, solver, abstol, reltol, u0switch, ... 
    relativecutoff, method, ampvec, widthvec, i)
% Compare wavefront position using either spectral or spatial method for
% different initial amplitudes and widths

% Define midpoint
xlength = xmax - xmin;
midpoint = xmin + (xlength / 2);


% Initialize cell to hold legend titles for plotting
legendvec = cell(length(ampvec)*length(widthvec),1);

% Form standardized tspan for plotting purposes
tmin = 0;
tspan = tmin:dt:tmax;


% Number of time steps, initialize plot vector
n = tmax/dt+1;
plotwithme = zeros(length(ampvec)*length(widthvec), n);

for l = 1:length(ampvec)
    for j = 1:length(widthvec)
        % Define width and amplitude from vector values
        amp0 = ampvec(l);
        widthsac = widthvec(j);
        width = midpoint*widthsac;
        
        % Rescale time as necessary
        scale = (amp0./(2*widthsac)).^2;
        dtrescaled = dt/scale;
        tmaxrescaled = tmax/scale;
        tplotrescaled = dtrescaled;
        tspanrescaled = tmin:dtrescaled:tmaxrescaled;
        
        % Set up initial data
        [N, h, xmin, xmax, x, ~, u0, ~, maximum] = ... 
            Setup(k, dtrescaled, tmaxrescaled, u0switch, amp0, width, pow);
        
        % Modify timestep based on tmaximum and tplot
        plotgap = round(tplotrescaled/dtrescaled);
        dtrescaled = tplotrescaled/plotgap;
        nplots = round(tmaxrescaled/tplotrescaled);

        % Solve PDE with function call
        [ t, u ] = SolvePDE( u0, PDE_name, solver, abstol, reltol, x, ... 
        tspanrescaled, deg_nonlinearity, N, h, amp0, visc );

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
        plotwithme((l-1)*length(widthvec)+j,:) ... 
            = (wavefront(:)-(pi+width))/amp0;
        
        % Fill legend vector
        legendvec{(l-1)*length(widthvec)+j} = strcat(num2str(amp0,'amplitude=%.1f, '),num2str(width,'width=%.2f'));

    end
end

% Plot various wavefront position vectors
colorvec = ['b','r','g','m','c','k'];
figure(i), hold on
for l = 1:(length(widthvec)*length(ampvec))
    plot(tspanrescaled,plotwithme(l,:),colorvec(mod(l-1,6)+1),'DisplayName',legendvec{l}) 
end
legend('show'),
title('Initial Data Amplitude/Width Wavefront Position Study'),
xlabel('t'),ylabel('Wavefront Position')

% Increment i
i = i+1;

end

