function [i, power] = LpResolutionStudy(kvec, dt, tmax, u0switch, amp0, width, pow, ... 
    PDE_name, solver, abstol, reltol, deg_nonlinearity, visc, i)
% Solve PDE at different resolutions, then create comparison plot of L^p
% norms with respect to time normalized by initial data norm and return
% powerlaw fit of L^p norm in time for each resolution

% Get length of tspan vector
tlength = tmax / dt + 1;

% Instantiate index, N, norm, and power vectors
Nvec = 2.^kvec;
maxnorm = zeros(1,length(Nvec));
power = zeros(1,length(kvec));

% Instantiate Lp norm matrices
Lpnorm = zeros(tlength,length(Nvec));
renormLpnorm = zeros(size(Lpnorm));

for l = 1:length(kvec)
   k = kvec(l);
   [N, h, xmin, xmax, x, tspan, u0, nplots, maximum] = ... 
    Setup(k, dt, tmax, u0switch, amp0, width, pow);

   [t, u] = SolvePDE( u0, PDE_name, solver, abstol, reltol, ...
    x, tspan, deg_nonlinearity, N, h, amp0, visc );

    p = 2;
    % Compute Lp norm at each time and save
    for j = 1:length(t)
        Lpnorm(j,l) = norm(u(j,:),p);
    end

    % Delete first fraction of Lp norm data for polyfit
    fitfraction = 0.8;
    fitcutoff = floor((1-fitfraction).*length(t));
    Lpfit = Lpnorm(fitcutoff+1:length(t),l);
    tspanfit = t(fitcutoff+1:length(t));
    
    % Compute and save powerlaw fit for Lp norm
    coeff = polyfit(log(tspanfit)+1,log(Lpfit),1);
    power(l) = coeff(1);
    
     % Lp norm normalized to initial data
     renormLpnorm(:,l) = Lpnorm(:,l)/norm(u0,p);

end

% Plot powerlaw coefficient solution by resolution, increment i
figure(i),
plot(kvec(:),power(kvec(:)-min(kvec)+1)), 
title('Powerlaw behavior for L^p norm'),xlabel('k (number of modes=2^k)'),
ylabel('Power')
set(gca,'xtick',kvec)
i = i+1;

% Plot L^p norms with respect to time at different resolutions
colorvec = ['b','r','g','m','c','k'];
figure(i), hold on
for l = 1:length((Nvec))
    plot(t,renormLpnorm(:,l),colorvec(l))
    legendInfo{l} = ['2^{' num2str(kvec(l)) '} modes'];
end
% Plot title
title('L^p norm Resolution Study')

% Plot legend
legend(legendInfo)

% Increment i
i = i+1;

end

