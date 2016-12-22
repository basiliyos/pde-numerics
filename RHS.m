function [dudt]  = RHS(t, u, N, epsilon, PDE_name, deg_nonlinearity)
% Compute right-hand side of various PDE



switch PDE_name
    case 'exponential'
        dudt = u;
    case 'characteristic'
        dudt = deriv(u,N,1);
    case 'heat equation'
        dudt = deriv(u,N,2);
    case 'free linear Schrodinger'
        % Parameters, potential function
        hbar = 1;
        m = 1;
        V = 0;  
        % Compute second spectral derivative
        w = deriv(u,N,2);
        dudt = 1i.*hbar.*w./(2.*m) - 1i.*V.*u./hbar; 
    case 'linear Schrodinger with potential'
        % Parameters, potential function
        hbar = 1;
        m = 1; 
        V = cos(2.*x).^2;
        % Compute second spectral derivative
        w = deriv(u,N,2);
        dudt = 1i.*hbar.*w./(2.*m) - 1i.*V.*u./hbar; 
    case 'KdV'
        w = deriv(u,N,3);
        v = deriv(u,N,1);
        vtrunc = truncate(v,N,deg_nonlinearity);
        utrunc = truncate(u,N,deg_nonlinearity);
        dudt = -w - 6.*utrunc.*vtrunc;    
    case 'Burgers'
        v = deriv(u,N,1);
        utrunc = truncate(u,N,deg_nonlinearity);
        vtrunc = truncate(v,N,deg_nonlinearity);
        dudt = utrunc.*vtrunc;
    case 'cubic focusing NLS'
        % Focusing
        mu = 1;
        
        % Compute second spectral derivative
        w = deriv(u,N,2);
        
        % Calculate dudt
        dudt = 1i.*w + 1i.*mu.*2.*u.^2.*conj(u);
        
        % Truncate dudt
        dudt = truncate(dudt,N,deg_nonlinearity);
    case 'cubic defocusing NLS'
        % Defocusing
        mu = -1;
                
        % Compute second spectral derivative
        w = deriv(u,N,2);
        
        % Calculate dudt
        dudt = 1i.*w + 1i.*mu.*2.*u.^2.*conj(u);
        
        % Truncate dudt
        dudt = truncate(dudt,N,deg_nonlinearity);        
    case 'modified DNLS'
        v = deriv(u,N,1);

        derivvec = (1i.*[0:N/2-1 0 -N/2+1:-1]').^2;
        Nmin = floor(N^(1/2)+1);
        Nmax = ceil(N - N^(1/2) + 1);
        derivvec(1:Nmin) = 0;
        derivvec(Nmax:N) = 0;
        viscosity = epsilon*ifft(derivvec.*fft(u));

        w = deriv(abs(u).^2.*v,N,1);
        dispersion = deriv(u,N,-1);
        dudt = 1i.*w - (1/2).*dispersion + viscosity;
    case 'DNLS'
        % Compute first spectral derivative
        v = deriv(u,N,1);

        % % Viscosity?
        % epsilon = 4*h;
        % epsilon = 0;

        % Add second-order numerical viscosity
        viscosity = epsilon*deriv(u,N,2);

        % Compute second spectral derivative
        w = deriv(abs(u).^2.*v,N,1);
        dudt = 1i.*w + viscosity;
        dudt = truncate(dudt,N,deg_nonlinearity);
    case 'DNLS pad'
        upad = 2 * [u(1:N/2); zeros(N+1,1); u((N/2)+2:N)];
        Npad = 2 * N;
        derivvecpad = (1i.*[0:Npad/2-1 0 -Npad/2+1:-1]');
        secondderivvecpad = (1i.*[0:Npad/2-1 0 -Npad/2+1:-1]').^2;
        
        uhatpad = fft(upad);
        vpad = ifft(derivvecpad .* uhatpad);
        
        viscositypad = ifft(epsilon .* secondderivvecpad .* uhatpad);
        
        dudtpad = 1i .* ifft(derivvecpad .* fft(abs(upad).^2 .* vpad)) ...
            + viscositypad;
        
        dudt = 0.5 * [dudtpad(1:N/2);0;dudtpad(3*N/2+2:Npad)];
    case 'DNLS parabolic regularization'
        bilaplacian = deriv(u,N,4);
        parreg = 0.0001;
        
        v = deriv(u,N,1);
        w = deriv(abs(u).^2.*v,N,1);
        dudt = 1i.*w - parreg.*bilaplacian;
end

% Print t to see progress
t

end



