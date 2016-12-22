function [ utrunc ] = truncate(u,N,k)
% % Eliminate Fourier modes of u for a nonlinearity of degree k

% p = portion of modes to be removed
p = 2/(k+1);

% Range of modes to be removed
Nmin = floor((p*N)/2) + 2;
Nmax = ceil(N-((p*N)/2));

% Truncate modes
uhat = fft(u);
uhat(Nmin:Nmax) = 0;
utrunc = ifft(uhat);

end



