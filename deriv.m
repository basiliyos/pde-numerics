function [ derivu ] = deriv( u,N,k )
% Compute kth spatial derivative of function u spectrally with N Fourier modes
    uhat = fft(u);
    derivvec = [0;(1i.*[1:N/2-1]').^k;0;(1i.*[-N/2+1:-1]').^k];
    derivuhat = derivvec.*uhat;
    derivu = ifft(derivuhat);  
end

