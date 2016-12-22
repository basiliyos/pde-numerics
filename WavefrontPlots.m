function [ i ] =  ... 
    WavefrontPlots( wavefrontposition, wavefrontindex, u, ux, t, i, D )
% Plot wavefront-dependent quantities

% Evaluate u and ux at wavefront
uatfront = zeros(length(t));
uxatfront = zeros(length(t));
for j = 1:length(t)
   uatfront(j) = u(j,wavefrontindex(j));
   uxatfront(j) = ux(j,wavefrontindex(j));
end

% Plot Re(|u|^2 u_x) at wavefront, increment i
figure(i),
plot(t,real((abs(uatfront).^2).*uxatfront)), 
title('Re(|u|^2 u_x) at wavefront')
i = i+1;

% Plot Im(|u|^2 u_x) at wavefront, increment i
figure(i),
plot(t,imag((abs(uatfront).^2).*uxatfront)), 
title('Im(|u|^2 u_x) at wavefront')
i = i+1;

% Plot |(|u|^2 u_x)| at wavefront, increment i
figure(i)
plot(t,abs(uxatfront)), 
title('|(|u|^2 u_x)| at wavefront')
i = i+1;

% Compute and plot wavefront speed, increment i
wavefrontspeed = timederiv(t,wavefrontposition');
figure(i)
plot(t,wavefrontspeed)
title('Wavefront speed'), xlabel('t'), ylabel('Speed of wavefront')



end

