function [ i ] = WaterfallPlots( u, ux, x, t, i )
% Make waterfall plots of various quantities involving solution u 

% Waterfall plot of Re(|u|^2 u_x) and i increment
figure(i)
waterfall(x,t,real(abs(u.^2).*ux)), title('Re(|u|^2 u_x)'), ... 
    xlabel('x'), ylabel('t')
i = i+1; 

% Waterfall plot of Im(|u|^2 u_x) and i increment
figure(i)
waterfall(x,t,imag(abs(u.^2).*ux)), title('Im(|u|^2 u_x)'), ... 
    xlabel('x'), ylabel('t')
i = i+1;

% Waterfall plot of |(|u|^2 u_x)| and i increment
figure(i)
waterfall(x,t,abs(abs(u.^2.*ux))), title('|(|u|^2u_x)|'), ... 
    xlabel('x'), ylabel('t')
i = i+1;

end

