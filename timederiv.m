function [ dudt ] = timederiv( tspan, u )
% Calculate derivative of function u with respect to variable t using
% Toeplitz matrix

% Find dt
dt = tspan(2) - tspan(1);

% Time derivative matrix
D = zeros(length(tspan),length(tspan));
for j = 1:length(tspan)-1
   D(j,j+1) = 1;
   D(j+1,j) = -1;
end
D(1,length(tspan)) = -1;
D(length(tspan),1) = 1;

% Normalization
D = ( 1/(2*dt) )*D;

% Calculate time derivative
dudt = D * u;


end

