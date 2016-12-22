function [ u_p, i ] = Lpnorm( u, ~, tspan, p, i )
% Compute and plot L^p norm of u as a function of time

% Initialize p-norm array
u_p = zeros(1,length(tspan));

% Calculate L^p norm at each time step
for j = 1:length(tspan)
    u_p(j) = norm(u(j,:),p);
end

% Plot L^p norm
figure(i), plot(tspan,u_p),axis([0 max(tspan) 0 u_p(1)*1.1])

% Increment i
i = i+1;

end

