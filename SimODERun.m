function [] = SimODERun( alpha, q0 )
% Set up and solve self-similar ODE and plot results

    
% Remaining initial data
p0 = 0;
s0 = 0;
u0 = [q0,p0,s0]';

% Xi grid setup
dxi = 0.0001;
ximax = 21*sqrt(2)*q0;
xispan = 0:dxi:ximax;

% Solve ODE and extract each column vector from solution
[xi,u] = ode15s(@SimRHS,xispan,u0,[],alpha);
q = u(:,1);
p = u(:,2);
s = u(:,3);

% Plot results
figure(1), subplot(4,1,1), plot(xi,q), title('q')
subplot(4,1,2), plot(xi,p), title('p')
subplot(4,1,3), plot(xi,s), title('s')
subplot(4,1,4), plot(xi,s./q), title('k')

     

end

