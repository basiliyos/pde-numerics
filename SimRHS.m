function [ dudt ] = SimRHS( xi,u,alpha )
% Solve RHS of similarity solution ODE system 
q = u(1);
p = u(2);
s = u(3);

dqdt = p;
dpdt = 2/q*(s.^2-(alpha+0.5)*xi*s-0.25*p.^2);
dsdt = -alpha-(p./q).*(s-0.5*(alpha+0.5)*xi);

dudt = [dqdt,dpdt,dsdt]';

xi
end

