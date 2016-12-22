function [ u0 ] = SelfSimInitialData( alpha, amp0, x, N )
% Construct initial data for u by solving ODE for solutions of the form
% u = t^{alpha} f(x/t^{alpha+1/2})


% Initial data vector
q0 = amp0;
p0 = 0;
s0 = 0;
v0 = [q0,p0,s0]';

% Xi grid
dxi = 0.0001;
ximax = 21*sqrt(2)*q0;
xispan = 0:dxi:ximax;

% Run ODE solver until breakdown time
try
    [xi,g] = ode23tb(@SimRHS,xispan,v0,[],alpha);
end

% Compute ximaxnew
ximaxnew = max(xi);
q = g(:,1);

% Find closest element to pi + ximaxnew in x
[xmin,indxmin] = min(abs(ximaxnew + pi - x));

% Subtract N/2 from index to so that x = pi corresponds to xi = 0
newxispanlength = indxmin - N/2;

% Create new grid, run similarity ODE 
newxispan = linspace(0,ximaxnew,newxispanlength);

% Run similarity ODE with new xi grid
[newxispan2,gnew] = ode23tb(@SimRHS,newxispan,v0,[],alpha);

% Extract columns of new solution, find derivative of u's phase
qnew = gnew(:,1);
snew = gnew(:,3);
phasederiv = snew ./ qnew;

% Find phase using trapezoidal integration
phase = cumtrapz(newxispan,phasederiv);

% Multiply phase by modulus to get similarity solution
u0segment = qnew .^ (1/2) .* exp(1i .* phase);

% Calculate necessary amount of padding
startlength = max(size(u0segment));
topad = (N - (startlength * 2)) / 2;

% Pad u0simsegment with zeros to form u0
u0 = [zeros(topad,1);wrev(u0segment);u0segment;zeros(topad,1)];

end

