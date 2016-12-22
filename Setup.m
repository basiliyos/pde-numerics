function [N, h, xmin, xmax, x, tspan, u0, nplots, maximum] = ... 
    Setup(k, dt, tmax, u0switch, amp0, width, pow)
% Basic setup of spatial and temporal vectors and initial data for
% 2*pi-periodic solutions

% Set up spatial grid, define midpoint, width, and initial max value
N = 2 ^ k;
xmin = 0;
xmax = 2 * pi;
h = (xmax - xmin) / N;
x = h*(1:N)';
midpoint = (xmax-xmin)/2;
maximum = amp0;

% Choose u0 based on input string
switch u0switch
    % Note: parabola is a special case of power with pow = 1
    case 'parabola'
        parabola = (width)^2-(midpoint-x).^2;
        u0 = amp0*max(0,parabola);
        maximum = width^2*amp0;
    case 'sim'
        alpha = -1/4;
        u0 = SelfSimInitialData(alpha, amp0, x, N);
        maximum = u0(N/2);
    case 'powerlaw'
        widthalt = width/(xmax-xmin);
        u0 = amp0*max(real(((width)^2-(x-midpoint).^2).^(pow)),0);
        u0(1:N*widthalt) = 0;
        u0(N*(1-(widthalt))+1:N) = 0;
        maximum = u0(N/2);   
    case 'cos'
        u0 = amp0*cos(x).^2;
        u0(1:N/4) = 0;
        u0(3*N/4+1:N) = 0;
    case 'sin exp'
        u0 = amp0*exp(-sin(x));        
    case 'abs cos'
        u0 = amp0*abs(cos(x+midpoint))+1;
        maximum = amp0 + 1;
    case 'sech'
        u0 = amp0*sech(x-midpoint).^2;
    case 'mollifier'
        u0(1:N/4) = 0.;
        a = 0.5 * midpoint;
        b = 1.5 * midpoint;
        xpoly = (x - (a)).*((b)-x);
        u0((N/4)+1:(3*N/4)-1) = exp(-1./xpoly((N/4)+1:(3*N/4)-1));
        u0(3*N/4:N) = 0.;
        u0 = amp0*u0;
    
end

% Set up tspan vector and define number of plots
tmin = 0;
tspan = tmin:dt:tmax;
nplots = round(tmax / dt);

end

